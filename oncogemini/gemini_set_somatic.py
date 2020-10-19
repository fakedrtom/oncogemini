#!/usr/bin/env python
from __future__ import absolute_import, print_function
from .gemini_constants import *
from . import GeminiQuery
from . import gemini_utils as utils

def tag_somatic_mutations(args):

    gq = GeminiQuery.GeminiQuery(args.db)

    # if purity is invoked, get the purity values
    if args.purity:
        query = "select name, purity from samples"
        purity = {}
        gq.run(query)
        utils.get_purity(gq, purity)
    
    # first query is to establish the patients in db and assign sample names 
    # to patient_ids
    query = "select patient_id, name, time from samples"
    gq.run(query)
    
    patients = []
    names = {}
    utils.get_names(gq,patients,names)
    patients = sorted(list(set(patients)))

    # iterate through each patient in the db
    for patient in patients:
        print("Processing patient " + patient)
        samples = 'All'
        samples = utils.get_samples(patient,names,samples)
        
        # second query is to sort samples as either normal or tumor
        # based on time column in samples table
        # query uses same arguments as first query
        query = "select patient_id, name, time from samples"
        gq.run(query)
        normal_samples = []
        tumor_samples = []
        timepoints = {}
        samples_tps = {}
        utils.sort_samples(gq,normal_samples,tumor_samples,timepoints,samples_tps,patient,samples)

        #check and make sure patient has normal samples
        if len(normal_samples) == 0:
            print('No normal samples for patient ' + patient)
            print('Unable to identify somatics')
            print('Skipping patient ' + patient)
            continue

        # final query will actually go through each of the variants
        # and see if the variant matches supplied somatic criteria
        query = "SELECT variant_id, chrom, start, end, ref, alt, gene, \
                        gts, gt_types, gt_ref_depths, gt_alt_depths, \
                        gt_alt_freqs, gt_depths, gt_quals \
                 FROM variants"
        gq.run(query)
        smp2idx = gq.sample_to_idx

        somatic_counter = 0
        somatic_v_ids = []

        if args.dry_run:
            print('Somatics for patient ' + patient)
            print('\t'.join(['variant_id', 'chrom', 'start','end', 'ref', \
                        'alt', 'gene', 'normal_GTs', 'tumor_GTs', 'normal_DPs', \
                        'tumor_DPs', 'norm_counts', 'tum_counts', 'normal_AFs', \
                        'tumor_AFs']))

        for row in gq:
            normDPs = []
            tumDPs = []
            norm_counts = []
            tum_counts = []
            normAFs = []
            tumAFs = []
            normGTs = []
            tumGTs = []
            depths = []
            quals = []
        
            # build lists of metrics for filtering
            for s in samples:
                smpidx = smp2idx[s]
                sampleDP = row['gt_depths'][smpidx]
                depths.append(sampleDP)
                sampleGQ = row['gt_quals'][smpidx]
                quals.append(sampleGQ)
                sample_count = row['gt_alt_depths'][smpidx]
                if args.purity:
                    sampleAF = utils.purityAF(row['gt_alt_freqs'][smpidx],purity[s])
                else:
                    sampleAF = row['gt_alt_freqs'][smpidx]
                sampleGT = row['gt_types'][smpidx]
                if s in normal_samples:
                    normDPs.append(sampleDP)
                    norm_counts.append(sample_count)
                    normAFs.append(sampleAF)
                    normGTs.append(sampleGT)
                if s in tumor_samples:
                    tumDPs.append(sampleDP)
                    tum_counts.append(sample_count)
                    tumAFs.append(sampleAF)
                    tumGTs.append(sampleGT)
            
            # start filtering
            if min(depths) < args.minDP or min(quals) < args.minGQ:
                continue
            if any(gt != HOM_REF for gt in normGTs):
                continue
            if all(gt == HOM_REF for gt in tumGTs):
                continue
            if min(normDPs) < args.normDP or max(norm_counts) > args.normCount or max(normAFs) > args.normAF:
                continue
            if any(dp >= args.tumDP for dp in tumDPs):
                if any(count >= args.tumCount for count in tum_counts):
                    if any(af >= args.tumAF for af in tumAFs):
                        somatic_counter += 1
                        somatic_v_ids.append((1, row['variant_id']))
                        
                        print('\t'.join(str(s) for s in [row['variant_id'], row['chrom'], row['start'], row['end'], row['ref'], \
                                                        row['alt'], row['gene'], normGTs, tumGTs, normDPs, tumDPs, norm_counts, \
                                                        tum_counts, normAFs, tumAFs]))

        if not args.dry_run:
            # establish a connection to the database
            from . import database
            import sqlalchemy as sql
            import sys

            conn, metadata = database.get_session_metadata(args.db)
            
            is_somatic = 'is_somatic_' + patient

            # alter the database by adding a new is_somatic
            alter_qry = "ALTER TABLE variants ADD " + is_somatic + " INTEGER DEFAULT 0"
            try:
                conn.execute(sql.text(alter_qry))
            except sql.exc.OperationalError:
                sys.stderr.write("WARNING: Column \"("
                             + is_somatic
                             + ")\" already exists in variants table. Overwriting values.\n")

            # reset values so that records don't retain old annotations.
            cursor = conn.bind.connect()
            cursor.execute("UPDATE variants SET " + is_somatic + " = 0")

            # now set the identified mutations to True.
            update_qry = "UPDATE variants SET " + is_somatic + " = 1 "
            update_qry += " WHERE variant_id IN (%s)"
            update_qry %= ",".join(str(x[1]) for x in somatic_v_ids)
            res = conn.execute(update_qry)
            assert res.rowcount == somatic_counter
            print("Identified and set", somatic_counter, "somatic mutations")

            # create an index on is_somatic
            cmd = 'create index var_som_idx on variants(' + is_somatic +')'
            try:
                conn.execute(cmd)
            except sql.exc.OperationalError:
                sys.stderr.write("WARNING: Index \"("
                                 + is_somatic
                                 + ")\" already exists for variants table. Overwriting index.\n")

            # save the results
            conn.commit()
        else:
            print("Would have identified and set", somatic_counter, "somatic mutations")

def set_somatic(parser, args):

    tag_somatic_mutations(args)
