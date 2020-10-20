#!/usr/bin/env python
from __future__ import absolute_import, print_function
from .gemini_constants import *
from . import GeminiQuery
from . import gemini_utils as utils

def tag_somatic_mutations(args):

    gq = GeminiQuery.GeminiQuery(args.db)

    # interpret parameters
    if args.minDP:
        minDP = args.minDP
    else:
        minDP = 0
    if args.minGQ:
        minGQ = args.minGQ
    else:
        minGQ = 0
    normFlag = 0
    if args.normAF:
        normAF = args.normAF
        normFlag = 1
    else:
        normAF = 0
    if args.normCount:
        normCount = args.normCount
        normFlag = 1
    else:
        normCount = 0
    if args.normDP:
        normDP = args.normDP
        normFlag = 1
    else:
        normDP = 0
    tumFlag = 0
    if args.tumAF:
        tumAF = args.tumAF
        tumFlag = 1
    else:
        tumAF = 0
    if args.tumCount:
        tumCount = args.tumCount
        tumFlag = 1
    else:
        tumCount = 0
    if args.tumDP:
        tumDP = args.tumDP
        tumFlag = 1
    else:
        tumDP = 0

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
            if min(depths) < minDP or min(quals) < minGQ:
                continue

            # normal samples filtering
            if normFlag == 0:
                if any(gt != HOM_REF for gt in normGTs):
                    continue
            elif normFlag == 1:
                if args.normDP:
                    if min(normDPs) < args.normDP: #or max(norm_counts) > args.normCount or max(normAFs) > args.normAF:
                        continue
                if args.normCount:
                    if max(norm_counts) > args.normCount:
                        continue
                if args.normAF:
                    if max(normAFs) > args.normAF:
                        continue

            # tumor samples filtering
            # if parameters are provided, have to make sure the same sample passes them
            if tumFlag == 0:
                if all(gt == HOM_REF for gt in tumGTs):
                    continue
                if any(gt == HET for gt in tumGTs) or any(gt == HOM_ALT for gt in tumGTs):
                    somatic_counter += 1
                    somatic_v_ids.append((1, row['variant_id']))
                    print('\t'.join(str(s) for s in [row['variant_id'], row['chrom'], row['start'], row['end'], row['ref'], \
                                                     row['alt'], row['gene'], normGTs, tumGTs, normDPs, tumDPs, norm_counts, \
                                                     tum_counts, normAFs, tumAFs]))
            elif tumFlag == 1:
                tum_passed = {}
                if args.tumDP:
                    passed = [i for i, x in enumerate(tumDPs) if x >= args.tumDP]
                    tum_passed['DP'] = passed
                if args.tumCount:
                    passed = [i for i, x in enumerate(tum_counts) if x >= args.tumCount]
                    tum_passed['Count'] = passed
                if args.tumAF:
                    passed = [i for i, x in enumerate(tumAFs) if x >= args.tumAF]
                    tum_passed['AF'] = passed
                if len(tum_passed) == 0:
                    continue
                else:
                    final = []
                    for i in tum_passed:
                        if len(final) == 0:
                            final = final + tum_passed[i]
                        else:
                            final = list(set(final) & set(tum_passed[i]))
                    if len(final) > 0:
                        somatic_counter += 1
                        somatic_v_ids.append((1, row['variant_id']))
                        print('\t'.join(str(s) for s in [row['variant_id'], row['chrom'], row['start'], row['end'], row['ref'], \
                                                        row['alt'], row['gene'], normGTs, tumGTs, normDPs, tumDPs, norm_counts, \
                                                        tum_counts, normAFs, tumAFs]))
#            if any(dp >= args.tumDP for dp in tumDPs):
#                if any(count >= args.tumCount for count in tum_counts):
#                    if any(af >= args.tumAF for af in tumAFs):
#                        somatic_counter += 1
#                        somatic_v_ids.append((1, row['variant_id']))                        
#                        print('\t'.join(str(s) for s in [row['variant_id'], row['chrom'], row['start'], row['end'], row['ref'], \
#                                                        row['alt'], row['gene'], normGTs, tumGTs, normDPs, tumDPs, norm_counts, \
#                                                        tum_counts, normAFs, tumAFs]))

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
