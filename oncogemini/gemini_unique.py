#!/usr/bin/env python
from __future__ import absolute_import
from . import GeminiQuery
from . import gemini_utils as utils
import sys

#Unique mutations are unique to the individual sample(s)
#Not shared with other samples in the database
#If multiple samples are specified
#Mutations are shared between those samples
#But absent from all others

def unique(parser, args):

    # create a connection to the database that was 
    # passed in as an argument via the command line
    gq = GeminiQuery.GeminiQuery(args.db)

    # get paramters from the args for filtering
    if args.patient is not None:
        patient = args.patient
    else:
        patient = 'none'
    if args.minDP is None:
        minDP = int(-1)
    else:
        minDP = int(args.minDP)
    if args.minGQ is None:
        minGQ = int(-1)
    else:
        minGQ = int(args.minGQ)
    if args.samples is None or args.samples == 'All':
        samples = 'All'
    else:
        samples = args.samples.split(',')
    if args.maxOthers is None:
        maxOthers = float(0)
    else:
        maxOthers = float(args.maxOthers)
    if args.increase is None:
        increase = float(0)
    else:
        increase = float(args.increase)
    if args.cancers is None:
        cancers = 'none'
    else:
        query = "pragma table_info(variants)"
        gq.run(query)
        utils.check_cancer_annotations(gq)
        cancers = args.cancers.split(',')
    if args.specific is None:
        sys.exit("No sample(s) specified with --specific, please provide sample(s)")
    else:
        specific = args.specific.split(',')
    if args.purity:
        query = "select name, purity from samples"
        purity = {}
        gq.run(query)
        utils.get_purity(gq, purity)

    # define sample search query
    query = "select patient_id, name, time from samples"

    # execute the sample search query
    gq.run(query)

    # designating which patient to perform the query on
    # if no patient is specified at the command line
    # and only 1 patient is present in the database
    # that patient will be used
    # also verify that patient is among possible patient_ids
    # sample names are saved to patient specific dict
    patients = []
    names = {}
    utils.get_names(gq,patients,names)
    patient = utils.get_patient(patient,patients)
    if args.somatic_only:
        is_somatic = 'is_somatic_' + patient
    samples = utils.get_samples(patient,names,samples)

    #Make sure the samples requested with --specific are present
    for s in specific:
        if s not in samples:
            sys.exit("Sample listed with --specific, " + s + ", is not found, check the sample manifest file for available samples")

    # iterate again through each sample and save which sample is the normal
    # non-normal, tumor sample names are saved to a list
    # establish which timepoints belong to which samples names
    # this is done for the specified --patient and --samples
    # designate the last and first time points
    gq.run(query)
    other_samples = []
    unique_samples = []
    timepoints = {}
    for row in gq:
        if row['patient_id'] == patient and row['name'] in samples:
            if row['name'] in specific:
                unique_samples.append(row['name'])
            elif row['name'] not in specific:
                other_samples.append(row['name'])
            if int(row['time']) not in timepoints:
                timepoints[int(row['time'])] = []
            timepoints[int(row['time'])].append(row['name'])
    
    # check arrays to see if samples have been added
    # if arrays are empty there is probably a problem in samples
    # check the ped file being loaded into the db
    if len(other_samples) == 0 and len(unique_samples) == 0:
        sys.exit("There are no samples; check the sample manifest file for proper format and loading")
    if len(other_samples) == 0 and len(unique_samples) > 0:
        sys.exit("There are no other samples to compare --specific samples to; check the sample manifest file for proper format and loading")
    if len(unique_samples) == 0:
        sys.exit("There are no --specific samples; check the sample manifest file for proper format and loading")

    # create a new connection to the database that includes the genotype columns
    # using the database passed in as an argument via the command line
    gq = GeminiQuery.GeminiQuery(args.db, include_gt_cols=True)
    
    # define the unique query
    if args.columns is not None:
        columns = args.columns
        if cancers != 'none':
            columns = args.columns + ",civic_gene_abbreviations,cgi_gene_abbreviations"
    else:
        columns = args.columns
    if args.filter is not None:
        filter = args.filter
        if args.somatic_only:
            filter = args.filter + 'and ' + is_somatic + '==1'
    else:
        filter = str(1)
        if args.somatic_only:
            filter = is_somatic + '==1'
    query = utils.make_query(columns,filter)

    # execute the unique query
    gq.run(query)

    # get the sample index numbers so we can get sample specific GT info (AFs, DPs, etc.)
    smp2idx = gq.sample_to_idx

    # print header and add the AFs of included samples and the calculated slope
    addHeader = []
    header = gq.header.split('\t')
    if cancers != 'none':
        addHeader.extend(header[:len(header)-2])
    else:
        addHeader.extend(header)
    for key in timepoints:
        for s in timepoints[key]:
            if s in specific:
                af = 'alt_AF.' + s
                addHeader.append(af)
                if args.purity:
                    raw = 'raw.alt_AF.' + s
                    addHeader.append(raw)
    print '\t'.join(addHeader)

    # iterate through each row of the unique results and print
    for row in gq:
        output = []
        out = str(row).split('\t')
        if cancers != 'none':
            output.extend(out[:len(out)-2])
        else:
            output.extend(out)
        otherAFs = []
        uniqAFs = []
        depths = []
        quals = []
        addEnd = []
        for key in timepoints:
            for s in timepoints[key]:
                if s in samples:
                    smpidx = smp2idx[s]
                    if args.purity:
                        sampleAF = utils.purityAF(row['gt_alt_freqs'][smpidx],purity[s])
                        rawAF = row['gt_alt_freqs'][smpidx]
                    else:
                        sampleAF = row['gt_alt_freqs'][smpidx]
                    if s in specific:
                        if sampleAF >= 0:
                            uniqAFs.append(sampleAF)
                        addEnd.append(str(sampleAF))
                        if args.purity:
                            addEnd.append(str(rawAF))
                    if s not in specific and sampleAF >= 0:
                        otherAFs.append(sampleAF)
                    sampleDP = row['gt_depths'][smpidx]
                    depths.append(sampleDP)
                    sampleGQ = row['gt_quals'][smpidx]
                    quals.append(sampleGQ)

        #check that requirements have been met
        if min(depths) < minDP or min(quals) < minGQ:
            continue
        if len(otherAFs) == 0 or len(uniqAFs) == 0:
            continue
        if len(otherAFs) > 0 and max(otherAFs) > maxOthers:
            continue
        if any(af <= (maxOthers + increase) for af in uniqAFs):
            continue

        # print results that meet the requirements
        # if args.cancer has been used, filter results to cancer matches
        # add selected sample AFs
        output.extend(addEnd)
        if cancers != 'none':
            abbrevs = str(row['civic_gene_abbreviations']).split(',') + str(row['cgi_gene_abbreviations']).split(',')
            include = 0
            for c in cancers:
                if c in abbrevs:
                    include += 1
            if include > 0:
                print '\t'.join(output)
        else:
            print '\t'.join(output)
