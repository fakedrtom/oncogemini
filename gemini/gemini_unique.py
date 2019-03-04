#!/usr/bin/env python
from __future__ import absolute_import

from . import GeminiQuery
#from . import sql_utils

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
        cancer_abbrevs = 0
        for row in gq:
            fields = str(row).rstrip('\n').split('\t')
            if fields[1] == 'civic_gene_abbreviations':
                cancer_abbrevs += 1
            if fields[1] == 'cgi_gene_abbreviations':
                cancer_abbrevs += 1
        if cancer_abbrevs == 0:
            raise NameError('No civic_gene_abbreviations or cgi_gene_abbreviations found in database, cannot use --cancers')
        cancers = args.cancers.split(',')
    if args.specific is None:
        raise NameError('No sample(s) specified with --specific, please provide sample(s)')
    else:
        specific = args.specific.split(',')

    # define sample search query                                                                                                                                                                  
    if args.purity:
        query = "select patient_id, name, time, purity from samples"
    else:
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
    purity = {}
    for row in gq:
        patients.append(row['patient_id'])
        if row['patient_id'] not in names:
            names[row['patient_id']] = []
        names[row['patient_id']].append(row['name'])
        if args.purity:
            purity[row['name']] = float(row['purity'])
    if args.patient is None and len(set(patients)) == 1:
        patient = patients[0]
    elif args.patient is None and len(set(patients)) > 1:
        raise NameError('More than 1 patient is present, specify a patient_id with --patient')
    if patient not in patients:
        raise NameError('Specified patient is not found, check the ped file for available patient_ids')

    # check that specified samples with --samples are present
    # otherwise all names for given patient from ped will asigned to samples list
    if samples != 'All':
        for sample in samples:
            if sample not in names[patient]:
                raise NameError('Specified samples, ' + sample + ', is not found')
    elif samples == 'All':
        samples = names[patient]
    #Also make sure the samples requested with --specific are present
    for s in specific:
        if s not in samples:
            raise NameError('Sample listed with --specific, ' + s + ', is not found, check the ped file for available samples')

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
#    endpoint = max(timepoints.keys())
#    startpoint = min(timepoints.keys())
#    times = sorted(timepoints.keys(), reverse=True)
    
    # check arrays to see if samples have been added
    # if arrays are empty there is probably a problem in samples
    # check the ped file being loaded into the db
    if len(other_samples) == 0 and len(unique_samples) == 0:
        raise NameError('There are no samples; check the ped file for proper format and loading')
    if len(other_samples) == 0 and len(unique_samples) > 0:
        raise NameError('There are no other samples to compare --specific samples to; check the ped file for proper format and loading')
    if len(unique_samples) == 0:
        raise NameError('There are no --specific samples; check the ped file for proper format and loading')

    # create a new connection to the database that includes the genotype columns
    # using the database passed in as an argument via the command line
    gq = GeminiQuery.GeminiQuery(args.db, include_gt_cols=True)
    
    # define the unique query
    if args.columns is not None:
        # the user only wants to report a subset of the columns
        if cancers == 'none':
            query = "SELECT " + args.columns + " FROM variants"
        elif cancers != 'none':
            query = "SELECT " + args.columns + ",civic_gene_abbreviations,cgi_gene_abbreviations FROM variants"
    else:
        # report the kitchen sink
        query = "SELECT * FROM variants"
    if args.filter is not None:
        # add any non-genotype column limits to the where clause
        query += " WHERE " + args.filter

    # execute the unique query (but don't do anything with the results)"
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
                        sampleAF = float(row['gt_alt_freqs'][smpidx]/purity[s])
                        rawAF = row['gt_alt_freqs'][smpidx]
                    else:
                        sampleAF = row['gt_alt_freqs'][smpidx]
                    if sampleAF > 1:
                        sampleAF = 1
                    if s in specific:
                        uniqAFs.append(sampleAF)
                        addEnd.append(str(sampleAF))
                        if args.purity:
                            addEnd.append(str(rawAF))
                    if s not in specific:
                        otherAFs.append(sampleAF)
                    sampleDP = row['gt_depths'][smpidx]
                    depths.append(sampleDP)
                    sampleGQ = row['gt_quals'][smpidx]
                    quals.append(sampleGQ)

        #check that requirements have been met
        if min(depths) < minDP or min(quals) < minGQ:
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
