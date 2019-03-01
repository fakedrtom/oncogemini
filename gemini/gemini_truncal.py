#!/usr/bin/env python
from __future__ import absolute_import

from . import GeminiQuery
#from . import sql_utils

# Truncal mutations are categorized as being absent 
# in a normal tissue sample, but present in all 
# subsequent tumor samples.
# Here we allow for a maximum allele frequency
# in the normal sample, but all tumor samples
# must have allele frequencies greater than that 
# maximum allowed normal allele frequency.
# If more separation from the normal allele frequency
# is desired, an increase parameter can be 
# specified 

def truncal(parser, args):

    # create a connection to the database that was 
    # passed in as an argument via the command line
    gq = GeminiQuery.GeminiQuery(args.db)

    # define sample search query
    if args.purity:
        query = "select patient_id, name, time, purity from samples"
    else:
        query = "select patient_id, name, time from samples"
    
    # execute the sample search query
    gq.run(query)

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
    if args.maxNorm is None:
        maxNorm = float(0)
    else:
        maxNorm = float(args.maxNorm)
    if args.increase is None:
        increase = float(0)
    else:
        increase = float(args.increase)

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
    
    # iterate again through each sample and save which sample is the normal
    # non-normal, tumor sample names are saved to a list
    # establish which timepoints belong to which samples names
    # this is done for the specified --patient and --samples
    # designate the last and first time points
    gq.run(query)
    normal_samples = []
    tumor_samples = []
    timepoints = {}
    for row in gq:
        if row['patient_id'] == patient and row['name'] in samples:
            if int(row['time']) == 0:
                normal_samples.append(row['name'])
            elif int(row['time']) > 0:
                tumor_samples.append(row['name'])
            if int(row['time']) not in timepoints:
                timepoints[int(row['time'])] = []
            timepoints[int(row['time'])].append(row['name'])
#    endpoint = max(timepoints.keys())
#    startpoint = min(timepoints.keys())
#    times = sorted(timepoints.keys(), reverse=True)
    
    # check arrays to see if samples have been added
    # if arrays are empty there is probably a problem in samples
    # check the ped file being loaded into the db
    if len(normal_samples) == 0 and len(tumor_samples) == 0:
        raise NameError('There are no samples; check the ped file for proper format and loading')
    if len(normal_samples) == 0:
        raise NameError('There are no normal samples; check the ped file for proper format and loading')
    if len(tumor_samples) == 0:
        raise NameError('There are no tumor samples; check the ped file for proper format and loading')

    # create a new connection to the database that includes the genotype columns
    # using the database passed in as an argument via the command line
    gq = GeminiQuery.GeminiQuery(args.db, include_gt_cols=True)
    
    # get from the args the maxNorm value
    #if args.maxNorm is None:
    #    maxNorm = str(0)
    #elif args.maxNorm is not None:
    #    maxNorm = args.maxNorm

    # define the truncal query
    if args.columns is not None:
        # the user only wants to report a subset of the columns
        query = "SELECT " + args.columns + " FROM variants"
    else:
        # report the kitchen sink
        query = "SELECT * FROM variants"
    if args.filter is not None:
        # add any non-genotype column limits to the where clause
        query += " WHERE " + args.filter
    # query = "select chrom, start, end, gt_alt_freqs, gt_types from variants where impact_severity !='LOW' and (max_evi =='A' or max_evi == 'B' or max_rating >= 4)"

    # create gt_filter command using saved sample info
#    if args.purity:
#        filter_cmd = ""
#        for sample in normal_samples:
#            filter_cmd += "gt_alt_freqs." + sample + " <= " + maxNorm + " and "
#        for sample in tumor_samples:
#            if sample == tumor_samples[len(tumor_samples)-1]:
#                filter_cmd += "gt_alt_freqs." + sample + " > " + str(float(maxNorm) + float(increase))
#                continue
#            filter_cmd += "gt_alt_freqs." + sample + " > " + str(float(maxNorm) + float(increase)) + " and "
#    else:
#        filter_cmd = ""
#        for sample in normal_samples:
#            filter_cmd += "gt_alt_freqs." + sample + " <= " + maxNorm + " and "
#        for sample in tumor_samples:
#            if sample == tumor_samples[len(tumor_samples)-1]:
#                filter_cmd += "gt_alt_freqs." + sample + " > " + str(float(maxNorm) + float(increase))
#                continue 
#            filter_cmd += "gt_alt_freqs." + sample + " > " + str(float(maxNorm) + float(increase)) + " and " 
#    gt_filter = filter_cmd

    # execute the truncal query (but don't do anything with the results)"
    gq.run(query)#, gt_filter)

    # get the sample index numbers so we can get sample specific GT info (AFs, DPs, etc.)
    smp2idx = gq.sample_to_idx

    # print header and add the AFs of included samples and the calculated slope
    addHeader = []
    for key in timepoints:
        for s in timepoints[key]:
            if s in samples:
                af = 'alt_AF.' + s
                addHeader.append(af)
                if args.purity:
                    raw = 'raw.alt_AF.' + s
                    addHeader.append(raw)
    print(gq.header) + "\t" + '\t'.join(addHeader)

    # iterate through each row of the truncal results and print
    for row in gq:
        normAFs = []
        tumsAFs = []
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
                    if s in normal_samples:
                        normAFs.append(sampleAF)
                    if s in tumor_samples:
                        tumsAFs.append(sampleAF)
                    sampleDP = row['gt_depths'][smpidx]
                    depths.append(sampleDP)
                    sampleGQ = row['gt_quals'][smpidx]
                    quals.append(sampleGQ)
                    addEnd.append(str(sampleAF))
                    if args.purity:
                        addEnd.append(str(rawAF))

        #check that requirements have been met
        if min(depths) < minDP or min(quals) < minGQ:
            continue
        if len(normAFs) > 0 and max(normAFs) > maxNorm:
            continue
        if any(af <= (maxNorm + increase) for af in tumsAFs):
            continue
#        if len(tumsAFs) > 0 and min(tumsAFs) < maxNorm + increase:
#            continue
        # print results that meet the requirements
        # add selected sample AFs
        print str(row) + "\t" + '\t'.join(addEnd)
