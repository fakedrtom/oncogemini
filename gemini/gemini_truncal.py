#!/usr/bin/env python
from __future__ import absolute_import

from . import GeminiQuery
from . import gemini_utils as utils
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
    if args.maxNorm is None:
        maxNorm = float(0)
    else:
        maxNorm = float(args.maxNorm)
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
#        cancer_abbrevs = 0
#        for row in gq:
#            fields = str(row).rstrip('\n').split('\t')
#            if fields[1] == 'civic_gene_abbreviations':
#                cancer_abbrevs += 1
#            if fields[1] == 'cgi_gene_abbreviations':
#                cancer_abbrevs += 1
 #       if cancer_abbrevs == 0:
 #           raise NameError('No civic_gene_abbreviations or cgi_gene_abbreviations found in database, cannot use --cancers')
        cancers = args.cancers.split(',')
    if args.purity:
        query = "select name, purity from samples"
        purity = {}
        gq.run(query)
        utils.get_purity(gq, purity)
#    else:
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
#    purity = {}
    utils.get_names(gq,patients,names)
    patient = utils.get_patient(patient,patients)
    if args.somatic_only:
        is_somatic = 'is_somatic_' + patient
    samples = utils.get_samples(patient,names,samples)
#    for row in gq:
#        patients.append(row['patient_id'])
#        if row['patient_id'] not in names:
#            names[row['patient_id']] = []
#        names[row['patient_id']].append(row['name'])
#        if args.purity:
#            purity[row['name']] = float(row['purity'])
#    if args.patient is None and len(set(patients)) == 1:
#        patient = patients[0]
#    elif args.patient is None and len(set(patients)) > 1:
#        raise NameError('More than 1 patient is present, specify a patient_id with --patient')
#    if patient not in patients:
#        raise NameError('Specified patient is not found, check the ped file for available patient_ids')

    # check that specified samples with --samples are present
    # otherwise all names for given patient from ped will asigned to samples list
#    if samples != 'All':
#        for sample in samples:
#            if sample not in names[patient]:
#                raise NameError('Specified samples, ' + sample + ', is not found')
#    elif samples == 'All':
#        samples = names[patient]
    
    # iterate again through each sample and save which sample is the normal
    # non-normal, tumor sample names are saved to a list
    # establish which timepoints belong to which samples names
    # this is done for the specified --patient and --samples
    # designate the last and first time points
    gq.run(query)
    normal_samples = []
    tumor_samples = []
    timepoints = {}
    samples_tps = {}
    utils.sort_samples(gq,normal_samples,tumor_samples,timepoints,samples_tps,patient,samples)
#    for row in gq:
#        if row['patient_id'] == patient and row['name'] in samples:
#            if int(row['time']) == 0:
#                normal_samples.append(row['name'])
#            elif int(row['time']) > 0:
#                tumor_samples.append(row['name'])
#            if int(row['time']) not in timepoints:
#                timepoints[int(row['time'])] = []
#            timepoints[int(row['time'])].append(row['name'])
#    endpoint = max(timepoints.keys())
#    startpoint = min(timepoints.keys())
#    times = sorted(timepoints.keys(), reverse=True)
    
    # check arrays to see if samples have been added
    # if arrays are empty there is probably a problem in samples
    # check the ped file being loaded into the db
#    if len(normal_samples) == 0 and len(tumor_samples) == 0:
#        raise NameError('There are no samples; check the ped file for proper format and loading')
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
#    if args.columns is not None:
        # the user only wants to report a subset of the columns
#        if cancers == 'none':
#            query = "SELECT " + args.columns + " FROM variants"
#        elif cancers != 'none':
#            query = "SELECT " + args.columns + ",civic_gene_abbreviations,cgi_gene_abbreviations FROM variants"
#    else:
        # report the kitchen sink
#        query = "SELECT * FROM variants"
#    if args.filter is not None:
        # add any non-genotype column limits to the where clause
#        query += " WHERE " + args.filter
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
    header = gq.header.split('\t')
    if cancers != 'none':
        addHeader.extend(header[:len(header)-2])
    else:
        addHeader.extend(header)
    for key in timepoints:
        for s in timepoints[key]:
            if s in samples:
                af = 'alt_AF.' + s
                addHeader.append(af)
                if args.purity:
                    raw = 'raw.alt_AF.' + s
                    addHeader.append(raw)
    print '\t'.join(addHeader)

    # iterate through each row of the truncal results and print
    for row in gq:
        output = []
        out = str(row).split('\t')
        if cancers != 'none':
            output.extend(out[:len(out)-2])
        else:
            output.extend(out)
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
                        sampleAF = utils.purityAF(row['gt_alt_freqs'][smpidx],purity[s])
                        rawAF = row['gt_alt_freqs'][smpidx]
                    else:
                        sampleAF = row['gt_alt_freqs'][smpidx]
#                    if sampleAF > 1:
#                        sampleAF = 1
                    if s in normal_samples and sampleAF >= 0:
                        normAFs.append(sampleAF)
                    if s in tumor_samples and sampleAF >= 0:
                        tumsAFs.append(sampleAF)
                    sampleDP = row['gt_depths'][smpidx]
                    depths.append(sampleDP)
                    sampleGQ = row['gt_quals'][smpidx]
                    quals.append(sampleGQ)
                    addEnd.append(str(sampleAF))
                    if args.purity:
                        addEnd.append(str(rawAF))

        # if there are now no values in normal or tumor list, skip variant
        if len(normAFs) == 0 or len(tumsAFs) == 0:
            continue

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
