#!/usr/bin/env python
from __future__ import absolute_import

from . import GeminiQuery
from . import gemini_utils as utils
#from . import sql_utils
from scipy import stats
import numpy as np

# Bottleneck mutations are categorized as exhibiting 
# increasing allele frequencies over time in multiple 
# tumor samples.
# Here we allow for a maximum allele frequency in the
# normal sample(s) (default is 0).
# If the slope of the allele frequencies across 
# timpoints meets the required slope
# the variant is returned.
# The slope of the tumor only allele frequencies
# must also be positive (not including normal AFs).
# Minimum endpoint frequencies can be specified.
# Minimum differences between first and last
# timepoints can be specified and required.
# Rather than process all samples, a selection
# of samples can be specified.

def bottleneck(parser, args):

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
    if args.maxNorm is None:
        maxNorm = float(0)
    else:
        maxNorm = float(args.maxNorm)
    if args.minSlope is None:
        minSlope = float(0.05)
    else:
        minSlope = float(args.minSlope)
    if args.minR is None:
        minR = float(0.5)
    else:
        minR = float(args.minR)
    if args.samples is None or args.samples.lower() == 'all':
        samples = 'All'
    else:
        samples = args.samples.split(',')
    if args.minEnd is None:
        minEnd = float(0)
    else:
        minEnd = float(args.minEnd)
    if args.endDiff is None:
        endDiff = float(0)
    else:
        endDiff = float(args.endDiff)
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
#        if cancer_abbrevs == 0:
#            raise NameError('No civic_gene_abbreviations or cgi_gene_abbreviations found in database, cannot use --cancers')
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
    patient = utils.get_names(gq,patients,names,patient)
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
#        if int(row['time']) == 0 and row['patient_id'] == patient and row['name'] in samples:
#            normal_samples.append(row['name'])
#        elif int(row['time']) > 0 and row['patient_id'] == patient and row['name'] in samples:
#            tumor_samples.append(row['name'])
#        if row['patient_id'] == patient:
#            if samples == 'All':
#                if int(row['time']) not in timepoints:
#                    timepoints[int(row['time'])] = []
#                timepoints[int(row['time'])].append(row['name'])
#            else:
#                if row['name'] in samples:
#                    if int(row['time']) not in timepoints:
#                        timepoints[int(row['time'])] = []
#                    timepoints[int(row['time'])].append(row['name'])
#    all_samples = normal_samples + tumor_samples
#    endpoint = max(timepoints.keys())
#    startpoint = min(timepoints.keys())
#    times = sorted(timepoints.keys(), reverse=True)
    
    # check arrays to see if samples have been added
    # if arrays are empty there is probably a problem in samples
    # check the ped file being loaded into the db
#    if len(normal_samples) == 0 and len(tumor_samples) == 0:
#        raise NameError('There are no samples; check the ped file for proper format and loading')
#    if len(normal_samples) == 0:
#        raise NameError('There are no normal samples; check the ped file for proper format and loading')
#    if len(tumor_samples) == 0:
#        raise NameError('There are no tumor samples; check the ped file for proper format and loading')
    
    # create a new connection to the database that includes the genotype columns
    # using the database passed in as an argument via the command line
    gq = GeminiQuery.GeminiQuery(args.db, include_gt_cols=True)
    
    # define the loh query
    query = utils.make_query(args.columns,args.filter,cancers)
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
#    filter_cmd = ""
#    count = 0
#    while(count < len(times)-1):
#        samplesA = timepoints[times[count]]
#        samplesB = timepoints[times[count+1]]
#        for a in samplesA:
#            for b in samplesB:
#                filter_cmd += "gt_alt_freqs." + a + " >= gt_alt_freqs." + b + " and "
#        count += 1
#    for sample in normal_samples:
#        filter_cmd += "gt_alt_freqs." + sample + " <= " + maxNorm + " and "
#    endpoint = timepoints[times[0]]
#    for last in endpoint:
#        filter_cmd += "gt_alt_freqs." + last + " > " + str(float(maxNorm) + float(endDiff)) + " and "
#    gt_filter = filter_cmd
#    if gt_filter.endswith(' and '):
#        gt_filter = gt_filter[:-5]

    # execute a new query to process the variants
#    gq.run(query, gt_filter)
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
            if s in samples:
                af = 'alt_AF.' + s
                addHeader.append(af)
                if args.purity:
                    raw = 'raw.alt_AF.' + s
                    addHeader.append(raw)
    addHeader.append('slope')
    addHeader.append('intercept')
    addHeader.append('r_value')
    print '\t'.join(addHeader)

    # iterate through each row of the query results
    # make sure that all args parameters are being met
    for row in gq:
        output = []
        out = str(row).split('\t')
        if cancers != 'none':
            output.extend(out[:len(out)-2])
        else:
            output.extend(out)
        normAFs = []
        tumsAFs = []
#        endAFs = []
#        startAFs = []
        timeAFs = {}
        depths = []
        quals = []
#        count = 0
#        x = []
        y = []
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
                    if s in normal_samples and sampleAF >= 0:
                        normAFs.append(sampleAF)
                    if s in tumor_samples and sampleAF >= 0:
                        tumsAFs.append(sampleAF)
                    if key not in timeAFs:
                        timeAFs[key] = []
                    if sampleAF >= 0 :
                        timeAFs[key].append(sampleAF)
#                    if key == endpoint and sampleAF > 0:
#                        endAFs.append(sampleAF)
#                    if key == startpoint and sampleAF > 0:
#                        startAFs.append(sampleAF)
#                    x.append(count)
                    if sampleAF >= 0:
                        y.append(sampleAF)
                    sampleDP = row['gt_depths'][smpidx]
                    depths.append(sampleDP)
                    sampleGQ = row['gt_quals'][smpidx]
                    quals.append(sampleGQ)
                    addEnd.append(str(sampleAF))
#                    count += 1
                    if args.purity:
                        addEnd.append(str(rawAF))
        endAFs = []
        for i in sorted(timeAFs.keys(), reverse=True):
            if len(timeAFs[i]) > 0:
                endAFs = timeAFs[i]
                break
        startAFs = []
        for i in sorted(timeAFs.keys()):
            if len(timeAFs[i]) > 0:
                startAFs = timeAFs[i]
                break

        #check that requirements have been met
        if min(depths) < minDP or min(quals) < minGQ:
            continue
        if len(normAFs) > 0 and max(normAFs) > maxNorm:
            continue
        if len(y) <= 1:
            continue
        x = range(0,len(y))
        slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
        # spearman is causing some sort of double scalar, division problem when --samples reduces the number of samples
#        spear, spear_p = stats.spearmanr(x,y, nan_policy="omit")
        addEnd.append(str(slope))
        addEnd.append(str(intercept))
        addEnd.append(str(r_value))
        if slope < minSlope or r_value < minR:
            continue
        if min(endAFs) < minEnd:
            continue
        if min(endAFs) - max(startAFs) < endDiff:
            continue

        # check that the slope of the non-normal sample allele frequencies
        # is positive
        if len(tumsAFs) > 1:
            tumsx = range(len(tumsAFs))
            slope, intercept, r_value, p_value, std_err = stats.linregress(tumsx,tumsAFs)
            if slope < 0:
                continue

        # print results that meet the requirements
        # if args.cancer has been used, filter results to cancer matches
        # add selected sample AFs and slope to the end of the line 
        output.extend(addEnd)
        if cancers != 'none':
            abbrevs = str(row['civic_gene_abbreviations']).split(',')  + str(row['cgi_gene_abbreviations']).split(',')
            include = 0
            for c in cancers:
                if c in abbrevs:
                    include += 1
            if include > 0:
                print '\t'.join(output)
        else:
            print '\t'.join(output)
