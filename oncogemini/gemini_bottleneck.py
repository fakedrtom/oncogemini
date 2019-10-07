#!/usr/bin/env python
from __future__ import absolute_import
from . import GeminiQuery
from . import gemini_utils as utils
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
        cancers = args.cancers.split(',')
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
    
    # create a new connection to the database that includes the genotype columns
    # using the database passed in as an argument via the command line
    gq = GeminiQuery.GeminiQuery(args.db, include_gt_cols=True)
    
    # define the bottleneck query
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

    # execute a new query to process the variants
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
        timeAFs = {}
        depths = []
        quals = []
        y = []
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
                    if s in normal_samples and sampleAF >= 0:
                        normAFs.append(sampleAF)
                    if s in tumor_samples and sampleAF >= 0:
                        tumsAFs.append(sampleAF)
                    if key not in timeAFs:
                        timeAFs[key] = []
                    if sampleAF >= 0 :
                        timeAFs[key].append(sampleAF)
                    if sampleAF >= 0:
                        y.append(sampleAF)
                    sampleDP = row['gt_depths'][smpidx]
                    depths.append(sampleDP)
                    sampleGQ = row['gt_quals'][smpidx]
                    quals.append(sampleGQ)
                    addEnd.append(str(sampleAF))
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
