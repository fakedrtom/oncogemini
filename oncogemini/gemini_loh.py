#!/usr/bin/env python
from __future__ import absolute_import
from . import GeminiQuery
from . import gemini_utils as utils
import sys

# LOH mutations are categorized as being heterozygous 
# in a normal tissue sample, but have risen homozygous 
# levels in the tumor samples.
# Here we allow for a maximum and minimum allele
# frequency in the normal sample to set the boundaries
# for the heterozygous calls in the normal. 
# A minimum allele frequency for the tumor samples is
# also used to define homozygous calls in tumor samples.

def loh(parser, args):

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
    if args.samples is None or args.samples.lower() == 'all':
        samples = 'All'
    else:
        samples = args.samples.split(',')
    if args.specific is None:
        somatic = 'none'
    else:
        somatic = args.specific
    if args.cancers is None:
        cancers = 'none'
    else:
        query = "pragma table_info(variants)"
        gq.run(query)
        utils.check_cancer_annotations(gq)
        cancers = args.cancers.split(',')
    if args.maxNorm is None:
        maxNorm = float(0.7)
    else:
        maxNorm = float(args.maxNorm)
    if args.minNorm is None:
        minNorm = float(0.3)
    else:
        minNorm = float(args.minNorm)
    if args.minTumor is None:
        minTumor = float(0.8)
    else:
        minTumor = float(args.minTumor)
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
    samples = utils.get_samples(patient,names,samples)
    if somatic != 'none' and somatic not in samples:
        sys.exit("Error: Specified sample name with --specific is not found, make sure a single sample only is provided and check the sample manifest file for available sample names")

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
    startpoint = min(timepoints.keys())

    # if only sample included with --specific is the first timepoint, --specific won't work
    if somatic != 'none':
        if samples_tps[somatic] == startpoint:
            sys.exit("Error: Specified sample with --specific is the first timepoint, specify a sample that has a preceding sample")

    # check arrays to see if samples have been added
    # if arrays are empty there is probably a problem in samples
    # check the ped file being loaded into the db
    if len(normal_samples) == 0 and somatic == 'none':
        sys.exit("Error: There are no normal samples; check the sample manifest file for proper format and loading")
    if len(tumor_samples) == 0 and somatic == 'none':
        sys.exit("Error: There are no tumor samples; check the sample manifest file for proper format and loading")

    # create a new connection to the database that includes the genotype columns
    # using the database passed in as an argument via the command line
    gq = GeminiQuery.GeminiQuery(args.db, include_gt_cols=True)
    
    # define the loh query
    if args.columns is not None:
        columns = args.columns
        if cancers != 'none':
            columns = args.columns + ",civic_gene_abbreviations,cgi_gene_abbreviations"
    else:
        columns = args.columns
    if args.filter is not None:
        filter = args.filter
    else:
        filter = str(1)
    query = utils.make_query(columns,filter)

    # execute a new query to process the variants
    gq.run(query)#, gt_filter)

    # get the sample index numbers so we can get sample specific GT info (AFs, DPs, etc.)
    smp2idx = gq.sample_to_idx

    # print header and add the AFs of included samples
    addHeader = []
    header = gq.header.split('\t')
    if cancers != 'none':
        addHeader.extend(header[:len(header)-2])
    else:
        addHeader.extend(header)
    if somatic == 'none':
        for key in timepoints:
            for s in timepoints[key]:
                if s in samples:
                    af = 'alt_AF.' + s
                    addHeader.append(af)
                    if args.purity:
                        raw = 'raw.alt_AF.' + s
                        addHeader.append(raw)
    elif somatic != 'none':
        preceding = samples_tps[somatic] - 1
        for s in timepoints[preceding]:
            if s in samples:
                af = 'alt_AF.' + s
                addHeader.append(af)
                if args.purity:
                    raw = 'raw.alt_AF.' + s
                    addHeader.append(raw)
        af = 'alt_AF.' + somatic
        addHeader.append(af)
        if args.purity:
            raw = 'raw.alt_AF.' + somatic
            addHeader.append(raw)
    print '\t'.join(addHeader)

    # iterate through each row of the results and print
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
        if somatic == 'none':
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
                        sampleDP = row['gt_depths'][smpidx]
                        depths.append(sampleDP)
                        sampleGQ = row['gt_quals'][smpidx]
                        quals.append(sampleGQ)
                        addEnd.append(str(sampleAF))
                        if args.purity:
                            addEnd.append(str(rawAF))
        elif somatic != 'none':
            preceding = samples_tps[somatic] - 1
            for s in timepoints[preceding]:
                smpidx = smp2idx[s]
                if args.purity:
                    sampleAF = utils.purityAF(row['gt_alt_freqs'][smpidx],purity[s])
                    rawAF = row['gt_alt_freqs'][smpidx]
                else:
                    sampleAF = row['gt_alt_freqs'][smpidx]
                if sampleAF >= 0:
                    normAFs.append(sampleAF)
                sampleDP = row['gt_depths'][smpidx]
                depths.append(sampleDP)
                sampleGQ = row['gt_quals'][smpidx]
                quals.append(sampleGQ)
                addEnd.append(str(sampleAF))
                if args.purity:
                    addEnd.append(str(rawAF))
            smpidx = smp2idx[somatic]
            if args.purity:
                sampleAF = utils.purityAF(row['gt_alt_freqs'][smpidx],purity[s])
                rawAF = row['gt_alt_freqs'][smpidx]
            else:
                sampleAF = row['gt_alt_freqs'][smpidx]
            if sampleAF >= 0:
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
        if len(normAFs) == 0 or len(tumsAFs) == 0:
            continue
        if any(af < minNorm or af > maxNorm for af in normAFs):
            continue
        if any(af < minTumor for af in tumsAFs):
            continue

        # print results that meet the requirements
        # if args.cancer has been used, filter results to cancer matches
        # add selected sample AFs
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
