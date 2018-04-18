#!/usr/bin/env python
from __future__ import absolute_import

from . import GeminiQuery
#from . import sql_utils

# Bottleneck mutations are categorized as exhibiting 
# increasing allele frequencies over time in multiple 
# tumor samples.
# Here we allow for a maximum allele frequency in the
# normal sample(s) (default is 0).
# Tumor samples are ordered by their timepoints and  
# each exhibits a higher allele frequency than its
# predecessor.
# A required amount of increase between timepoints
# can be specified (default is 0).
# Additionally a total end difference between the 
# normal samples and last timepoint can be specified
# (default is 0).

def bottleneck(parser, args):

    # create a connection to the database that was 
    # passed in as an argument via the command line
    gq = GeminiQuery.GeminiQuery(args.db)

    # get paramters from the args for filtering
    if args.patient is not None:
        patient = args.patient
    if args.maxNorm is None:
        maxNorm = str(0)
    else:
        maxNorm = args.maxNorm
    if args.increase is None:
        increase = str(0)
    else:
        increase = args.increase
    if args.endDiff is None:
        endDiff = str(0)
    else:
        endDiff = args.endDiff

    # define sample search query
    query = "select patient_id, name, time from samples"

    # execute the sample search query
    gq.run(query)

    # designating which patient to perform the query on
    # if no patient is specified at the command line
    # and only 1 patient is present in the database
    # that patient will be used
    # also verify that patient is among possible patient_ids
    patients = []
    for row in gq:
        patients.append(row['patient_id'])
    if args.patient is None and len(set(patients)) == 1:
        patient = patients[0]
    elif args.patient is None and len(set(patients)) > 1:
        raise NameError('More than 1 patient is present, specify a patient_id with --patient')
    if patient not in patients:
        raise NameError('Specified patient is not found, check the ped file for available patient_ids')
    
    # iterate again through each sample and save which sample is the normal
    # non-normal sample names are saved to a list
    # establish which timepoint is the last
    gq.run(query)
    normal_samples = []
    other_samples = []
    timepoints = {}
    for row in gq:
        if int(row['time']) == 0 and row['patient_id'] == patient:
            normal_samples.append(row['name'])
        elif int(row['time']) > 0 and row['patient_id'] == patient:
            other_samples.append(row['name'])
        if row['patient_id'] == patient:
            if int(row['time']) not in timepoints:
                timepoints[int(row['time'])] = []
            timepoints[int(row['time'])].append(row['name'])
    endpoint = max(timepoints.keys())
    times = sorted(timepoints.keys(), reverse=True)

    # check arrays to see if samples have been added
    # if arrays are empty there is probably a problem in samples
    # check the ped file being loaded into the db
    if len(normal_samples) == 0:
        raise NameError('There are no normal samples; check the ped file for proper format and loading')
    if len(other_samples) == 0:
        raise NameError('There are no tumor samples; check the ped file for proper format and loading')

    # create a new connection to the database that includes the genotype columns
    # using the database passed in as an argument via the command line
    gq = GeminiQuery.GeminiQuery(args.db, include_gt_cols=True)
    
    # define the loh query
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
    filter_cmd = ""
    count = 0
    while(count < len(times)-1):
        samplesA = timepoints[times[count]]
        samplesB = timepoints[times[count+1]]
        for a in samplesA:
            for b in samplesB:
                filter_cmd += "gt_alt_freqs." + a + " >= gt_alt_freqs." + b + " and "
        count += 1
    for sample in normal_samples:
        filter_cmd += "gt_alt_freqs." + sample + " <= " + maxNorm + " and "
    endpoint = timepoints[times[0]]
    for last in endpoint:
        for sample in normal_samples:
            filter_cmd += "gt_alt_freqs." + last + " > " + "gt_alt_freqs." + sample + " and "
#    for sample in other_samples:
#        if sample == other_samples[len(other_samples)-1]:
#            filter_cmd += "gt_alt_freqs." + sample + " > " + minTumor
#            continue 
#        filter_cmd += "gt_alt_freqs." + sample + " > " + minTumor + " and " 
    gt_filter = filter_cmd
    if gt_filter.endswith(' and '):
        gt_filter = gt_filter[:-5]

    # execute the truncal query (but don't do anything with the results)"
    gq.run(query, gt_filter)

    # iterate through each row of the truncal results and print
    for row in gq:
        print(row)
