#!/usr/bin/python

import sys
import os
from collections import namedtuple

if len(sys.argv) < 4:
    print("Usage: %s <signals_file> <sex_qc_file> <chromosome_boundaries (BED)>" % sys.argv[0])
    sys.exit(0)

# Read in sample sexes
male_keywords = {"Male","male","MALE","M","m","1"}
female_keywords = {"Female","female","FEMALE","F","f","2"}

sample_sexes = {}
read_header=False
with open(sys.argv[2]) as sex_qc:
    for line in sex_qc:
        tokens = line.rstrip().split("\t")
        sample = tokens[0]
        sex = tokens[len(tokens)-1]
        if not read_header and sample == "SampleID":
            read_header=True
            continue
        read_header=True

        if sex in male_keywords:
            sample_sexes[sample] = "Male"
        elif sex in female_keywords:
            sample_sexes[sample] = "Female"
        else:
            print("Unrecognized sex argument in %s: '%s'" % (sys.argv[2], sex))
            sample_sexes[sample] = "."

# Read in chromosome boundary definitions BED file
Chromosome = namedtuple('Chromosome', 'start end centromere_start centromere_end')
chr_bounds = {}
with open(sys.argv[3]) as chr_boundaries:
    for line in chr_boundaries:
        tokens = line.split('\t')
        chrom = tokens[0]
        cent_start_exon = int(tokens[1])
        cent_end_exon = int(tokens[2])
        start = 1
        end = int(tokens[3])
    
        chr_tuple = Chromosome(start, end, cent_start_exon, cent_end_exon)
        chr_bounds[chrom] = chr_tuple

# Parse aggregated signals
Event = namedtuple('Event', 'tier event start end ab_frac sup1 sup2 sup3')
sample_chrs = {}
with open(sys.argv[1]) as signals:
    for line in signals:
        # Parse tokens from each signal line
        tokens  = line.split("\t")
        sample  = tokens[0]
        chrom   = tokens[1]
        start   = int(tokens[2])
        end     = int(tokens[3])
        event   = tokens[4]
        tier    = int(tokens[5])
        ab_frac = float(tokens[6])
        sup1    = float(tokens[7])
        sup2    = float(tokens[8])
        sup3    = None;
        # Last metrics is float or a sex string. Convert sex to 1/2 (M/F) float encoding
        if event == "Coverage":
            if str(tokens[8]) in male_keywords:
                sup3 = 1.0
            elif str(tokens[8]) in female_keywords:
                sup3 = 2.0
        elif event != "LocalMedAB":
            sup3 = float(tokens[9])

        # Build tuple of event metrics
        #event_tuple = (tier, event, start, end, ab_frac, sup1, sup2, sup3)
        event_tuple = Event(tier, event, start, end, ab_frac, sup1, sup2, sup3)

        # Insert event into sample->chromosome->list of events data structure
        if sample in sample_chrs:
            if chrom in sample_chrs[sample]:
                sample_chrs[sample][chrom].append( event_tuple )
            else:
                sample_chrs[sample][chrom] = [ event_tuple ]
        else:
            sample_chrs[sample] = {}
            sample_chrs[sample][chrom] = [ event_tuple ]

# Parse data structure, select/rank anomalies through integrated call metrics
for s in sample_chrs:
    is_male = 0
    if sample_sexes[s] in male_keywords:
        is_male = 1

    for c in sample_chrs[s]:
        chr_tuple = chr_bounds[c]
        # Number of callable bases on the chromosome for normalizing allele balance chromosome fractions
        chr_bases = chr_tuple.end - chr_tuple.start
        chr_bases_nocentromere = chr_bases - (chr_tuple.centromere_end - chr_tuple.centromere_start)

        # Assess copy variant vs. copy neutral first
        copy_neutral = 1    # 1=copy neutral, 0=gain/loss
        copy_tier = 4       # Tier (confidence level) of copy neutral vs. copy variant (4 = no event detected, <4 = tier of coverage event)
        copy_fraction = 0   # Chromosomal dosage change estimate (0=neutral)
        for e in sample_chrs[s][c]:
            if e.event == "Coverage":
                copy_tier = e.tier
                copy_fraction = e.ab_frac
                copy_neutral = 0
                if c == "X" and copy_tier == 3 and abs(e.ab_frac > 0.05):
                    copy_tier = 2
                
        # Given possible copy-state, interpret allele-balance events
        # ROH events trump hetAB events, process them first
        roh_tier = 4        # Minimum tier for ROH (4 = no event)
        start_roh = -1      # Start/end of all ROH events on this chrom. (there may be multiple)
        end_roh = -1
        roh_frac_numerator = 0
        for e in sample_chrs[s][c]:
            if e.event == "ROH":
                if e.tier < roh_tier:
                    roh_tier = e.tier
                if e.start < start_roh or start_roh < 0:
                    start_roh = e.start
                if e.end > end_roh or end_roh < 0:
                    end_roh = e.end
                roh_frac_numerator += (e.end - e.start + 1)
        roh_fraction = 0
        if roh_tier < 4 and start_roh > 0 and end_roh > 0:
            roh_fraction = float(roh_frac_numerator) / float(chr_bases_nocentromere)

        # Next get het AB events (local and chromosome-side)
        # Also check for overlapping ROH, which nullifies the local het AB events
        chrom_ab_tier = 4
        chrom_ab_fit  = -1
        chrom_ab      = -1
        local_ab_tier = 4
        local_ab_starts = []
        local_ab_ends = []
        local_ab_fracs = []
        local_abs = []
        for e in sample_chrs[s][c]:
            if e.event == "LocalMedAB":
                # Require no overlapping ROH event
                if roh_tier == 4 or not (start_roh <= e.end and end_roh >= e.start):
                    if local_ab_tier > e.tier:
                        local_ab_tier = e.tier
                    local_ab_starts.append(e.start)
                    local_ab_ends.append(e.end)
                    local_ab_fraction = (float(e.end-e.start) / float(chr_bases_nocentromere) * ((0.45 - e.ab_frac) / 0.45))
                    local_ab_fracs.append(local_ab_fraction)
                    local_abs.append(e.ab_frac)
                else:
                    # Ignore the event, but expand ROH boundaries if AB event overlaps and extends further
                    if start_roh > 0 and e.start < start_roh and e.end >= start_roh:
                        roh_fraction += float(start_roh - e.start) / float(chr_bases_nocentromere)
                        start_roh = e.start
                    if end_roh > 0 and e.end > end_roh and e.start <= end_roh:
                        roh_fraction += float(e.end - end_roh) / float(chr_bases_nocentromere)
                        end_roh = e.end
                    #if start_roh > 0 and end_roh > 0:
                    #    roh_fraction = float(end_roh - start_roh) / float(chr_bases_nocentromere)
            elif e.event == "ChromMedAB":
                chrom_ab_tier = e.tier
                chrom_ab_fit = e.sup3
                chrom_ab = e.ab_frac

        # Integrate everything to make a final call
        # Interpret chromosome X in males independently of autosomes
        primary_events = {}
        supporting_events = {}
        combined_tier = 4
        ab_fraction = 0.0
        if c == "X" and is_male:
            if copy_tier < 4:
                primary_events['Coverage'] = (copy_tier, copy_fraction)
            if chrom_ab_tier < 4:
                ab_fraction = 1.0 + ((chrom_ab_fit - chrom_ab) / (chrom_ab_fit - 0.3333))
                primary_events['ChromMedAB'] = (chrom_ab_tier, ab_fraction)

            combined_tier = min(copy_tier, chrom_ab_tier)
            if combined_tier >= 2 and copy_tier < 4 and chrom_ab_tier < 4:
                combined_tier = 1
        elif c == "Y":
            primary_events['Coverage'] = (copy_tier, copy_fraction)
            combined_tier = copy_tier
        else:
            # Autosomal call or female chrX
            # Step 1: Decide if marginal coverage calls are copy neutral or not
            if not copy_neutral:
                if copy_tier == 2 and \
                   (roh_tier > 1 or copy_fraction > 0) and \
                   (roh_fraction < 0.1 or copy_fraction > 0) and \
                   chrom_ab_tier > 2 and \
                   local_ab_tier > 2 and \
                   (c != "X" or (chrom_ab_tier > 3 and local_ab_tier > 3)):
                        copy_neutral = 1

                elif copy_tier == 3 and \
                     (roh_tier > 1 or roh_fraction < 0.1) and \
                     chrom_ab_tier > 1 and \
                     local_ab_tier > 1 and \
                     (c != "X" or (chrom_ab_tier > 2 and local_ab_tier > 2)):
                        copy_neutral = 1

            # Step 2: assess AB metrics given assumption of chromosomal gain/loss
            if not copy_neutral:
                if copy_fraction < 0: # Chromosomal loss
                    combined_tier = copy_tier
                    primary_events['Coverage'] = (copy_tier, copy_fraction)
                    if roh_tier <= 2:
                        ab_fraction -= roh_fraction
                        combined_tier = min(combined_tier, 2)
                        if roh_fraction > 0.15 or abs(roh_fraction - abs(copy_fraction)) < 0.2: # Large ROH, or ROH vs dosage fraction estimates within 20%
                            combined_tier = 1
                            primary_events['ROH'] = (roh_tier, roh_fraction)
                        else:
                            supporting_events['ROH'] = (roh_tier, roh_fraction)
                    if local_ab_tier <= 2: #Note: local_ab_... lists have already been ROH-overlap filtered
                        mosaic = "yes"
                        for i in xrange(len(local_ab_fracs)):
                            ab_fraction -= local_ab_fracs[i]
                            if local_abs[i] < 0.05:
                                mosaic = "uncertain"
                        if local_ab_tier <= chrom_ab_tier:
                            primary_events['LocalMedAB'] = (local_ab_tier, ab_fraction, mosaic)
                        else:
                            supporting_events['LocalMedAB'] = (local_ab_tier, ab_fraction, mosaic)
                        combined_tier = min(combined_tier, local_ab_tier)
                    if chrom_ab_tier <= 2 and chrom_ab >= 0:
                        try:
                            chr_fraction = -1.0*((chrom_ab_fit - chrom_ab) / chrom_ab_fit)
                        except:
                            print("Division by 0 %s %s %s" % (s, c, chrom_ab_fit))
                            continue
                        mosaic = "yes"
                        if chrom_ab < 0.05:
                            mosaic = "uncertain"
                        if chrom_ab_tier < min(roh_tier, local_ab_tier):
                            if chr_fraction < ab_fraction:
                                ab_fraction = chr_fraction
                            if chr_fraction < -0.15 or abs(chr_fraction - abs(copy_fraction)) < 0.2:
                                combined_tier = 1
                            primary_events['ChromMedAB'] = (chrom_ab_tier, chr_fraction, mosaic)
                        else:
                            supporting_events['ChromMedAB'] = (chrom_ab_tier, chr_fraction, mosaic)
                        combined_tier = min(combined_tier, chrom_ab_tier)
                    if combined_tier == 2 and \
                        ( (abs(ab_fraction) > 0.15 and abs(copy_fraction) > 0.15) or \
                          (c == "X" and copy_tier <= 2 and (chrom_ab_tier <= 2 or roh_tier <= 2)) ):
                            combined_tier = 1
                    elif combined_tier == 3:
                        n_tier3_metrics = 0
                        if copy_tier <= 3:      n_tier3_metrics += 1
                        if roh_tier <= 3:       n_tier3_metrics += 1
                        if chrom_ab_tier <= 3:  n_tier3_metrics += 1
                        if local_ab_tier <= 3:  n_tier3_metrics += 1
                        if n_tier3_metrics < 2:
                            combined_tier = 4   # Require at least 2 tier3 or better calls to make a reportable combined tier 3 event

                else: # Chromosomal gain
                    combined_tier = copy_tier
                    primary_events['Coverage'] = (copy_tier, copy_fraction)
                    if roh_tier == 1:
                        supporting_events['ROH'] = (roh_tier, roh_fraction) # Shouldn't be relevant for a gain, but add it to the report still
                    if local_ab_tier <= 2:
                        mosaic = "yes"
                        for i in xrange(len(local_ab_fracs)):
                            # fractions are computed relative to 0.45-0 range (i.e. deletion), recalculate for duplication
                            # = * 0.45 * 1/(0.4500-0.3333) = *3.865 such that ab=0.3333 -> 1
                            ab_fraction += local_ab_fracs[i]*3.865
                            if local_abs[i] > 0.31 and local_abs[i] < 0.34: # This will be a minimum AB over the event, thus biased lower
                                mosaic = "no"
                            elif (local_abs[i] > 0.23 and local_abs[i] < 0.26) or \
                                 (local_abs[i] > 0.29 and local_abs[i] < 0.36):     # Expand band, include possible tetrasomies
                                mosaic = "uncertain"
                        if local_ab_tier <= chrom_ab_tier:
                            primary_events['LocalMedAB'] = (local_ab_tier, ab_fraction, mosaic)
                        else:
                            supporting_events['LocalMedAB'] = (local_ab_tier, ab_fraction, mosaic)
                        combined_tier = min(combined_tier, local_ab_tier)
                    if chrom_ab_tier <= 2 and chrom_ab >= 0:
                        chr_fraction = (chrom_ab_fit - chrom_ab) / (chrom_ab_fit - 0.3333)
                        mosaic = "yes"
                        mosaic_error = abs(chrom_ab - 0.3333)
                        if mosaic_error < 0.01:
                            mosaic = "no"
                        elif mosaic_error < 0.02 or abs(chrom_ab - 0.25) < 0.01:
                            mosaic = "uncertain"
                        if chrom_ab_tier < local_ab_tier:
                            if chr_fraction > ab_fraction:
                                ab_fraction = chr_fraction
                            if chr_fraction > 0.15 or ab_fraction > 0.15 or abs(chr_fraction - abs(copy_fraction)) < 0.2:
                                combined_tier = 1
                            primary_events['ChromMedAB'] = (chrom_ab_tier, chr_fraction, mosaic)
                        else:
                            supporting_events['ChromMedAB'] = (chrom_ab_tier, chr_fraction, mosaic)
                    if combined_tier == 2 and abs(ab_fraction) > 0.15 and abs(copy_fraction) > 0.15:
                        combined_tier = 1
                    elif combined_tier == 3:
                        n_tier3_metrics = 0
                        if copy_tier <= 3:      n_tier3_metrics += 1
                        if roh_tier <= 1:       n_tier3_metrics += 1    # for gains, allow tier 1 ROH events to be counted in the aggregate
                        if chrom_ab_tier <= 3:  n_tier3_metrics += 1
                        if local_ab_tier <= 3:  n_tier3_metrics += 1
                        if n_tier3_metrics < 2:
                            combined_tier = 4   # Require at least 2 tier3 or better calls to make a reportable combined tier 3 event

            else: # Copy neutral
                combined_tier = 4
                if copy_tier <= 3:
                    supporting_events['Coverage'] = (copy_tier, copy_fraction)
                    combined_tier = copy_tier

                if roh_tier <= 2:
                    ab_fraction += roh_fraction
                    primary_events['ROH'] = (roh_tier, roh_fraction)
                    if roh_fraction > 0.25:
                        combined_tier = max(roh_tier, 1)
                    elif roh_fraction > 0.15:
                        roh_tier = 2
                        combined_tier = roh_tier
                    else:
                        roh_tier = 3
                        combined_tier = roh_tier
                if local_ab_tier <= 2: #Note: local_ab_... lists have already been ROH-overlap filtered
                    for i in xrange(len(local_ab_fracs)):
                        ab_fraction += local_ab_fracs[i]
                    if local_ab_tier <= chrom_ab_tier:
                        primary_events['LocalMedAB'] = (local_ab_tier, ab_fraction)
                    else:
                        supporting_events['LocalMedAB'] = (local_ab_tier, ab_fraction)
                    combined_tier = min(combined_tier, local_ab_tier)
                if chrom_ab_tier <= 2 and chrom_ab >= 0:
                    try:
                        chr_fraction = (chrom_ab_fit - chrom_ab) / chrom_ab_fit
                    except:
                        print("Division by 0 %s %s %s" % (s, c, chrom_ab_fit))
                        continue
                    if chrom_ab_tier < min(roh_tier, local_ab_tier):
                        if chr_fraction > ab_fraction:
                            ab_fraction = chr_fraction
                        if chr_fraction > 0.15:
                            combined_tier = 1
                        primary_events['ChromMedAB'] = (chrom_ab_tier, chr_fraction)
                    else:
                        supporting_events['ChromMedAB'] = (chrom_ab_tier, chr_fraction)
                    combined_tier = min(combined_tier, chrom_ab_tier)
                if combined_tier == 2 and ab_fraction > 0.15:
                    combined_tier = 1
                elif combined_tier == 3:
                    n_tier3_metrics = 0
                    if copy_tier <= 3:      n_tier3_metrics += 1
                    if roh_tier <= 3:       n_tier3_metrics += 2 # Reports all tier 3 ROH events
                    if chrom_ab_tier <= 3:  n_tier3_metrics += 1
                    if local_ab_tier <= 3:  n_tier3_metrics += 1
                    if n_tier3_metrics < 2:
                        combined_tier = 4   # Require at least 2 tier3 or better calls to make a reportable combined tier 3 event

        # Step 3: Make final call
        if combined_tier <= 3:
            dosage = "neutral"
            if "Coverage" in primary_events:
                dosage = "loss"
                if copy_fraction > 0:
                    dosage = "gain"
                if copy_tier >= 3:
                    dosage = "%s_uncertain" % dosage
            sex = sample_sexes[s]
            scale = "uncertain"
            mosaic = "uncertain"
            if "ROH" in primary_events:
                if abs(ab_fraction) < 0.85:
                    scale = "partial"
                elif abs(ab_fraction) < 0.95:
                    scale = "whole_uncertain"
                else:
                    scale = "whole"
                if dosage == "gain":
                    scale = "%s_uncertain" % scale
                else:
                    mosaic = "no"
            else:
                if "ChromMedAB" in primary_events:
                    scale = "whole"
                    mosaic = "yes"
                    if "LocalMedAB" in supporting_events:
                        scale = "%s_uncertain" % scale
                    if dosage == "neutral":
                        if len(primary_events['ChromMedAB']) > 2 and \
                            primary_events['ChromMedAB'][2] == "uncertain": # ChrABs < 0.05
                                mosaic = "uncertain"
                    elif abs(copy_fraction) > 0.85 and (abs(ab_fraction) > 0.85 or c == "X"):
                        mosaic = "no"
                elif "LocalMedAB" in primary_events:
                    scale = "partial"
                    mosaic = "uncertain"
                    if "ChromMedAB" in supporting_events:
                        scale = "%s_uncertain" % scale
                    if dosage == "neutral":
                        mosaic = "yes" # Assumes a lack of ROH
                    elif len(primary_events['LocalMedAB']) > 2:
                        mosaic = primary_events['LocalMedAB'][2]
            if dosage == "neutral" and copy_tier < 4:
                dosage = "neutral_uncertain"

            for e in sample_chrs[s][c]:
                if e.event in primary_events and e.tier <= primary_events[e.event][0]:
                    print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tPrimary\t%s\t%s" % (s, c, e.start, e.end, e.event, copy_fraction, ab_fraction, dosage, scale, mosaic, e.ab_frac, primary_events[e.event][1], sex, e.tier, combined_tier)
                elif e.event in supporting_events and e.tier <= supporting_events[e.event][0]:
                    print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tSupporting\t%s\t%s" % (s, c, e.start, e.end, e.event, copy_fraction, ab_fraction, dosage, scale, mosaic, e.ab_frac, supporting_events[e.event][1], sex, e.tier, combined_tier)
                else:
                    print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tTertiary\t%s\t%s" % (s, c, e.start, e.end, e.event, copy_fraction, ab_fraction, dosage, scale, mosaic, e.ab_frac, ".", sex, e.tier, combined_tier)












