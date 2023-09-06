#!/bin/bash

while getopts "b:a:p:" arg; do
  case $arg in
    b) bed=$OPTARG;;
    a) bam=$OPTARG;;
    p) prefix=$OPTARG;;
  esac
done

: ${bed?"No bed file provided through -b"}
: ${bam?"No bam file provided through -a"}
: ${prefix?"No prefix provided through -p"}

# merge overlapping amplicons
bedtools merge -i $bed > merged.bed
bed=merged.bed

# process bam
	# OUTPUT FILE
	echo -e "sample\treads\tmapped\ton.target\tdetected_targets\tcov_uniformity" > $prefix.mapping.stats.tsv

	# NO. OF READS AND PROPORTION OF MAPPED READS
        # run samtools
        stats=$(samtools flagstat $bam)
        total_reads=$(echo $stats | grep -oP '^\d+')
        mapped=$(echo $stats | grep -oP '(\d+\.\d+)%' | head -n1 | sed 's/%//')

        # ON-TARGET RATE
        # run bedtools intersect
        intersect_file="${prefix}.intersect.bam"
        bedtools intersect -a $bam -b $bed -u -sorted > $intersect_file
        # calculate on-target rate
        reads_on_target=$(samtools view -F4 $intersect_file | cut -f1 | sort | uniq | wc -l)
        total_mapped_reads=$(samtools view -F4 $bam | cut -f1 | sort | uniq | wc -l)
        on_target_rate=$(echo "scale=4; $reads_on_target/$total_mapped_reads" | bc)
        on_target_rate=$(echo "$on_target_rate * 100" | bc)
        on_target_rate=$(printf "%.2f" $on_target_rate)
	# remove intermediate files
	rm $intersect_file

        # DETECTED TARGETS
        # run bedtools
	genome_file="${prefix}.genome"
	samtools view -H ${bam} | grep -P "@SQ\tSN:" | sed 's/@SQ\tSN://' | sed 's/\tLN:/\t/' > $genome_file
        targets_file="${prefix}.targets.bed"
        bedtools intersect -a $bed -b $bam -c -sorted -g $genome_file > ${targets_file}.unfiltered
        awk '$4 !=0' ${targets_file}.unfiltered > $targets_file
        rm ${targets_file}.unfiltered
	rm $genome_file
        # calculate detected targets proportion
        total_targets=$(cat $bed | wc -l)
        detected_targets=$(cat $targets_file | wc -l)
        targets_detected=$(echo "scale=4; $detected_targets/$total_targets" | bc)
        targets_detected=$(echo "$targets_detected * 100" | bc)
        targets_detected=$(printf "%.2f" $targets_detected)
	rm ${targets_file}

	# COVERAGE UNIFORMITY
	# calculate mean read depth of all amplicons in the panel
	mean_read_depth=$(samtools depth -a $bam | awk '{sum+=$3} END {print sum/NR}')
	# calculate coverage uniformity
	threshold_depth=$(echo "scale=4; $mean_read_depth * 0.2" | bc)
	samtools depth -a $bam | awk -v threshold="$threshold_depth" '$3 >= threshold {print $1"\t"$2"\t"$2}' > temp.bed
	amplicons_above_threshold=$(bedtools merge -i temp.bed | bedtools intersect -a $bed -b stdin | wc -l)
	rm temp.bed
	amplicon_count=$(wc -l < $bed)
	# calculate cov uniformity
	cov_uniformity=$(echo "scale=4; $amplicons_above_threshold/$amplicon_count" | bc)
	cov_uniformity=$(echo "$cov_uniformity * 100" | bc)
	cov_uniformity=$(printf "%.2f" $cov_uniformity)

        # append results to the tsv file
        echo -e "${prefix}\t${total_reads}\t${mapped}\t${on_target_rate}\t${targets_detected}\t${cov_uniformity}" >> $prefix.mapping.stats.tsv

	# create a plot with the info
	Rscript PLOT_Flomics_QC_mapping.R $prefix.mapping.stats.tsv


rm $bed # remove the merged bed file and leave input bed intact

#####################
