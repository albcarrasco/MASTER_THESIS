#!/bin/bash

while getopts "v:p:" arg; do
	case $arg in
	   v) vcf_file=$OPTARG;;
	   p) sample=$OPTARG;;
	esac
done

: ${vcf_file?"No vcf file provided through -v"}
: ${sample?"No prefix provided through -p"}

# process VCF
    
	# remove germline and pon variants
	echo "remove..."
	bcftools filter -e'FILTER~"germline"' $vcf_file | bcftools filter -e'FILTER~"panel_of_normals"' -Oz > $sample.somatic.vcf.gz
	tabix -p vcf $sample.somatic.vcf.gz

	# normalize the vcf to be annotated
	echo "normalize..."
	bcftools norm -m -any $sample.somatic.vcf.gz -Oz > $sample.norm.vcf.gz
	tabix -p vcf $sample.norm.vcf.gz

	# annotate using COSMIC vcfs
	echo "annotate..."
	bcftools annotate -a !{cosmic_vcf} -c ID,INFO $sample.norm.vcf.gz -Oz > $sample.ann.vcf.gz
	tabix -p vcf $sample.ann.vcf.gz

	# remove variants without annotation in COSMIC
	echo "remove..."
	bcftools filter -e'ID="."' $sample.ann.vcf.gz -Oz > $sample.cosmic.vcf.gz
	tabix -p vcf $sample.cosmic.vcf.gz

	# annotate using CIVIC vcfs
	echo "annotate..."
	bcftools annotate -a !{civic_vcf} -c INFO/GN,INFO/VT $sample.cosmic.vcf.gz -Oz > $sample.cosmic.civic.vcf.gz
	tabix -p vcf $sample.cosmic.civic.vcf.gz

	# create the VARIANTS tsv file
	echo "creating variants tsv..."
	echo -e "Gene\tLEGACY_MUTATION_ID\tmolecular_profile\tCDS_mutation\tAA_mutation\t%_cfDNA_or_Amplification\tNumber_of_Cancer_Samples_(COSMIC)" > PRE_VARIANTS.tsv
	bcftools query -f'%INFO/GENE\t%INFO/LEGACY_ID\t%INFO/GN %INFO/VT\t%INFO/CDS\t%INFO/AA\t[%AF]\t%INFO/CNT\n' $sample.cosmic.civic.vcf.gz | awk 'BEGIN{FS=OFS="\t"} {for (i=1; i<=NF; i++) gsub(/_.*/, "", $i)} 1' >> PRE_VARIANTS.tsv
	python3 SOURCE_annotate.py  !{cosmic_cmc} !{cosmic_freqs} PRE_VARIANTS.tsv
	{ head -1 VARIANTS.tsv && tail +2 VARIANTS.tsv | sort -t$'\t' -k10,10 -k7,7nr; } > ${sample}.ONCOZOOM_VARIANTS.tsv
	rm PRE_VARIANTS.tsv
	rm VARIANTS.tsv

	# create the EVIDENCE tsv file
	echo "creating evidence tsv..."
	python3 EVIDENCE_create.py !{civic_tsv} ${sample}.ONCOZOOM_VARIANTS.tsv ${sample} 

	# create the RELEVANT VARIANTS tsv file
	echo "creating the relevant variants tsv..."
	python3 python/RELEVANT_SELECTION.py ${sample}.ONCOZOOM_VARIANTS.tsv ${sample}.ONCOZOOM_EVIDENCE.tsv ${sample}

	# rm tmp vcfs
	rm $sample.cosmic.civic.vcf.gz*
	rm $sample.cosmic.vcf.gz*
	rm $sample.ann.vcf.gz*
	rm $sample.norm.vcf.gz*
	rm $sample.somatic.vcf.gz*
done
##################################
	
