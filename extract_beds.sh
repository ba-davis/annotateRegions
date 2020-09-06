#!/bin/bash

# TODO: clean up intro bed filee output
#          - need GeneID of intron
#          - strand? not necessary for stran symmetrical cpg methylation, but for future uses

# From a gtf file, obtain 7 bed files for various gene features:
#   1. a gene bed file (6col bed file of gene regions)
#   2. a gene-level promoter bed file
#   3. a gene-level TSS bed file
#   4. an exon bed file
#   5. an intron bed file
#   6. a 5'UTR bed file
#   7. a 3'UTR bed file

#--------------------------------------------------------------------------------------------------------------------#

# User-specified input variables
my_gtf=$1                                # input gtf file
chr_lens=$2                              # input chromosome lengths text file (used to get promoters, tss, and introns)
prom_size=$3                             # input size of promoters in bp (upstream from TSS)

# Script variables
name=${my_gtf%.*}                        # get basename for gene feature bed files
tmp_bed=tmp.bed                          # tmp file
tmp_bed2=tmp2.bed                        # tmp file
gene_bed=$name.genes.bed                 # set name of gene bed file
prom_bed=$name.promoters.$prom_size.bed  # set name of promoter bed file
exon_bed=$name.exon.bed                  # set name of exon bed file
intron_bed=$name.intron.bed              # set name of intron bed file
fiveprime_bed=$name.fiveprimeutr.bed     # set name of 5' UTR bed file
threeprime_bed=$name.threeprimeutr.bed   # set name of 3' UTR bed file
tss_bed=$name.tss.bed                    # set name of TSS bed file

#--------------------------------------------------------------------------------------------------------------------#

# sort the gtf file
cat $my_gtf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k4,4n -k5,5n"}' > my_gtf.sorted.gtf
# sort the chr_lens.txt file
cat $chr_lens | sort -k1,1 -k2,2n > chr_lens_sorted.txt

#-------------------#
# GET GENE BED FILE #
#-------------------#
style="gene"
awk -v s_style=$style 'BEGIN{FS=OFS="\t"} $3==s_style {print}' $my_gtf | awk '{gsub(/\"|\;/,"",$10)}1' | awk '{print $1 "\t" $4 "\t" $5 "\t" $10 "\t" $6 "\t" $7}' > $tmp_bed
# sort the gene bed file with bedtools sort
bedtools sort -i $tmp_bed > $tmp_bed2
rm $tmp_bed
# subtract 1 from start col (2nd field)
awk '{print $1 "\t" ($2 - 1) "\t" $3 "\t" $4 "\t" $5 "\t" $6}' $tmp_bed2 > $gene_bed
rm $tmp_bed2


#----------------------------------------------#
# GET PROMOTER BED FILE from the gene bed file #
#----------------------------------------------#
bedtools flank -i $gene_bed -g $chr_lens -l $prom_size -r 0 -s > $tmp_bed
# bedtools sort
bedtools sort -i $tmp_bed > $prom_bed
rm $tmp_bed

#------------------#
# GET TSS BED FILE #
#------------------#
bedtools flank -i $gene_bed -g $chr_lens -l 1 -r 0 -s | awk 'BEGIN{FS=OFS="\t"} ($6=="+"){$2=$2+1} ($6=="+"){$3=$3+1} ($6=="-"){$2=$2-1} ($6=="-"){$3=$3-1} {print $0}' > $tmp_bed
# bedtools sort
bedtools sort -i $tmp_bed > $tss_bed

#-------------------#
# GET EXON BED FILE #
#-------------------#
style='exon'
# print only lines where 3rd field in the defined style
awk -v s_style=$style 'BEGIN{FS=OFS="\t"} $3==s_style {print}' $my_gtf | awk '{gsub(/\"|\;/,"",$10)}1' | awk '{print $1 "\t" $4 "\t" $5 "\t" $10 "\t" $6 "\t" $7}' > $tmp_bed
# sort the gene bed file with bedtools sort
bedtools sort -i $tmp_bed > $tmp_bed2
rm $tmp_bed
# subtract 1 from start col (2nd field)
awk '{print $1 "\t" ($2 - 1) "\t" $3 "\t" $4 "\t" $5 "\t" $6}' $tmp_bed2 > $exon_bed
rm $tmp_bed2

#---------------------#
# GET INTRON BED FILE #
#---------------------#
# get the complement of every feature in the gtf file (can be considered intergenic regions)
# use the sorted gtf file and sorted chr_lens.txt file
bedtools complement -i my_gtf.sorted.gtf -g chr_lens_sorted.txt > intergenic_sorted.bed
# get introns
# considered as the complement of the intergenic and exon regions
# this works because the other features in the gtf file (CDS, UTR, Start/Stop codon) are all sub-intervals of exon
# thus, everything that is not intergenic or exon must be intron
cut -f1,2,3 $exon_bed > exon.bed # select first 3 columns from exon bed file
bedtools complement -i <(cat exon.bed intergenic_sorted.bed | sort -k1,1 -k2,2n) -g chr_lens_sorted.txt > $intron_bed
rm intergenic_sorted.bed
rm exon.bed
# get introns..currently performed via genomation

#--------------------#
# GET 5'UTR BED FILE #
#--------------------#
style="five_prime_utr"
# print only lines where 3rd field in the defined style
awk -v s_style=$style 'BEGIN{FS=OFS="\t"} $3==s_style {print}' $my_gtf | awk '{gsub(/\"|\;/,"",$10)}1' | awk '{print $1 "\t" $4 "\t" $5 "\t" $10 "\t" $6 "\t" $7}' > $tmp_bed
# sort the gene bed file with bedtools sort
bedtools sort -i $tmp_bed > $tmp_bed2
rm $tmp_bed
# subtract 1 from start col (2nd field)
awk '{print $1 "\t" ($2 - 1) "\t" $3 "\t" $4 "\t" $5 "\t" $6}' $tmp_bed2 > $fiveprime_bed
rm $tmp_bed2

#--------------------#
# GET 3'UTR BED FILE #
#--------------------#
style="three_prime_utr"
# print only lines where 3rd field in the defined style
awk -v s_style=$style 'BEGIN{FS=OFS="\t"} $3==s_style {print}' $my_gtf | awk '{gsub(/\"|\;/,"",$10)}1' | awk '{print $1 "\t" $4 "\t" $5 "\t" $10 "\t" $6 "\t" $7}' > $tmp_bed
# sort the gene bed file with bedtools sort
bedtools sort -i $tmp_bed > $tmp_bed2
rm $tmp_bed
# subtract 1 from start col (2nd field)
awk '{print $1 "\t" ($2 - 1) "\t" $3 "\t" $4 "\t" $5 "\t" $6}' $tmp_bed2 > $threeprime_bed
rm $tmp_bed2

rm my_gtf.sorted.gtf
rm chr_lens_sorted.txt
