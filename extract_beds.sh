#!/bin/bash

# From a gtf file, obtain several bed files for various gene features:
#   1. a gene bed file (6col bed file of gene regions)
#   2. a gene-level promoter bed file
#   3. a gene-level TSS bed file
#   4. an exon bed file
#   5. a gene-level first exon bed file
#   6. an intron bed file
#   7. a gene-level first intron bed file
#   8. a 5'UTR bed file
#   9. a 3'UTR bed file

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
first_exon_bed=$name.first_exon.bed      # set name of first exon bed file
intron_bed=$name.intron.bed              # set name of intron bed file
first_intron_bed=$name.first_intron.bed  # set name of first intron bed file
fiveprime_bed=$name.fiveprimeutr.bed     # set name of 5' UTR bed file
threeprime_bed=$name.threeprimeutr.bed   # set name of 3' UTR bed file
tss_bed=$name.tss.bed                    # set name of TSS bed file

#--------------------------------------------------------------------------------------------------------------------#

# sort the gtf file
cat $my_gtf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k4,4n -k5,5n"}' > my_gtf.sorted.gtf
# sort the chr_lens.txt file
cat $chr_lens | sort -k1,1 -k2,2n > chr_lens_sorted.txt

#--------------------#
# CHECK GTF FEATURES #
#--------------------#
# check to see if the above features exist in the provided input gtf file
cut -f3 $my_gtf | sort | uniq > $name.key

if grep -Fiq "five_prime_utr" $name.key
then
    has_five_prime_utr=true
else
    has_five_prime_utr=false
fi

if grep -Fiq "three_prime_utr" $name.key
then
    has_three_prime_utr=true
else
    has_three_prime_utr=false
fi
rm $name.key

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
echo "Created bed file for genes."

#----------------------------------------------#
# GET PROMOTER BED FILE from the gene bed file #
#----------------------------------------------#
bedtools flank -i $gene_bed -g $chr_lens -l $prom_size -r 0 -s > $tmp_bed
# bedtools sort
bedtools sort -i $tmp_bed > $prom_bed
rm $tmp_bed
echo "Created bed file for promoters $prom_size bp upstream from TSS."

#------------------#
# GET TSS BED FILE #
#------------------#
bedtools flank -i $gene_bed -g $chr_lens -l 1 -r 0 -s | awk 'BEGIN{FS=OFS="\t"} ($6=="+"){$2=$2+1} ($6=="+"){$3=$3+1} ($6=="-"){$2=$2-1} ($6=="-"){$3=$3-1} {print $0}' > $tmp_bed
# bedtools sort
bedtools sort -i $tmp_bed > $tss_bed
echo "Created bed file for Transcription Start Sites."

#-------------------#
# GET EXON BED FILE #
#-------------------#
style='exon'
# print only lines where 3rd field in the defined style
awk -v s_style=$style 'BEGIN{FS=OFS="\t"} $3==s_style {print}' $my_gtf | awk '{gsub(/\"|\;/,"",$10)}1' | awk '{print $1 "\t" $4 "\t" $5 "\t" $10 "\t" $6 "\t" $7}' > $tmp_bed
# sort the gene bed file with bedtools sort
bedtools sort -i $tmp_bed > $tmp_bed2
rm $tmp_bed
# remove any duplicate lines, keeping one unique entry per duplicate
awk '!a[$0]++' $tmp_bed2 > $tmp_bed
rm $tmp_bed2
# subtract 1 from start col (2nd field)
awk '{print $1 "\t" ($2 - 1) "\t" $3 "\t" $4 "\t" $5 "\t" $6}' $tmp_bed > $exon_bed
rm $tmp_bed
echo "Created bed file for exons."

#-------------------------#
# GET FIRST EXON BED FILE #
#-------------------------#
# first, separate the exon bed file into 2 files based on strand
awk '$6 == "+" { print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' $exon_bed > $name.plus.bed
awk '$6 == "-" { print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' $exon_bed > $name.minus.bed
# for the exons on the plus strand, we want to keep the first exon entry per gene
awk '!a[$4]++' $name.plus.bed > $name.plus2.bed
# for the exons on the minus strand, we want to keep the last exon entry per gene
tac $name.minus.bed | awk '!a[$4]++' | tac > $name.minus2.bed
# combine the first exons from both strands into one file
cat $name.plus2.bed $name.minus2.bed > $name.first_exons.bed
# sort it with bedtools
bedtools sort -i $name.first_exons.bed > $first_exon_bed
# remove tmp files
rm $name.plus.bed
rm $name.minus.bed
rm $name.plus2.bed
rm $name.minus2.bed
rm $name.first_exons.bed
echo "Created bed file for first exons."

#---------------------------------#
# OLD METHOD: GET INTRON BED FILE #
#---------------------------------#
## get the complement of every feature in the gtf file (can be considered intergenic regions)
## use the sorted gtf file and sorted chr_lens.txt file
#bedtools complement -i my_gtf.sorted.gtf -g chr_lens_sorted.txt > intergenic_sorted.bed
## get introns
## considered as the complement of the intergenic and exon regions
## this works because the other features in the gtf file (CDS, UTR, Start/Stop codon) are all sub-intervals of exon
## thus, everything that is not intergenic or exon must be intron
#cut -f1,2,3 $exon_bed > exon.bed # select first 3 columns from exon bed file
#bedtools complement -i <(cat exon.bed intergenic_sorted.bed | sort -k1,1 -k2,2n) -g chr_lens_sorted.txt > $intron_bed
## We need to add columns 4,5,6 to intron bed. Get this info from the 6-col exon bed file
## add a key (chr_start) as first column to intron.bed file (3 col)
#awk 'BEGIN{FS=OFS="\t"} {print $1"_"$2, $0}' $intron_bed > new_intron.bed
## add a key (chr_end) as first column to exon.bed file (6 col)
#awk 'BEGIN{FS=OFS="\t"} {print $1"_"$3, $0}' $exon_bed > new_exon.bed
## match and create final 6-col intron bed file
#awk 'BEGIN{FS=OFS="\t"} NR==FNR {a[$1]=$5;b[$1]=$6;c[$1]=$7;next}{print $0 "\t" a[$1] "\t" b[$1] "\t" c[$1]}' new_exon.bed new_intron.bed | cut -f2- > intron.bed
## For some gtf files, there will be a "gene" entry with no exon entries
## This causes the entire gene to appear in the intron bed file
## For these special cases, we will remove them from the intron bed file
## Remove any line with a blank value for the strand field
#awk '$6!=""' intron.bed > $intron_bed
## remove tmp files
#rm intron.bed
#rm intergenic_sorted.bed
#rm exon.bed
#rm new_intron.bed
#rm new_exon.bed
#echo "Created bed file for introns."

#---------------------#
# GET INTRON BED FILE #
#---------------------#
# split up the genes.bed file into plus and minus strand genes
awk '$6 == "+" { print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' $gene_bed > $name.gene.plus.bed
awk '$6 == "-" { print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' $gene_bed > $name.gene.minus.bed
# Take the complement of each strand gene bed file, using chr lengths
bedtools complement -i $name.gene.plus.bed -g chr_lens_sorted.txt > $name.plus.intergenic_sorted.bed
bedtools complement -i $name.gene.minus.bed -g chr_lens_sorted.txt > $name.minus.intergenic_sorted.bed
# split the exon bed file into plus and minus strands (only keep 3 columns)
awk '$6 == "+" { print $1"\t"$2"\t"$3}' $exon_bed > $name.exon.plus.bed
awk '$6 == "-" { print $1"\t"$2"\t"$3}' $exon_bed > $name.exon.minus.bed
# Get introns for each strand
#   considered as the complement of the intergenic and exon regions
#   this works because the other features in the gtf file (CDS, UTR, Start/Stop codon) are all sub-intervals of exon
#   thus, everything that is not intergenic or exon must be intron
bedtools complement -i <(cat $name.exon.plus.bed $name.plus.intergenic_sorted.bed | sort -k1,1 -k2,2n) -g chr_lens_sorted.txt > $name.intron.plus.bed
bedtools complement -i <(cat $name.exon.minus.bed $name.minus.intergenic_sorted.bed | sort -k1,1 -k2,2n) -g chr_lens_sorted.txt > $name.intron.minus.bed
# Add columns 4,5,6 from exon file to intron file
#   add a key (chr_start) as first column to intron.bed file (3 col)
awk 'BEGIN{FS=OFS="\t"} {print $1"_"$2, $0}' $name.intron.plus.bed > $name.new_intron.plus.bed
awk 'BEGIN{FS=OFS="\t"} {print $1"_"$2, $0}' $name.intron.minus.bed > $name.new_intron.minus.bed
#   add a key (chr_end) as first column to exon.bed file (6 col)
awk 'BEGIN{FS=OFS="\t"} {print $1"_"$3, $0}' $exon_bed > $name.new_exon.bed
# match and create final 6-col intron bed file
awk 'BEGIN{FS=OFS="\t"} NR==FNR {a[$1]=$5;b[$1]=$6;c[$1]=$7;next}{print $0 "\t" a[$1] "\t" b[$1] "\t" c[$1]}' $name.new_exon.bed $name.new_intron.plus.bed | cut -f2- > $name.intron.plus.bed
awk 'BEGIN{FS=OFS="\t"} NR==FNR {a[$1]=$5;b[$1]=$6;c[$1]=$7;next}{print $0 "\t" a[$1] "\t" b[$1] "\t" c[$1]}' $name.new_exon.bed $name.new_intron.minus.bed | cut -f2- > $name.intron.minus.bed
# combine the plus and minus intron bed files
cat $name.intron.plus.bed $name.intron.minus.bed > $name.introns.bed
# sort it with bedtools
bedtools sort -i $name.introns.bed > $name.introns_sorted.bed
# remove entries with no introns
#   For some gtf files, there will be a "gene" entry with no exon entries
#   This causes the entire gene to appear in the intron bed file
#   For these special cases, we will remove them from the intron bed file
#   Remove any line with a blank value for the strand field
awk '$6!=""' $name.introns_sorted.bed > $intron_bed
# remove tmp files
rm $name.gene.plus.bed
rm $name.gene.minus.bed
rm $name.plus.intergenic_sorted.bed
rm $name.minus.intergenic_sorted.bed
rm $name.exon.plus.bed
rm $name.exon.minus.bed
rm $name.intron.plus.bed
rm $name.intron.minus.bed
rm $name.new_intron.plus.bed
rm $name.new_intron.minus.bed
rm $name.new_exon.bed
rm $name.introns.bed
rm $name.introns_sorted.bed
echo "Created bed file for introns."

#---------------------------#
# GET FIRST INTRON BED FILE #
#---------------------------#
# first, separate the intron bed file into 2 files based on strand
awk '$6 == "+" { print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' $intron_bed > $name.plus.bed
awk '$6 == "-" { print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' $intron_bed > $name.minus.bed
# for the introns on the plus strand, we want to keep the first intron entry per gene
awk '!a[$4]++' $name.plus.bed > $name.plus2.bed
# for the introns on the minus strand, we want to keep the last intron entry per gene
tac $name.minus.bed | awk '!a[$4]++' | tac > $name.minus2.bed
# combine the first introns from both strands into one file
cat $name.plus2.bed $name.minus2.bed > $name.first_introns.bed
# sort it with bedtools
bedtools sort -i $name.first_introns.bed > $first_intron_bed
# remove tmp files
rm $name.plus.bed
rm $name.minus.bed
rm $name.plus2.bed
rm $name.minus2.bed
rm $name.first_introns.bed
echo "Created bed file for first introns."

#--------------------#
# GET 5'UTR BED FILE #
#--------------------#
if $has_five_prime_utr
then
  style="five_prime_utr"
  # print only lines where 3rd field in the defined style
  awk -v s_style=$style 'BEGIN{FS=OFS="\t"} $3==s_style {print}' $my_gtf | awk '{gsub(/\"|\;/,"",$10)}1' | awk '{print $1 "\t" $4 "\t" $5 "\t" $10 "\t" $6 "\t" $7}' > $tmp_bed
  # sort the gene bed file with bedtools sort
  bedtools sort -i $tmp_bed > $tmp_bed2
  rm $tmp_bed
  # subtract 1 from start col (2nd field)
  awk '{print $1 "\t" ($2 - 1) "\t" $3 "\t" $4 "\t" $5 "\t" $6}' $tmp_bed2 > $fiveprime_bed
  rm $tmp_bed2
  echo "Created bed file for Five Prime UTRs."
else
    echo "No five_prime_utr features found in input gtf file."
fi

#--------------------#
# GET 3'UTR BED FILE #
#--------------------#
if $has_three_prime_utr
then
  style="three_prime_utr"
  # print only lines where 3rd field in the defined style
  awk -v s_style=$style 'BEGIN{FS=OFS="\t"} $3==s_style {print}' $my_gtf | awk '{gsub(/\"|\;/,"",$10)}1' | awk '{print $1 "\t" $4 "\t" $5 "\t" $10 "\t" $6 "\t" $7}' > $tmp_bed
  # sort the gene bed file with bedtools sort
  bedtools sort -i $tmp_bed > $tmp_bed2
  rm $tmp_bed
  # subtract 1 from start col (2nd field)
  awk '{print $1 "\t" ($2 - 1) "\t" $3 "\t" $4 "\t" $5 "\t" $6}' $tmp_bed2 > $threeprime_bed
  rm $tmp_bed2
  echo "Created bed file for Three Prime UTRs."
else
  echo "No three_prime_utr features found in input gtf file."
fi

rm my_gtf.sorted.gtf
rm chr_lens_sorted.txt
