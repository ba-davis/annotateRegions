# annotateRegions

### Dependencies
[BEDtools](https://bedtools.readthedocs.io/en/latest/) installed and accessible in PATH variable   
R

### extract_beds
The script 'extract_beds.sh' creates the following gene feature bed files from the provided gtf:
1. gene bed file
2. promoter bed file
3. tss.bed file
4. exon.bed file
5. intron.bed file
6. 5'UTR bed file
7. 3'UTR bed file

This script utilizes [BEDtools](https://bedtools.readthedocs.io/en/latest/)

Inputs:
- gtf file
- chr_lens.txt file
- size of promoters in bp (upstream from TSS)


TODO:
add sanity checks at the beginning to check if the features exist in the input gtf file
if the feature isn't found, output a warning and do not attempt to create that feature bed file

extract first intron bed file
