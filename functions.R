

library(dplyr)
library(writexl)
#library(plotrix)
library(ggplot2)

# Strategy for use of functions:
# Run bedtools closest on the input target regions bed file and the given gene feature bed file            (bedtools_closest)
# Parse the bedtools closest output results to print stats on overlapping features and obtain clean output (clean_overlaps)
# Perform the above 2 steps for each gene feature, merge them all into one table                           (annotate_overlaps)
# Add gene info (matching GeneID to gene name and info)                                                    (add_gene_info)
# Add result info (such as DMR info per input target region)                                               (add_result_info)
# Sort results by unqID (DMR_ID)                                                                           (sort_ids)
# Export results to an excel file                                                                          (export_excel)

# TODO:
# add functionality to work with bed6 input file, and to consider strandedness when finding overlaps
# make certain features optional (not all gtfs have 5'UTR and 3'UTR)


# bedtools sort function
# bed1: path to bed file to sort
# outputs a sorted bed file in the same directory as input bed file
bedtools_sort <- function(bed1) {

  # create name for sorted bed file
  out <- paste0(gsub('.{4}$', '', bed1), ".annotR_sorted.bed")

  # don't use scientific notation when writing out
  options(scipen =99)

  # create the command string and call the command using system()
  command=paste("bedtools sort", "-i", bed1, ">" , out, sep=" ")
  #cat(command,"\n")
  try(system(command))

  return(out)
}

#bedtools closest function
# bed1: path to target regions bed file to be annotated
# bed2: path to bed file of gene part regions, used for annotating
bedtools_closest <- function(bed1, bed2) {

  # create temp output file
  out=tempfile()
  options(scipen =99) # not to use scientific notation when writing out

  # create the command string and call the command using system()
  command=paste("bedtools closest", "-a", bed1, "-b", bed2, "-d", ">" , out, sep=" ")
  #cat(command,"\n")
  try(system(command))

  # read in the output tmp file
  res=read.table(out,header=F)

  unlink(out)
  return(res)
}

# Clean bedtools closest output
# df: output of the above bedtools closest function
# ID1: column number of unqID for target regions
# ID2: column number of GeneID
# dis: column number of "distance"
# feature: type of features we are checking for overlaps of
clean_overlaps <- function(df, ID1, ID2, dis, feature) {

  # set colname of unqID column
  colnames(df)[ID1] <- "unqID"

  # how many target regions
  n_targs <- length(unique(df[ ,ID1]))

  # how many target regions overlap a feature (distance of 0)
  df_zero <- df[df[[dis]] == "0", ]
  n_overlap <- length(unique(df_zero[[ID1]]))
  perc_overlap <- round((n_overlap / n_targs)*100, 1)

  print(paste0(n_overlap, " out of ", n_targs, " target regions ", "(", perc_overlap, "%) overlap an ", feature, "."))

  # How many target features overlap a feature belonging to multiple genes
  # remove duplicate rows (where both unqID and GeneID is the same, keeps one)
  df_zero_rmdup <- df_zero[!duplicated(df_zero[c(ID1,ID2)]),]
  # keep rows where unqID appears more than once (and must have a different GeneID per above)
  check <- as.data.frame(df_zero_rmdup %>% group_by(unqID) %>% filter(n() > 1))

  #perc_multiple <- round((nrow(check) / n_overlap)*100, 1)
  
  print(paste0("Of these ", n_overlap, " target regions, ", length(unique(check$unqID)), " overlap ", feature, "s from more than one gene."))

  # remove any rows where the same values are in the unqIDa and GeneID cols
  df_rmdup <- df[!duplicated(df[c(ID1,ID2)]),]

  return(df_rmdup)
}


# function to obtain cleaned bedtools closest output per feature
# infile: path to sorted target regions bed file
# feature_bed: path to feature bed file
# feature_name: name of feature
# unqID_col: column number of unique name for target regions from bedtools closest output
# geneID_col: column number of feature GeneID from bedtools closest output
# dis_col: column number of distance from bedtools closest output
overlap_feature <- function(infile, feature_bed, feature_name, unqID_col=4, geneID_col=8, dis_col=11, TSS=FALSE) {

  # Execute bedtools closest on target regions bed file and feature bed file
  closest_res <- bedtools_closest(infile, feature_bed)
  colnames(closest_res)[geneID_col] <- paste0(feature_name, ".GeneID")

  # Produce the cleaned res table (remove duplicates where unqID and GeneID are the same)
  res_clean <- clean_overlaps(df=closest_res,
                              ID1=unqID_col,
                              ID2=geneID_col,
                              dis=dis_col,
                              feature=feature_name
  )

  # If the distance column is not zero (not overlapping) set the GeneID to "NA"
  # We want this for all features other than TSS
  if (TSS==FALSE) {
    res_clean[,geneID_col] <- as.character(res_clean[,geneID_col])
    res_clean[,geneID_col][res_clean[,dis_col] != 0] <- "NA"
    # Remove Unnecessary Columns
    res_final <- res_clean[ ,c(unqID_col,geneID_col)]
  }
  else if (TSS==TRUE) {
    res_clean[,geneID_col] <- as.character(res_clean[,geneID_col])
    # Remove Unnecessary Columns
    res_final <- res_clean[ ,c(unqID_col,geneID_col,dis_col)]
    colnames(res_final)[3] <- "TSS.dist"
  }

  return(res_final)
}


# Function to produce annotation table containing IDs of target regions and overlapping exons, introns, promoters, closest TSS
annotate_overlaps <- function(infile, bed4=TRUE, gene_bed6, exon_bed6, first_exon_bed6, intron_bed6, first_intron_bed6, promoter_bed6, tss_bed6, fiveprime_bed6=NULL, threeprime_bed6=NULL) {

  # sort the input bed file and store new sorted bed file
  # this ensures the input bed file is sorted, and this new sorted bed file will be used for bedtools functions
  new_bed <- bedtools_sort(infile)

  # if bed4=TRUE (default, input bed file is 4col), then set col numbers of bedtools closest output
  if (bed4==TRUE) {
    unqID_col <- 4
    geneID_col <- 8
    dis_col <- 11
  }

  #--------------------------------------------------------------#
  # Find Overlaps of input regions with Required Feature Regions #
  #--------------------------------------------------------------#
  # find overlaps with genes
  gene_df <- overlap_feature(new_bed, gene_bed6, "Gene", unqID_col, geneID_col, dis_col)
  # find overlaps with exons
  exon_df <- overlap_feature(new_bed, exon_bed6, "Exon", unqID_col, geneID_col, dis_col)
  # find overlaps with first exons
  first_exon_df <- overlap_feature(new_bed, first_exon_bed6, "First_Exon", unqID_col, geneID_col, dis_col)
  # find overlaps with introns
  intron_df <- overlap_feature(new_bed, intron_bed6, "Intron", unqID_col, geneID_col, dis_col)
  # find overlaps with first_introns
  first_intron_df <- overlap_feature(new_bed, first_intron_bed6, "First_Intron", unqID_col, geneID_col, dis_col)
  # find overlaps with promoters
  prom_df <- overlap_feature(new_bed, promoter_bed6, "Promoter", unqID_col, geneID_col, dis_col)
  # find overlaps with TSSes
  tss_df <- overlap_feature(new_bed, tss_bed6, "TSS", unqID_col, geneID_col, dis_col, TSS=TRUE)

  #------------------------------------#
  # Optional, 5'UTR and 3'UTR Overlaps #
  #------------------------------------#
  # if a 5'UTR bed6 file is provided, find overlaps
  if (!(is.null(fiveprime_bed6))) {
    fiveprime_df <- overlap_feature(new_bed, fiveprime_bed6, "Five_Prime_UTR", unqID_col, geneID_col, dis_col)
  }

  # if a 3'UTR bed6 file is provided, find overlaps
  if (!(is.null(threeprime_bed6))) {
    threeprime_df <- overlap_feature(new_bed, threeprime_bed6, "Three_Prime_UTR", unqID_col, geneID_col, dis_col)
  }

  #return(threeprime_df)

  #-------#
  # MERGE #
  #-------#
  # merge the multiple feature overlap dfs into one df
  # dfs: gene_df
  #      exon_df
  #      first_exon_df
  #      intron_df
  #      first_intron_df
  #      prom_df
  #      tss_df
  #      fiveprime_df
  #      threeprime_df
  # Merge the required gene feature overlaps, and add an intergenic column
  # should result in a df with 9 columns:
  #  <unqID><Gene.GeneID><Exon.GeneID><First_Exon.GeneID><Intron.GeneID><First_Intron.GeneID><Promoter.GeneID><TSS.GeneID><Intergenic>
  merge1 <- merge(gene_df, exon_df, by="unqID", all.x=T, sort=T)
  merge2 <- merge(merge1, first_exon_df, by="unqID", all.x=T, sort=T)  # 37 rows
  merge3 <- merge(merge2, intron_df, by="unqID", all.x=T, sort=T)       # 43 rows
  merge4 <- merge(merge3, first_intron_df, by="unqID", all.x=T, sort=T) # 43 rows
  merge5 <- merge(merge4, prom_df, by="unqID", all.x=T, sort=T)         # 43 rows
  merge6 <- merge(merge5, tss_df, by="unqID", all.x=T, sort=T)          # 43 rows
  # add an intergenic column if exon, intron, and promoter are all NA
  merge6$Intergenic <- "NA"
  merge6$Intergenic <- with(merge6, ifelse(Exon.GeneID == 'NA' & Intron.GeneID == 'NA' & Promoter.GeneID == 'NA', 'Intergenic', 'NA'))

  # If no 5'UTR or 3'UTR files were provided, then this is the final df
  # add else if's to include either 5'UTR, 3'UTR, or both
  if (is.null(fiveprime_bed6) & is.null(threeprime_bed6)) {
    # place Intergenic column before the TSS columns
    merge6 <- merge6[ ,c(1,2,3,4,5,6,7,10,8,9)]
    return(merge6)
  }
  else if (!(is.null(fiveprime_bed6)) & is.null(threeprime_bed6)) {
    merge7 <- merge(merge6, fiveprime_df, by="unqID", all.x=T, sort=T)
    merge7 <- merge7[ ,c(1,2,3,4,5,6,7,11,8,9,10)]
    return(merge7)
  }
  else if (!(is.null(threeprime_bed6)) & is.null(fiveprime_bed6)) {
    merge7 <- merge(merge6, threeprime_df, by="unqID", all.x=T, sort=T)
    merge7 <- merge7[ ,c(1,2,3,4,5,6,7,11,8,9,10)]
    return(merge7)
  }
  else {
    merge7 <- merge(merge6, fiveprime_df, by="unqID", all.x=T, sort=T)
    merge8 <- merge(merge7, threeprime_df, by="unqID", all.x=T, sort=T)
    merge8 <- merge8[ ,c(1,2,3,4,5,6,7,11,12,8,9,10)]
    return(merge8)
  }
}


# Function to add gene info to the annotation table
# usually biomart info
# first column must be the GeneIDs
# assumes headers are present
# five_prime_utr: defaults to False, set to true if it's present
# three_prime_utr: defaults to False, set to true if it's present
add_gene_info <- function(df, gene_info_file, five_prime_utr=FALSE, three_prime_utr=FALSE) {
  # read in the gene info
  gene_dat <- read.delim(gene_info_file, header=T)
  colnames(gene_dat)[1] <- "GeneID"

  # remove any duplicate gene ids
  gene_dat <- gene_dat[!duplicated(gene_dat$GeneID), ]

  #---------------------------------
  # MERGE to progressively add gene info to feature overlap GeneID's
  #---------------------------------
  # merge on Gene.GeneID
  colnames(gene_dat) <- paste("Gene",colnames(gene_dat), sep=".")
  mydat1 <- merge(df, gene_dat, by="Gene.GeneID", all.x=TRUE)

  # merge on exon gene ids
  colnames(gene_dat) <- gsub("Gene.", "Exon.", colnames(gene_dat))
  colnames(gene_dat)[1] <- "Exon.GeneID"
  mydat2 <- merge(mydat1, gene_dat, by="Exon.GeneID", all.x=TRUE)

  # merge on first exon gene ids
  colnames(gene_dat) <- paste0("First_", colnames(gene_dat))
  mydat3 <- merge(mydat2, gene_dat, by="First_Exon.GeneID", all.x=TRUE)

  # merge on intron gene ids
  colnames(gene_dat) <- gsub("First_Exon", "Intron", colnames(gene_dat))
  mydat4 <- merge(mydat3, gene_dat, by="Intron.GeneID", all.x=TRUE)

  # merge on first intron gene ids
  colnames(gene_dat) <- paste0("First_", colnames(gene_dat))
  mydat5 <- merge(mydat4, gene_dat, by="First_Intron.GeneID", all.x=TRUE)

  # merge on promoter gene ids
  colnames(gene_dat) <- gsub("First_Intron", "Promoter", colnames(gene_dat))
  mydat6 <- merge(mydat5, gene_dat, by="Promoter.GeneID", all.x=TRUE)

  # merge on TSS gene ids
  colnames(gene_dat) <- gsub("Promoter", "TSS", colnames(gene_dat))
  mydat7 <- merge(mydat6, gene_dat, by="TSS.GeneID", all.x=TRUE)

  # if there is no 5' or 3'UTR provided, this is the final df
  if (five_prime_utr==FALSE & three_prime_utr==FALSE) {
    # organize columns
    mydat7 <- mydat7[ ,c(8,7,6,5,4,3,2,9,1,10,11:ncol(mydat7))]
    # sort by unqID
    mydat7 <- mydat7[order(mydat7$unqID), ]
    return(mydat7)
  }
  
  # if 5'UTR is true
  if (five_prime_utr==TRUE & three_prime_utr==FALSE) {
    # merge on five prime utr gene ids
    colnames(gene_dat) <- gsub("TSS", "Five_Prime_UTR", colnames(gene_dat))
    mydat8 <- merge(mydat7, gene_dat, by="Five_Prime_UTR.GeneID", all.x=TRUE)
    # organize columns
    mydat8 <- mydat8[ ,c(9,8,7,6,5,4,3,1,11,2,10,12:ncol(mytab8))]
    # sort by unqID
    mydat8 <- mydat8[order(mydat8$unqID), ]
    return(mydat8)
  }
  # if 3'UTR is true
  if (five_prime_utr==FALSE & three_prime_utr==TRUE) {
  # merge on three prime utr gene ids
  colnames(gene_dat) <- gsub("TSS", "Three_Prime_UTR", colnames(gene_dat))
  mydat8 <- merge(mydat7, gene_dat, by="Three_Prime_UTR.GeneID", all.x=TRUE)
  # organize columns
  mydat8 <- mydat8[ ,c(9,8,7,6,5,4,3,1,11,2,10,12:ncol(mytab8))]
  # sort by unqID
  mydat8 <- mydat8[order(mydat8$unqID), ]
  return(mydat8)
  }
  # if 5'UTR and 3'UTR are true
  if (five_prime_utr==TRUE & three_prime_utr==TRUE) {
    # merge on five prime utr gene ids
    colnames(gene_dat) <- gsub("TSS", "Five_Prime_UTR", colnames(gene_dat))
    mydat8 <- merge(mydat7, gene_dat, by="Five_Prime_UTR.GeneID", all.x=TRUE)
    # merge on three prime utr gene ids
    colnames(gene_dat) <- gsub("Five", "Three", colnames(gene_dat))
    mydat9 <- merge(mydat8, gene_dat, by="Three_Prime_UTR.GeneID", all.x=TRUE)
    # organize columns
    mydat9 <- mydat9[ ,c(10,9,8,7,6,5,4,2,1,12,3,11,13:ncol(mydat9))]
    # sort by unqID
    mydat9 <- mydat9[order(mydat9$unqID), ]
    return(mydat9)
  }
}


# function to sort unqID named like "DMR_1"
sort_ids <- function(res) {
  # read in the results file which has desired order of unqIDs in first column
  targ_bed <- read.delim(res, header=T)

  # create a mapping of the unqID from the results bed file to it's ordered number
  # assumes its already in the desired order
  sort_df <- targ_bed
  sort_df$sort_order <- c(1:nrow(sort_df))
  sort_df <- sort_df[ ,c(1,ncol(sort_df))]
  colnames(sort_df)[1] <- "unqID"

  return(sort_df)
}

# FUNCTION to add result info to the annot df
# the result info may be from a differential analysis and contain further info about target regions
#   1st column is unqID (or equivalent, function changes colname)
#   if supplying coords, have as cols 2,3,4
#   2nd column must have colname of 'chr'
add_result_info <- function(df, result_file, UCSC=T) {
  # read in the result file
  res <- read.delim(result_file, header=T)
  colnames(res)[1] <- "unqID"

  #if (keep_chrs) {
    # remove the 0-based coords from the current annot.df
  #  df[ ,c(2,3,4)] <- NULL
    # merge on unqID
  #  mydat <- merge(res, df, by="unqID", all.x=TRUE)
  #}

  # add result info to annot df with merge
  mydat <- merge(res, df, by="unqID", all.x=TRUE)

  # sort unqID
  key <- sort_ids(result_file)
  mydat2 <- merge(mydat, key, by="unqID", all.x=T, sort=T)
  mydat2 <-mydat2[order(mydat2$sort_order), ]
  # remove sort column
  mydat2 <- mydat2[ ,-ncol(mydat2)]

  if (UCSC) {
    mydat2$chr <- paste0("chr", mydat2$chr)
  }

  # change <NA> values to "NA" to show up in excel
  # replace factors as characters
  i <- sapply(mydat2, is.factor)
  mydat2[i] <- lapply(mydat2[i], as.character)
  mydat2[is.na(mydat2)] <- "NA"

  return(mydat2)
}

# FUNCTION to export the results as excel file
#   df: annot.df (full output, gene info)
#   name: desired comparison/output file name prefix
#   remove.NA: default TRUE, remove DMRs which are NA in each feature
#   type: default "full", type of annot.df, further direction
export_excel <- function(df, name, remove.NA=TRUE, type="full") {

  # Split up into separate dfs based on gene feature type for separate sheets of excel file
  df.diff <- df[!duplicated(df$unqID), c(1:7)] # diff results (unique target region entries only) 
  df1 <- df[ ,c(1:16)]            # summary
  df2 <- df[ ,c(1:8,17)]          # genes
  df3 <- df[ ,c(1:7,9,18)]        # exons
  df4 <- df[ ,c(1:7,10,19)]       # first_exons
  df5 <- df[ ,c(1:7,11,20)]       # introns
  df6 <- df[ ,c(1:7,12,21)]       # first_introns
  df7 <- df[ ,c(1:7,13,22)]       # promoters
  df8 <- df[ ,c(1:7,14:16,23)]    # intergenic
  
  #df5 <- df[ ,c(1:7,11,31:35)]    # fiveprimeutr
  #df6 <- df[ ,c(1:7,12,36:40)]    # threeprimeutr
  #df7 <- df[ ,c(1:7,13:15,41:45)] # intergenic

  # remove rows from each feature type df with duplicate unqID and featureID
  df2 <- df2[!duplicated(df2[c("unqID","Gene.GeneID")]),]
  df3 <- df3[!duplicated(df3[c("unqID","Exon.GeneID")]),]
  df4 <- df4[!duplicated(df4[c("unqID","First_Exon.GeneID")]),]
  df5 <- df5[!duplicated(df5[c("unqID","Intron.GeneID")]),]
  df6 <- df6[!duplicated(df6[c("unqID","First_Intron.GeneID")]),]
  df7 <- df7[!duplicated(df7[c("unqID","Promoter.GeneID")]),]
  df8 <- df8[!duplicated(df8[c("unqID","TSS.GeneID")]),]

  # for each feature df, remove rows with NA value in that feature (if parameter remove.NA is TRUE)
  if (remove.NA) {
    df2 <- df2[df2$Gene.GeneID != "NA", ]
    df3 <- df3[df3$Exon.GeneID != "NA", ]
    df4 <- df4[df4$First_Exon.GeneID != "NA", ]
    df5 <- df5[df5$Intron.GeneID != "NA", ]
    df6 <- df6[df6$First_Intron.GeneID != "NA", ]
    df7 <- df7[df7$Promoter.GeneID != "NA", ]
    df8 <- df8[df8$Intergenic != "NA", ]
  }

  # combine into a list of dfs
  mydat_list <- list(df.diff, df1, df2, df3, df4, df5, df6, df7, df8)
  names(mydat_list) <- c("diff results", "annot summary", "gene", "exon", "first exon", "intron", "first intron", "promoter", "intergenic")

  # export excel file
  write_xlsx(mydat_list, path=paste0(name, ".annot.xlsx"))
}
