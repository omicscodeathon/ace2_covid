#!/usr/bin/Rscript

# Take 'all' htseq-count results and melt them in to one big dataframe

#Adapted from: https://wiki.bits.vib.be/index.php/NGS_RNASeq_DE_Exercise.4

# required packages
library(tibble)

# # where are we?
cntdir <- here::here("C:/Users/user/Desktop/marion/ACE2/Counts")
pat <- ".counts.txt"
hisat2.all <- list.files(path = cntdir,
                         pattern = pat,
                         all.files = TRUE,
                         recursive = FALSE,
                         ignore.case = FALSE,
                         include.dirs = FALSE)

# we choose the 'all' series
myfiles <- hisat2.all
DT <- list()

# read each file as array element of DT and rename the last 2 cols
# we created a list of single sample tables
for (i in 1:length(myfiles) ) {
  infile = paste(cntdir, myfiles[i], sep = "/")
  DT[[myfiles[i]]] <- read.table(infile, header = F, stringsAsFactors = FALSE)
  cnts <- gsub("(.*).counts.txt", "\\1", myfiles[i])
  colnames(DT[[myfiles[i]]]) <- c("ID", cnts)
}

# merge all elements based on first ID columns
reads_count <- DT[[myfiles[1]]]

# inspect
#head(reads_count)

# we now add each other table with the ID column as key
for (i in 2:length(myfiles)) {
  y <- DT[[myfiles[i]]]
  z <- merge(reads_count, y, by = c("ID"))
  reads_count <- z
}

# ID column becomes rownames
rownames(reads_count) <- reads_count$ID
reads_count <- reads_count[,-1]

# add total counts per sample
reads_count <- rbind(reads_count,
                                 tot.counts=colSums(reads_count))

# inspect and look at the top row names!
#head(reads_count)

#tail(reads_count)

####################################
# take summary rows to a new table
# ( not starting with Tb and tmp with invert=TRUE )

# transpose table for readability

reads_count_summary <- reads_count[grep("^__", rownames(reads_count), 
                                                    perl=TRUE, invert=FALSE), ]

# review
#reads_count_summary

# transpose table
t(reads_count_summary)

# write summary to file
write.csv(reads_count_summary,
          file = here::here("C:/Users/user/Desktop/marion/ACE2/reads_count_summary.csv"),
          row.names = TRUE)

####################################
# take all data rows to a new table
reads_count <- reads_count[grep("__", rownames(reads_count), perl=TRUE, invert=TRUE), ]

# inspect final merged table
#head(reads_count, 3)

# colnames(reads_count) <- gsub("(_.*$)", "", colnames(reads_count))

# write data to files
saveRDS(reads_count, file = here::here("C:/Users/user/Desktop/marion/ACE2/reads_count.RDS"))


# reads_count <- rownames_to_column(reads_count,"transcript_id")
write.csv(reads_count, 
          file = here::here("C:/Users/user/Desktop/marion/ACE2/reads_count.csv"))

# cleanup intermediate objects
rm(y, z, i, DT, reads_count)


