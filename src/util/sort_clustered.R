#! /usr/bin/Rscript

# A script which reads variants in the format:
# gene    contig  pos     ref     alt     effect  isof    change  other genotypes...
#
# and a file with the new sort order of genotypes in the format:
# gene
# contig
# pos
# ref
# alt
# effect
# isof
# change
# other
# genotypes...

# The output is the same as in the first input file
# but with the genotypes sorted according as in the second file

args <- commandArgs(trailingOnly = TRUE)

data_filename=args[1];
sort_filename=args[2];

cat(paste("Reference file: ", data_filename, "\n"), file=stderr());
cat(paste("Sorted samples file: ", sort_filename, "\n"), file=stderr());

data_csv=read.csv(data_filename, sep="\t", header=TRUE); # read genotyping data

sort_csv=scan(sort_filename, character()) # read desired final order of samples

sorted_data_csv <- data_csv[,sort_csv]

# output
write.table(sorted_data_csv, stdout(), quote=FALSE, sep="\t", row.names=FALSE)

cat("END\n", file=stderr());

## END