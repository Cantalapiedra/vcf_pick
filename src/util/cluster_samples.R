#! /usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)

data_filename=args[1];

cat(paste("Reference file: ", data_filename, "\n"), file=stderr());

data_csv=read.csv(data_filename, sep="\t", header=TRUE);

data_csv_d = dist(t(data_csv));
data_csv_h = hclust(data_csv_d);

sorted_samples = colnames(data_csv)[data_csv_h$order];

print_sample <- function(x){
    cat(paste(x, "\n"), file=stdout());
}
void = sapply(sorted_samples, print_sample);

cat("END\n", file=stderr());

## END