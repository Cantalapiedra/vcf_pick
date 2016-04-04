vcf_pick
========

Tools to parse VCF for specific purposes

##### vcf_filter

A script to parse a large VCF file to a smaller one, which will be faster to parse with other scripts.

Parameters:
- -c FILENAME: only contigs/chr present in this file will be output.
- -v FILENAME: only variants (contig/chr position) present in this file will be output.
- -g FILENAME: only genes present in this file will be output.
- -i FILENAME: only isoforms present in this file will be output.
- -s FILENAME: only samples present in this file will be output. Note that sample names should correspond to sample names to be output if using the "-t" option (see below).
- -t FILENAME: translation from VCF sample names to sample names to be output.

##### parse_contigs_variants

A script which outputs polymorphic variants for each contig/chr position.

`parse_contigs_variants.py test_02_snpeff.vcf -c test_contigs.list`

Parameters:
- -H FILENAME: allows to specify a file with the header fields of the VCF (to be used in case the header is not present in the VCF file.
- -c FILENAME: only contigs/chr present in this file will be output.
- -v FILENAME: only variants (contig/chr position) present in this file will be output.
- -s FILENAME: only samples present in this file will be output. Note that sample names should correspond to sample names to be output if using the "-t" option (see below).
- -t FILENAME: translation from VCF sample names to sample names to be output.
- --contigs_info FILENAME: additional fields to be output along with contigs/chr.
- -f OUTPUT_FORMAT: either summary, detail or tabular.
- --max_missing FLOAT: maximum missing genotypes to output a variant (default: 1.0).
- --max_heteros FLOAT: maximum heterozygous genotypes to output a variant (it behaves differently when changing heterozygous variants to missing variants: see --het_to_miss) (default: 1.0).
- --maf: minimum allele frequency to output a variant (default: 0.0).

Flags:
- --het_to_miss: a flag changing the behaviour when dealing with heterozygous and missing genotypes. If the flag is present and percentage of heterozygous genotypes in a variant is greater than "--max_heteros", those heterozygous genotypes will be changed to missing genotypes before assessing "--max_missing" value.
- -m: Output monomorphic variants along with polymorphic ones.
- -e: Show effects (snpEff) of variants on genes.
- -b: Output both alleles of each genotype instead of a single symbol (e.g.: 0/1 instead of "h").
- -n: Assign numeric values to the genotypes. Ignored if -b flag is present.
- -k: cluster samples using genotypes as input for R hclust.

##### parse_gene_variants

A script which outputs polymorphic variants for each gene/isoform.

`./parse_gene_variants.py test_02_snpeff.vcf -i test_isoforms.list`

Parameters:
- -H FILENAME: allows to specify a file with the header fields of the VCF (to be used in case the header is not present in the VCF file.
- -g FILENAME: only genes present in this file will be output.
- -i FILENAME: only isoforms present in this file will be output.
- -v FILENAME: only variants (contig/chr position) present in this file will be output.
- -s FILENAME: only samples present in this file will be output. Note that sample names should correspond to sample names to be output if using the "-t" option (see below).
- -t FILENAME: translation from VCF sample names to sample names to be output.
- --contigs_info FILENAME: additional fields to be output along with contigs/chr.
- --genes_info FILENAME: additional fields to be along with genes.
- -f OUTPUT_FORMAT: either summary, detail or tabular.
- --max_missing FLOAT: maximum missing genotypes to output a variant (default: 1.0).
- --max_heteros FLOAT: maximum heterozygous genotypes to output a variant (it behaves differently when changing heterozygous variants to missing variants: see --het_to_miss) (default: 1.0).
- --maf: minimum allele frequency to output a variant (default: 0.0).

Flags:
- --het_to_miss: a flag changing the behaviour when dealing with heterozygous and missing genotypes. If the flag is present and percentage of heterozygous genotypes in a variant is greater than "--max_heteros", those heterozygous genotypes will be changed to missing genotypes before assessing "--max_missing" value.
- -m: Output monomorphic variants along with polymorphic ones.
- -b: Output both alleles of each genotype instead of a single symbol (e.g.: 0/1 instead of "h").
- -n: Assign numeric values to the genotypes. Ignored if -b flag is present.
- -k: cluster samples using genotypes as input for R hclust.

File formats:

VCF_FILE: currently, only format VCFv4.1 with snpEffv4.0b has been tested.
- H FILENAME must contain a VCF formatted header, including the initial symbol ("#").
- -c FILENAME, -g FILENAME, -i FILENAME, -s FILENAME: must have a single column, each row representing an accession.
- -v FILENAME: must have 2 columns, each row containing contig/chr and position, tab separated.
- -t FILENAME: must contain 2 tab separated columns, first one corresponding to the sample name in the VCF file, second one corresponding to the sample name to be output.
- --contigs_info FILENAME: should contain the contig/chr name in the first column, followed by additional fields to be output, tab separated. Please, be sure that all the rows in this file have exactly the same number of columns.
- --genes_info FILENAME: should contain the gene/isoform name in the first column, followed by additional fields to be output, tab separated. Please, be sure that all the rows in this file have exactly the same number of columns.
- -f OUTPUT_FORMAT:
  - summary: output only the variants, without alleles or genotypes.
  - detail: below each variant, a row for each allele is shown, along with the genotypes of each allele.
  - tabular: along with each variant, columns are added to include the alleles of each sample.
In general, fields in the output are separated to each other by tabs.

##### parse_effects_genes

It simply outputs the number of occurrences of each kind of SNP effect for each isoform.

`./parse_effects_genes.py test/test_02_snpeff.vcf`

##### check_polymorph

It is a script which outputs only polymorphic variants between two specified lines.
Parameters are:
- VCF input file
- First line/genotype to be compared
- Second line/genotype to be compared
- Output heterozygous variants (yes/no)

Examples, filtering out heterozygous variants:

`./check_polymorph.sh test/test_01.vcf Parent2 Parent42 no`
