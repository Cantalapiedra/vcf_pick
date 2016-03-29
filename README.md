vcf_pick
========

Tools to parse VCF for specific purposes

##### check_polymorph

It is a script which outputs only polymorphic variants between two specified lines.
Parameters are:
- VCF input file
- First line/genotype to be compared
- Second line/genotype to be compared
- Output heterozygous variants (yes/no)

Examples, filtering out heterozygous variants:

`./check_polymorph.sh test/test_01.vcf Parent2 Parent42 no`

##### parse_effects_genes

It simply outputs the number of occurrences of each kind of SNP effect for each isoform.

`./parse_effects_genes.py test/test_02_snpeff.vcf`

##### parse_genes_variants

It allows to output variants and effects found for each gene, in three output formats:
- Summary: only the variants.
- Detail: the genotypes of each allele are appended to each of the variants.
- Tabular: variants along the alleles for all the genotypes

`./parse_genes_variants.py test_02_snpeff.vcf -H test_02_snpeff.vcf.header -i isof_list.01.tab -s samples.translation -f tabular`
TODO: test parameters
