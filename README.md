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

`./check_polymorph.sh test_01.vcf Parent2 Parent42 no`

##### parse_snpeff_genes

It simply outputs the number of occurences of a kind of SNP effect for each isoform.

`./parse_snpeff_genes.py test_02_snpeff.vcf`
