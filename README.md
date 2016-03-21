vcf_pic
=======

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
