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

##### parse_contigs_variants

As the previous script, but provided of a list of contigs (chr) instead of genes.
It does not depend on having snpEff annotation.

##### parse_variants

As the previous ones, but provided of a list of variants (contig or chr, position).
It does not depend on having snpEff annotation.

The three previous scripts can be provided of additional parameters:
- A list of samples to output, excluding the ones not in such list.
- A list of tuples to translate the name from the VCF file (often a codename) to another name (a common sample name).
- Three percentage thresholds to filter out variants:
  - Maximum heterozygous genotypes allowed.
  - Maximum missing data genotypes allowed.
  - Minimum allele frequency of each allele allowed.
In addition, monomorphic variants are excluded from the output, except if the "-m" flag is included in the options.
Furthermore, parse_contigs_variants and parse_variants can be run with the "-e" option (to show snpEff annotation if present in the VCF file) or without it (to hide the snpEff annotation).
Finally, all three scripts can produce three output formats (specified along the option "-f"):
- "summary": output only the variants, without alleles and genotypes.
- "detail": below each variant, a row for each allele is shown, along with the genotypes of each allele.
- "tabular": along each variant, a column is added to include the allele of each genotype.

