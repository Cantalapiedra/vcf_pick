# vcf_pick
Tools to parse VCF for specific purposes

check_polymorph is a script which outputs only polymorphic variants between two specified lines

no heteros allowed
./check_polymorph.sh test_01.vcf Parent2 Parent42 no

allow heteros
./check_polymorph.sh test_01.vcf Parent1 Parent42 yes
