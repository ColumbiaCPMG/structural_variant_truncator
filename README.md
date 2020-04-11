# structural_variant_truncator

Converts PennCNV structural variants from spanning multiple regions to one SV per region. Typically this is used for converting microarray called SVs to exome only SVs.

Input files are a PennCNV formatted list of SVs and a 3 column regions file listing chromosome, start and end.

Script detects tsv vs. csv file formats and the chromosome regardless if 'chr' prefix is used.