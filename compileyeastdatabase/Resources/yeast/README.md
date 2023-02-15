# barcode_details.json

An example of a barcode details configuration schema

# chr_map

this is a csv file with the following columns:

- refseq
- igenomes
- ensembl
- ucsc
- mitra
- seqlength
- numbered

where except for seqlength, which stores the length of each sequence, 
the other columns are various naming conventions for the yeast chromosomes. 
Note that there are additional plasmids in the yeast genome used for alignment 
from the mitra lab AND this genome has a deleted gene (adh1). The deletion in 
the fasta doesn't affect coordinates (replaced by Ns? not sure, but this 
was checked at one point. There is no documentation in hte mitra resources 
on how this fasta was created).

The purpose of this file is to provide a way of mapping between various 
reference files

# minus_adh1_2015.qbed

This is a background file from 2015 -- it is the background file that has 
been used in the Brent lab. However, the Mitra lab has recommended the use 
of the other background file, dSir4. Note that it is unclear if the 
genome would be the same for the two -- not sure if adh1 should be replaced 
in the genome and sir4 removed, for instance in order to accurately use the 
dSir4 background instead of the adh1 background

# orf_coding_all_R61-1-1_20080606.promoter_-700bpto0bp_with_ucsc_seqnames_common_names_coord_corrected_systematic.bed7

This is the promoter file created by Yiming. The original file had a mix of 
common and systematic names, and the start/end coordinates were different 
depending on the strand. This file has been modified so that the start/end 
coordinates adhere to typical bed format standards, and the common and 
systematic names have been separated into their own columns. Promoters are 
defined as 700bp upstream to the TSS, though how the TSS are found is 
unclear beyond what is stated in the DTO methods section (ie -- this file 
is unreproducible, was formatted in a very quesitonable way to begin with, 
and the promoter regions themselves when examined have a number of 
questionable choices)

# regions_not_orf_from_mitra.bed

This is the Mitra pipeline "promoter" definition file -- not really 
intentionally promoters, but rather "not orf" regions of the yeast genome

# S288C_dSir4_Background.qbed

This is an updated background data file which the Mitra lab recommends using 
in preference to the 2015 adh1 background data. However, in the brent lab 
this change has not been made as of 2023.