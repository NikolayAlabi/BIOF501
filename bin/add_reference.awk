#! /usr/bin/awk -f
BEGIN { FS = "\t"; OFS = FS }
# Print header info as-is
/^##/ { print; next }
# Add sample named "Reference" to the list of samples
/^#CHROM/ { print $0"\tReference"; next }
# Add homozygous reference allele to every locus
{ print $0"\t0:255,0,10" }
