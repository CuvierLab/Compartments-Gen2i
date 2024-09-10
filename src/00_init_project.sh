# Dowload genome related files

# ce11.fa
wget --content-disposition https://ftp.ensembl.org/pub/release-99/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna_sm.toplevel.fa.gz
mkdir -p bioDB/ce11/genome/
gunzip Caenorhabditis_elegans.WBcel235.dna_sm.toplevel.fa.gz
mv Caenorhabditis_elegans.WBcel235.dna_sm.toplevel.fa bioDB/ce11/genome/

