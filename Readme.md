#GFF3 Manipulator
GFF3 (General Feature Format version 3) files are annotation files that, when paired with the appropriate genome file, can provide a wealth of information, such as the number of specific domains present in various segments of the genome (genes, CDS, RNA, etc.). This repository offers a collection of Python scripts that can be useful for various tasks required in genomic analysis.

#Conversion from GTF to GFF3
GTF (General Transfer Format) files are another type of annotation file and, along with GFF3, are among the most widely used formats in genomics. Since this repository focuses on the use of GFF3 annotation files, a script is provided for converting GTF files to GFF3 format.

#Adding Introns
Many existing annotation files lack intron coordinates. Introns are crucial for several reasons: they play significant roles in gene expression regulation, alternative splicing, and the evolution of new transcipts. Moreover, the absence of intron information can affect the visualization and analysis of genomic data in tools like Geneious or Omicsbox. This script adds intron coordinates to GFF3 files, enhancing the annotation's completeness.

#Extracting Sequences
Databases typically provide FASTA files for each species, containing gene, transcript, and/or CDS sequences. However, if such files are missing or an intron FASTA file is needed, this script uses the coordinates from the GFF3 file to extract the sequences and create a corresponding FASTA file.

#Domain Search
Domains are fundamental units of protein structure and function. They are essential for understanding protein interactions, functions, and evolutionary history. This script extracts the number of matches for specified domains in each sequence provided in the FASTA file. When used in conjunction with the sequence extraction script, it can identify and locate domains of interest throughout the genome.
