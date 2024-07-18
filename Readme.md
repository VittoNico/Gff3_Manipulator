## GFF3 Manipulator
GFF3 (General Feature Format version 3) files are annotation files that, when paired with the appropriate genome file, can provide a wealth of information, such as the number of specific domains present in various segments of the genome (genes, CDS, RNA, etc.). This repository offers a collection of Python scripts that can be useful for various tasks required in genomic analysis. Each script requires that you provide the coordinate of the file you want to utilize. Use a Text manipulator program to doing that

# Conversion from GTF to GFF3
GTF (General Transfer Format) files are another type of annotation file and, along with GFF3, are among the most widely used formats in genomics. Since this repository focuses on the use of GFF3 annotation files, a script is provided for converting GTF files to GFF3 format.
```python
import re

def gtf_to_gff3(input_file, output_file):
    gff3_lines = []
    attributes_pattern = re.compile(r'(\w+) "([^"]+)"')

    with open(input_file, 'r') as infile:
        for line in infile:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            feature_type = fields[2]
            start, end = int(fields[3]), int(fields[4])
            strand = fields[6]
            attributes = dict(attributes_pattern.findall(fields[8]))

            start += 1  # GTF is 1-based, GFF3 is 0-based, adjust accordingly

            gff3_line = [fields[0], fields[1], feature_type, str(start), str(end), '.', strand, '.']

            attributes_str = ';'.join([f"{key}={value}" for key, value in attributes.items()])
            gff3_line.append(attributes_str)

            gff3_lines.append('\t'.join(gff3_line))

    with open(output_file, 'w') as outfile:
        outfile.write("##gff-version 3\n")
        outfile.write('\n'.join(gff3_lines))

if __name__ == "__main__":
    input_gtf = "/path/to/Annotation_file.gtf"
    output_gff3 = "/path/to/Annotation_file.gff3"
    gtf_to_gff3(input_gtf, output_gff3)
```

# Adding Introns
Many existing annotation files lack intron coordinates. Introns are crucial for several reasons: they play significant roles in gene expression regulation, alternative splicing, and the evolution of new transcipts. Moreover, the absence of intron information can affect the visualization and analysis of genomic data in tools like Geneious or Omicsbox. This script adds intron coordinates to GFF3 files, enhancing the annotation's completeness.

```python
from collections import defaultdict

def add_introns_to_gff3(gff3_file, output_file):
    gene_exons = defaultdict(list)

    with open(gff3_file, 'r') as gff3_handle:
        for line in gff3_handle:
            if line.startswith('#'):
                continue

            columns = line.strip().split('\t')
            feature_type = columns[2]
            start = int(columns[3])
            end = int(columns[4])
            attributes = dict(item.split('=') for item in columns[8].split(';'))

            if feature_type == 'exon':
                gene_id = attributes.get('Parent')
                gene_exons[gene_id].append((start, end))

    with open(output_file, 'w') as output_handle:
        with open(gff3_file, 'r') as gff3_handle:
            for line in gff3_handle:
                output_handle.write(line)
                if line.startswith('#'):
                    continue

                columns = line.strip().split('\t')
                feature_type = columns[2]
                attributes = dict(item.split('=') for item in columns[8].split(';'))

                if feature_type == 'gene':
                    gene_id = attributes.get('ID')
                    exons = sorted(gene_exons.get(gene_id, []))
                    introns = [(exons[i][1] + 1, exons[i + 1][0] - 1) for i in range(len(exons) - 1)]

                    for intron_start, intron_end in introns:
                        intron_line = '\t'.join([columns[0], columns[1], 'intron', str(intron_start), str(intron_end), '.', columns[6], '.', f"Parent={gene_id}"])
                        output_handle.write(intron_line + '\n')

if __name__ == "__main__":
    gff3_file = "/path/to/Annotation_file.gff3"
    output_file = "/path/to/Annotation_file_with_intron.gff3"
    add_introns_to_gff3(gff3_file, output_file)
```

# Extracting Sequences
Databases typically provide FASTA files for each species, containing gene, transcript, and/or CDS sequences. However, if such files are missing or an intron FASTA file is needed, this script uses the coordinates from the GFF3 file to extract the sequences and create a corresponding FASTA file.

```python
from Bio import SeqIO

def extract_gene_sequences(genome_file, gff3_file, output_file):
    genome = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))

    sequences = []

    with open(gff3_file, "r") as gff3_handle:
        for line in gff3_handle:
            if line.startswith("#"):
                continue

            columns = line.strip().split("\t")
            feature_type = columns[2]
            start = int(columns[3]) - 1  # Convert to 0-based
            end = int(columns[4])
            attributes = dict(item.split("=") for item in columns[8].split(";"))

            seq_id = columns[0]
            seq = genome[seq_id].seq[start:end]

            header = f">{attributes['ID']}_{feature_type}_{start + 1}:{end}"
            sequences.append(f"{header}\n{seq}")

    with open(output_file, "w") as output_handle:
        output_handle.write("\n".join(sequences))

if __name__ == "__main__":
    genome_file = "TAIR10_chr_all.fasta"
    gff3_file = "/path/to/Annotation_file_with_intron.gff3"
    output_file = "/path/to/all_sequences.fasta"
    extract_gene_sequences(genome_file, gff3_file, output_file)
```

# Domain Search
Domains are fundamental units of protein structure and function. They are essential for understanding protein interactions, functions, and evolutionary history. This script extracts the number of matches for specified domains in each sequence provided in the FASTA file. When used in conjunction with the sequence extraction script, it can identify and locate domains of interest throughout the genome.

```python
from Bio import SeqIO

def find_motif(sequence, motif):
    count = 0
    length = len(motif)
    for i in range(len(sequence) - length + 1):
        if sequence[i:i+length] == motif:
            count += 1
    return count

def main(fasta_file, motif, output_file):
    with open(fasta_file, "r") as handle, open(output_file, "w") as outfile:
        for record in SeqIO.parse(handle, "fasta"):
            sequence = str(record.seq)
            motif_count = find_motif(sequence, motif)
            result_line = f"Sequence: {record.id}, Motif Count: {motif_count}\n"
            print(result_line.strip())
            outfile.write(result_line)

if __name__ == "__main__":
    fasta_file = "/path/to/full_sequences.fasta"
    motif_to_find = "AGAAG"
    output_file = "/path/to/motif_counts.txt"
    main(fasta_file, motif_to_find, output_file)

```
