#-- import modules --#
import io
import sys
import argparse
import re
from Bio import SeqIO

#-- functions --#
def extract_exons_and_positions(seq: str) -> tuple:
    """
    Extract exons and their positions from a sequence.
    Parameters:
        -- seq: the DNA sequence
    Returns:
        --tuple: (full cDNA sequence, skipped cDNA sequence, list of exon positions)
    """
    exons = []
    exons_pos = []

    for match in re.finditer(r'[A-Z]+', str(seq)):
        exons.append(match.group())
        start = match.start() + 1  # 1-based
        end = match.end()
        exons_pos.append((start, end))

    full_exon_seq = ''.join(exons)

    if len(exons) >= 3:
        mid = len(exons) // 2
        skip_exon_seq = ''.join(exons[:mid] + exons[mid+1:])

    return full_exon_seq, skip_exon_seq, exons_pos

#-- main execution --#
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Create a reference with exon only (may include barcode).", allow_abbrev = False)
    parser.add_argument("-r", "--reference", type = str, required = True, help = "Reference FASTA file")
    parser.add_argument("-l", "--lib_type",  type = str, required = True, help = "Library type (random_intron, random_exon, muta_exon)")
    parser.add_argument("-p", "--prefix",    type = str, required = True, help = "Output prefix")
    
    args, unknown = parser.parse_known_args()

    if unknown:
        print(f"Error: Unrecognized arguments: {' '.join(unknown)}", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    valid_lib_types = {"random_intron", "random_exon", "muta_exon"}
    if args.lib_type not in valid_lib_types:
        sys.exit(f"Error: Invalid --lib_type '{args.lib_type}'. Must be one of: {', '.join(valid_lib_types)}")

    output_fasta = f"{args.prefix}.exon_ref.fasta"
    output_positions = f"{args.prefix}.exon_pos.tsv"
    
    with open(output_fasta, "w") as fasta_out, open(output_positions, "w") as tsv_out:
        fasta_records = SeqIO.parse(args.reference, "fasta")

        if args.lib_type == "random_intron":
            record = next(fasta_records)
            full_exon_seq, skip_exon_seq, exons_pos = extract_exons_and_positions(record.seq)
            fasta_out.write(f">exon_inclusion_random\n{full_exon_seq}\n")
            fasta_out.write(f">exon_skipping_random\n{skip_exon_seq}\n")
            for i, (start, end) in enumerate(exons_pos, 1):
                tsv_out.write(f"random\tE{i}\t{start}\t{end}\n")
        else:
            for record in SeqIO.parse(args.reference, "fasta"):
                full_exon_seq, skip_exon_seq, exons_pos = extract_exons_and_positions(record.seq)
                fasta_out.write(f">exon_inclusion_{record.id}\n{full_exon_seq}\n")
                fasta_out.write(f">exon_skipping_{record.id}\n{skip_exon_seq}\n")
                for i, (start, end) in enumerate(exons_pos, 1):
                    tsv_out.write(f"{record.id}\tE{i}\t{start}\t{end}\n")
