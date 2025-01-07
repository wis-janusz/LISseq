"""A CLI utility to remove 5'LTR and polyA tails from LISseq reads.

Reads a fastq file contating raw reads, removes HIV 5'LTR and polyA sequences
and saves as a new fastq file.

Typical usage example:
  python LISseq_clean.py input_file.fastq output_file.fastq
"""

import argparse
import sys
import regex
from Bio import SeqIO

def parse_args(arg_list: list[str] | None):
  parser = argparse.ArgumentParser()
  parser.add_argument("input_file", type=str, help="Path to input fastq file.")
  parser.add_argument("output_file", type=str, help="Path to output fastq file.")
  parser.add_argument("-ltr", type=str, default="GGAGTGAATTAGCCCTTCCA")
  return parser.parse_args(arg_list)


def main(arg_list: list[str] | None = None):
    args = parse_args(arg_list)
    in_seq = SeqIO.parse(args.input_file, format="fastq")
    out_seq = []
    for seq in in_seq:
       if "N" not in seq.seq:
        new_seq = seq[len(args.ltr):]
        polyA_pos = new_seq.seq.find("AAAAA")
        if polyA_pos != -1:
           new_seq = new_seq[:polyA_pos]

        out_seq.append(new_seq)
    SeqIO.write(out_seq, args.output_file, format="fastq")


if __name__ == "__main__":
    main()