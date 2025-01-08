"""A CLI utility to remove 5'LTR and polyA tails from LISseq reads.

Reads a fastq file contating raw reads, removes HIV 5'LTR and polyA sequences
and saves as a new fastq file.

Typical usage example:
  python LISseq_clean.py input_file.fastq output_file.fastq
"""

import argparse
import pathlib
import gzip
from Bio import SeqIO, SeqRecord, Seq


def _parse_args(arg_list: list[str] | None):
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", type=str, help="Path to input directory.")
    parser.add_argument("output_dir", type=str, help="Path to output directory.")
    parser.add_argument("-ltr", type=str, default="GGAGTGAATTAGCCCTTCCA")
    parser.add_argument("-A", type=int, default=5)
    parser.add_argument("-q", type=int, default=20)
    return parser.parse_args(arg_list)


def _find_fqgz(dir: str):
    in_path = pathlib.Path(dir)
    return in_path.glob("*.fq.gz")


def _parse_fqgz(input_file_path: pathlib.Path):
    with gzip.open(input_file_path, mode="rt") as input_file_handle:
        in_seq = SeqIO.parse(input_file_handle, format="fastq")
        try:
            for read in in_seq:
                yield read
        except ValueError:
            return


def _clean_read(
    read: SeqRecord.SeqRecord, ltr: str, polyAlen: int, min_quality: int
) -> SeqRecord.SeqRecord:
    new_read = read[len(ltr) :]
    polyA_pos = new_read.seq.find(polyAlen * "A")
    if polyA_pos != -1:
        new_read = new_read[:polyA_pos]
    if (
        "N" not in new_read.seq
        and min(new_read.letter_annotations["phred_quality"]) >= min_quality
    ):
        return new_read
    else:
        return SeqRecord.SeqRecord(seq=Seq.Seq(""), id="")


def _save_clean_reads(read_list: list, out_dir):
    SeqIO.write(read_list, out_dir, format="fastq")


def main(arg_list: list[str] | None = None):
    args = _parse_args(arg_list)
    in_fqgzs = _find_fqgz(args.input_dir)
    for infile in in_fqgzs:
        raw_reads = _parse_fqgz(infile)
        nice_reads = []
        for read in raw_reads:
            nice_read = _clean_read(read, args.ltr, args.A, args.q)
            if nice_read.id != "":
                nice_reads.append(nice_read)
        if len(nice_reads) > 0:
            _save_clean_reads(nice_reads, f"{args.output_dir}/{infile.stem}_clean.fq")


if __name__ == "__main__":
    main()
