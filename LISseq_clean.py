"""A CLI utility to clean up LISseq reads.

Reads fq.gz files contating raw reads from the input folder, removes HIV 5'LTR 
and polyA sequences, filters by quality and length of remaining sequence 
and saves as new fq files in output folder.

Typical usage example:
  python LISseq_clean.py input_folder output_folder
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
    parser.add_argument("-l", type=int, default=15)
    return parser.parse_args(arg_list)


def _find_fqgz(dir: str) -> list:
    in_path = pathlib.Path(dir)
    return list(in_path.glob("**/*.fq.gz"))


def _parse_fqgz(input_file_path: pathlib.Path):
    with gzip.open(input_file_path, mode="rt") as input_file_handle:
        in_seq = SeqIO.parse(input_file_handle, format="fastq")
        try:
            for read in in_seq:
                yield read
        except ValueError:
            return


def _clean_read(
    read: SeqRecord.SeqRecord, ltr: str, polyAlen: int, min_quality: int, min_len: int
) -> SeqRecord.SeqRecord:
    new_read = read[len(ltr) :]
    polyA_pos = new_read.seq.find(polyAlen * "A")
    if polyA_pos != -1:
        new_read = new_read[:polyA_pos]
    if len(new_read.seq) >= min_len:
        if (
            "N" not in new_read.seq
            and min(new_read.letter_annotations["phred_quality"]) >= min_quality
        ):
            return new_read
    return SeqRecord.SeqRecord(seq=Seq.Seq(""), id="")


def _save_clean_reads(read_list: list, out_dir):
    SeqIO.write(read_list, out_dir, format="fastq")


def _cleanup(args):
    in_fqgzs = _find_fqgz(args.input_dir)
    print(f"Found {len(in_fqgzs)} .fq.gz files in {args.input_dir}")
    counter = 1
    for infile in in_fqgzs:
        print(f"Processing file {counter} of {len(in_fqgzs)}.", end="\r", flush=True)
        raw_reads = _parse_fqgz(infile)
        nice_reads = []
        for read in raw_reads:
            nice_read = _clean_read(read, args.ltr, args.A, args.q, args.l)
            if nice_read.id != "":
                nice_reads.append(nice_read)
        if len(nice_reads) > 0:
            _save_clean_reads(nice_reads, f"{args.output_dir}/{infile.stem}_clean.fq")
        counter += 1


def main(arg_list: list[str] | None = None):
    args = _parse_args(arg_list)
    _cleanup(args)


if __name__ == "__main__":
    main()
