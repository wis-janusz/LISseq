"""A CLI utility to map LISseq reads.

Reads fq.gz files contating raw reads from the input directory, removes HIV 5'LTR 
and polyA sequences, filters by quality and length of remaining sequence
then maps reads to reference genome using bowtie2
and finally outputs integration sites.

Typical usage example:
  python LISseq_clean.py input_directory output_directory
"""

import argparse
import pathlib
import gzip
import shutil
import subprocess
import pysam
import requests
import os
import pandas as pd
from Bio import SeqIO, SeqRecord, Seq

ENSEMBL_SERVER = "http://rest.ensembl.org"


def _parse_args(arg_list: list[str] | None):
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", type=str, help="Path to input directory.")
    parser.add_argument("output_dir", type=str, help="Path to output directory.")
    parser.add_argument("-ltr", type=str, default="GGAGTGAATTAGCCCTTCCA", help="Sequence of the part of HIV 5' LTR present in raw reads.")
    parser.add_argument("-save_all_loci", action="store_true", help="If True, will save a table with all mapped loci for each file.")
    parser.add_argument("-A", type=int, default=5, help="Minimum length of streches of A to be removed from reads.")
    parser.add_argument("-q", type=int, default=20, help="Minimum sequencing quality of reads and minimum alingment quality of mapping.")
    parser.add_argument("-l", type=int, default=20, help="Minumim length of cleaned reads to be considered for mapping.")
    parser.add_argument("-idx", type=str, default="GRCh38_noalt_as")
    return parser.parse_args(arg_list)


def _find_fqgz(dir: str) -> list[pathlib.Path]:
    in_path = pathlib.Path(dir)
    if in_path.exists() == False:
        raise FileNotFoundError("Please provide a correct input directory.")
    if in_path.is_dir() == False:
        raise NotADirectoryError("Please provide a correct input directory.")
    
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
    if ltr in read.seq:
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


def _save_clean_reads(read_list: list, out_dir, filename):
    out_path = pathlib.Path(out_dir)
    if out_path.exists() == False:
        raise FileNotFoundError("Please provide a correct output directory.")
    if out_path.is_dir() == False:
        raise NotADirectoryError("Please provide a correct output directory.")
    
    SeqIO.write(read_list, f"{out_dir}/{filename}_clean.fq", format="fastq")


def _cleanup_reads(args):
    IS_dict = {}
    in_fqgzs = _find_fqgz(args.input_dir)
    print(f"Found {len(in_fqgzs)} fq.gz files in {args.input_dir}")
    counter = 1
    for infile in in_fqgzs:
        print(f"Processing file {counter} of {len(in_fqgzs)}.", end="\r", flush=True)
        sample_name = infile.stem.split("_")[0]
        raw_reads = _parse_fqgz(infile)
        nice_reads = []
        raw_reads_counter = 0
        for read in raw_reads:
            raw_reads_counter += 1
            nice_read = _clean_read(read, args.ltr, args.A, args.q, args.l)
            if nice_read.id != "":
                nice_reads.append(nice_read)
        if len(nice_reads) > 0:
            _save_clean_reads(nice_reads, args.output_dir, infile.stem)
        if raw_reads_counter > 0:
            IS_dict[sample_name] = {"raw_reads":raw_reads_counter, "filtered_reads":len(nice_reads)}
        counter += 1
    return IS_dict


def _find_fq(dir: str) -> list[pathlib.Path]:
    in_path = pathlib.Path(dir)
    return list(in_path.glob("**/*.fq"))


def _check_bowtie2(path:str = None):
    if shutil.which("bowtie2", path=path) != None:
        return True
    else:
        raise FileNotFoundError("bowtie2 not found. Please install it and add it to PATH.")


def _check_genome_index(idx_dir):
    idx_list = list(pathlib.Path(idx_dir).glob("*.bt2"))
    if len(idx_list) > 0:
        return idx_list[0].name.split(".")[0]
    else:
        raise FileNotFoundError(f"No bowtie2-indexed genome found in {idx_dir}, please download or build one.")


def _bowtie_map(clean_fq_file:pathlib.Path, idx_dir:str, idx_name:str, out_dir:str):
    out_file = f"{out_dir}/{clean_fq_file.stem}.sam"
    map_cmd = [
        "bowtie2",
        "-x",
        f"{idx_dir}/{idx_name}",
        "-U",
        "-q",
        clean_fq_file,
        "-S",
        out_file,
    ]
    result = subprocess.run(map_cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise Exception(result.stderr)
    return out_file

def _read_sam_to_df(sam_file:str) -> pd.DataFrame:
    reads = pysam.AlignmentFile(sam_file, mode="r")
    read_dict = {"read":[], "flag":[],"chr":[], "pos":[],"Q":[], "seq":[]}
    for read in reads:
        read_dict["read"].append(read.query_name)
        read_dict["flag"].append(read.flag)
        read_dict["chr"].append(read.reference_name)
        read_dict["pos"].append(read.reference_start)
        read_dict["Q"].append(read.mapping_quality)
        read_dict["seq"].append(read.query_sequence)
    mappings_df = pd.DataFrame(read_dict).set_index("read")
    return mappings_df


def _extract_IS(mappings_df:pd.DataFrame, Q) -> pd.DataFrame:
    filtered_df = mappings_df[(mappings_df["Q"] >= Q)]
    loci = filtered_df.groupby(["chr","pos","flag"]).agg(depth = ('Q', "count"), mean_q = ("Q", "mean"), seq = ("seq", "first"))
    loci = loci[loci["depth"] >= 1000]
    return loci, len(filtered_df)


def _map_IS(args, IS_dict):
    in_fqs = _find_fq(args.output_dir)
    flag_to_strand = {16:"-",0:"+"}
    if _check_bowtie2() == True:
        idx_name = _check_genome_index(args.idx)
        counter = 1
        for clean_fq in in_fqs:
            sample_name = clean_fq.stem.split("_")[0]
            print(f"Processing file {counter} of {len(in_fqs)}.", end="\r", flush=True)
            alignments = _bowtie_map(clean_fq, args.idx, idx_name, args.output_dir)
            mappings_df = _read_sam_to_df(alignments)
            loci, total_mapped = _extract_IS(mappings_df, args.q)
            loci = loci.sort_values("depth", ascending=False)
            if args.save_all_loci == True:
                loci.to_csv(f"{args.output_dir}/{clean_fq.stem}.csv")
            best_is = loci[loci["depth"] == loci["depth"].max()].reset_index()
            if len(best_is) > 0:
                IS_dict[sample_name]["chr"] = best_is.iloc[0,0]
                IS_dict[sample_name]["pos"] = best_is.iloc[0,1]
                IS_dict[sample_name]["strand"] = flag_to_strand[best_is.iloc[0,2]]  
                IS_dict[sample_name]["total_mapped"] = total_mapped
                IS_dict[sample_name]["mapped_to_locus"] = best_is.iloc[0,3]
                IS_dict[sample_name]["mapping_quality"] = round(best_is.iloc[0,4],2)
                IS_dict[sample_name]["seq"] = best_is.iloc[0,5]
            else:
                IS_dict[sample_name]["chr"] = 0
                IS_dict[sample_name]["pos"] = 0
                IS_dict[sample_name]["strand"] = ""     
                IS_dict[sample_name]["total_mapped"] = total_mapped
                IS_dict[sample_name]["mapped_to_locus"] = 0
                IS_dict[sample_name]["mapping_quality"] = 0   
                IS_dict[sample_name]["seq"] = ""
            counter += 1

def _get_gene(chr:int, pos:int) -> str:
    server = ENSEMBL_SERVER
    overlap_endpoint = f"/overlap/region/human/{chr}:{pos}-{pos}?feature=gene"

    overlap_response = requests.get(server+overlap_endpoint, headers={"Content-Type":"application/json"})
    if overlap_response.status_code == 200:
        genes = overlap_response.json()
        if len(genes) > 0:
            gene_data = genes[0]
            gene_symbol = gene_data.get("external_name", "id")
            return gene_symbol
        else:
            return "intergenic"
    else:
        return None


def main(arg_list: list[str] | None = None):
    args = _parse_args(arg_list)
    if pathlib.Path(args.output_dir).exists() == False:
        os.makedirs(args.output_dir)
    IS_dict = _cleanup_reads(args)
    print("\nFinished cleaning up reads.")
    print("Mapping and extracting IS.")
    _map_IS(args, IS_dict)
    IS_df = pd.DataFrame.from_dict(IS_dict).T
    IS_df["gene"] = IS_df.apply(lambda row: _get_gene(row["chr"], row["pos"]), axis=1)
    IS_df = IS_df[["chr", "pos","strand","gene","raw_reads","filtered_reads","total_mapped","mapped_to_locus","mapping_quality","seq"]]
    IS_df.to_csv(f"{args.output_dir}/integration_sites.csv")
    print("\nDone.")


if __name__ == "__main__":
    main()
