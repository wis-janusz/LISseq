"""A CLI utility to map LISseq reads.

Reads fq.gz files contating raw reads from the input directory, removes HIV 5'LTR
and polyA sequences, filters by quality and length of remaining sequence
then maps reads to reference genome using bowtie2
and finally outputs integration sites.

Typical usage example:
  python LISseq_clean.py input_directory output_directory
"""

import argparse
import gzip
import os
import pathlib
import shutil
import subprocess

import pandas as pd
import pysam
import requests
from Bio import SeqIO, SeqRecord, Seq

ENSEMBL_SERVER = "http://rest.ensembl.org"


def parse_args(arg_list: list[str] | None):
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", type=str, help="Path to input directory.")
    parser.add_argument("output_dir", type=str, help="Path to output directory.")
    parser.add_argument(
        "--save_all_loci",
        action="store_true",
        help="If True, will save a table with all mapped loci for each file. Default is False.",
    )
    parser.add_argument(
        "--idx",
        type=str,
        default="GRCh38_noalt_as",
        help="Specifies the directory containing reference genome index. Default is 'GRCh38_noalt_as'.",
    )
    parser.add_argument(
        "--ltr",
        type=str,
        default="GGAGTGAATTAGCCCTTCCA",
        help="Sequence of the part of HIV 5' LTR present in raw reads. Default is 'GGAGTGAATTAGCCCTTCCA'.",
    )
    parser.add_argument(
        "--ltrmax",
        type=int,
        default=1,
        help="Maximum number of mismatches allowed when looking for the LTR sequence. Default is 1.",
    )
    parser.add_argument(
        "-A",
        type=int,
        default=5,
        help="Minimum length of stretches of A to be removed from reads. Default is 5.",
    )
    parser.add_argument(
        "-q",
        type=int,
        default=20,
        help="Minimum sequencing quality of reads and minimum alignment quality of mapping. Default is 20.",
    )
    parser.add_argument(
        "-l",
        type=int,
        default=20,
        help="Minimum length of cleaned reads to be considered for mapping. Default is 20.",
    )
    return parser.parse_args(arg_list)


def _find_fqgz(dir: str) -> list[pathlib.Path]:
    in_path = pathlib.Path(dir)
    if not in_path.exists():
        raise FileNotFoundError("Please provide a correct input directory.")
    if not in_path.is_dir():
        raise NotADirectoryError("Please provide a correct input directory.")

    return list(in_path.glob("**/*.fq.gz"))


def _parse_fqgz(input_file_path: pathlib.Path):
    with gzip.open(input_file_path, mode="rt") as input_file_handle:
        in_seq = SeqIO.parse(input_file_handle, format="fastq")
        try:
            for read in in_seq:
                yield read
        except ValueError:
            # If the input fq.gz file is empty or has wrong format, just ignore it and return None
            return


def _has_ltr(read_seq: Seq.Seq, ltr: str, max_mismatches: int = 1) -> bool:
    """
    Check if a given sequence contains a specified long terminal repeat (LTR)
    with a maximum number of allowed mismatches.

    Args:
        read_seq (Seq.Seq): The sequence to be checked.
        ltr (str): The long terminal repeat sequence to search for.
        max_mismatches (int, optional): The maximum number of mismatches allowed. Defaults to 1.

    Returns:
        bool: True if the LTR is found within the allowed mismatches, False otherwise.
    """
    for i in range(len(read_seq) - len(ltr) + 1):
        window = read_seq[i : i + len(ltr)]
        mismatches = sum(1 for a, b in zip(window, ltr) if a != b)
        if mismatches <= max_mismatches:
            return True
    return False


def _clean_read(
    read: SeqRecord.SeqRecord,
    ltr: str,
    polyAlen: int,
    min_quality: int,
    min_len: int,
    max_mismatches: int = 1,
) -> SeqRecord.SeqRecord:
    """
    Cleans a sequencing read by removing specified LTR sequences and poly-A tails,
    and filtering based on quality and length.

    Args:
        read (SeqRecord.SeqRecord): The sequencing read to be cleaned.
        ltr (str): The LTR sequence to be removed from the start of the read.
        polyAlen (int): The minimum length of the poly-A tail to be removed.
        min_quality (int): The minimum quality score of any position in the read.
        min_len (int): The minimum length required for the read after cleaning.
        max_mismatches (int, optional): The maximum number of mismatches allowed in the LTR sequence. Defaults to 1.

    Returns:
        SeqRecord.SeqRecord: The cleaned sequencing read. If the read does not meet the criteria, an empty SeqRecord is returned.
    """
    if _has_ltr(read.seq, ltr, max_mismatches):
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


def _save_clean_reads(read_list: list, out_dir: str, filename: str):
    out_path = pathlib.Path(out_dir)
    if out_path.exists() == False:
        raise FileNotFoundError("Please provide a correct output directory.")
    if out_path.is_dir() == False:
        raise NotADirectoryError("Please provide a correct output directory.")

    SeqIO.write(read_list, f"{out_dir}/temp/{filename}_clean.fq", format="fastq")


def cleanup_reads(args: argparse.Namespace) -> dict:
    """
    Cleans up reads from input fastq.gz files and saves the cleaned reads to the output directory.
    Ignores files with 0 raw reads, but not with 0 reads passing filtering criteria.

    Args:
        args: An argparse.Namespace object containing the following attributes:
            - input_dir (str): Directory containing input fastq.gz files.
            - output_dir (str): Directory to save cleaned reads.
            - ltr (str): Long terminal repeat sequence to clean from reads.
            - A (int): Minimum number of 'A' bases required in the read.
            - q (int): Minimum quality score required for the read.
            - l (int): Minimum length required for the read.
            - ltrmax (int): Maximum length of the long terminal repeat sequence.

    Returns:
        dict: A dictionary of raw and filtered reads number for each sample; values are dictionaries with the following keys:
            - "raw_reads" (int): Number of raw reads processed.
            - "filtered_reads" (int): Number of reads that passed the filtering criteria.
    """
    read_no_dict = {}
    in_fqgzs = _find_fqgz(args.input_dir)
    print(f"Found {len(in_fqgzs)} fq.gz files in {args.input_dir}.")
    counter = 1
    for infile in in_fqgzs:
        print(f"Processing file {counter} of {len(in_fqgzs)}.", end="\r", flush=True)
        sample_name = infile.stem.split("_")[0]
        raw_reads = _parse_fqgz(infile)
        nice_reads = []
        raw_reads_counter = 0
        for read in raw_reads:
            raw_reads_counter += 1
            nice_read = _clean_read(read, args.ltr, args.A, args.q, args.l, args.ltrmax)
            if nice_read.id != "":
                nice_reads.append(nice_read)
        if raw_reads_counter > 0:
            _save_clean_reads(nice_reads, args.output_dir, infile.stem)
            read_no_dict[sample_name] = {
                "raw_reads": raw_reads_counter,
                "filtered_reads": len(nice_reads),
            }
        else:
            print(f"File {infile} contains no reads!")
        counter += 1
    return read_no_dict


def _find_fq(dir: str) -> list[pathlib.Path]:
    in_path = pathlib.Path(f"{dir}/temp/")
    return list(in_path.glob("**/*.fq"))


def _check_bowtie2(path: str = None):
    if shutil.which("bowtie2", path=path) != None:
        return True
    else:
        raise FileNotFoundError(
            "bowtie2 not found. Please install it and add it to PATH."
        )


def _check_genome_index(idx_dir):
    idx_list = list(pathlib.Path(idx_dir).glob("*.bt2"))
    if len(idx_list) > 0:
        return idx_list[0].name.split(".")[0]
    else:
        raise FileNotFoundError(
            f"No bowtie2-indexed genome found in {idx_dir}, please download or build one."
        )


def _bowtie_map(clean_fq_file: pathlib.Path, idx_dir: str, idx_name: str, out_dir: str):
    out_file = os.path.join(out_dir, "temp", f"{clean_fq_file.stem}.sam")
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


def _read_sam_to_df(sam_file: str) -> pd.DataFrame:
    """
    Reads a SAM file and converts it into a pandas DataFrame.

    Parameters:
    sam_file (str): The path to the SAM file to be read.

    Returns:
    pd.DataFrame: A DataFrame containing the SAM file data with columns:
        - 'read': The query name of the read.
        - 'flag': The bitwise flag.
        - 'chr': The reference sequence name.
        - 'pos': The 1-based leftmost mapping position.
        - 'Q': The mapping quality.
        - 'seq': The query sequence with a 5'LTR tag.

    Raises:
    ValueError: If the file is not a valid SAM alignment file.
    """
    try:
        with pysam.AlignmentFile(sam_file, mode="r") as reads:
            read_dict = {
                "read": [],
                "flag": [],
                "chr": [],
                "pos": [],
                "Q": [],
                "seq": [],
            }
            for read in reads:
                read_dict["read"].append(read.query_name)
                read_dict["flag"].append(read.flag)
                read_dict["chr"].append(read.reference_name)

                if read.flag == 16:
                    read_dict["pos"].append(
                        read.reference_start + len(read.query_sequence)
                    )
                    read_dict["seq"].append(read.query_sequence + "-5'LTR")
                else:
                    read_dict["pos"].append(read.reference_start + 1)
                    read_dict["seq"].append("5'LTR-" + read.query_sequence)
                read_dict["Q"].append(read.mapping_quality)

            mappings_df = pd.DataFrame(read_dict).set_index("read")
        return mappings_df
    except ValueError:
        with open(sam_file, mode="r") as in_file:
            print(in_file.read())
        raise ValueError("File is not a valid SAM alignment file.")


def _extract_IS(mappings_df: pd.DataFrame, Q) -> tuple[pd.DataFrame, int]:
    """
    Extract integration sites (IS) from a DataFrame of mappings based on a quality threshold.

    This function filters the input DataFrame to include only rows where the quality score (Q)
    is greater than or equal to the specified threshold. It then groups the filtered data by
    chromosome, position, and flag, and calculates the depth (count of reads), mean quality score,
    and the first sequence for each group. Finally, it filters the grouped data to include only
    groups with a depth of at least 1000.

    Args:
        mappings_df (pd.DataFrame): A DataFrame containing mapping information with columns
                                    "chr" (chromosome), "pos" (position), "flag", "Q" (quality score),
                                    and "seq" (sequence).
        Q (int): The quality score threshold for filtering the mappings.

    Returns:
        tuple[pd.DataFrame, int]: A tuple containing:
            - A DataFrame with the grouped and filtered integration sites, including columns
              "depth" (count of reads), "mean_q" (mean quality score), and "seq" (first sequence).
            - An integer representing the total number of reads mapped to the reference genome.
    """
    filtered_df = mappings_df[(mappings_df["Q"] >= Q)]
    loci = filtered_df.groupby(["chr", "pos", "flag"]).agg(
        depth=("Q", "count"), mean_q=("Q", "mean"), seq=("seq", "first")
    )
    loci = loci[loci["depth"] >= 1000]
    return loci, len(filtered_df)


def map_IS(args: argparse.Namespace, read_no_dict: dict) -> dict:
    """
    Maps integration sites (IS) from cleaned sequencing reads and returns a dictionary with mapping details.

    Args:
        args (argparse.Namespace): Command-line arguments containing parameters such as output directory, index, and quality threshold.
        read_no_dict (dict): Dictionary containing the number of raw and filtered reads for each sample.

    Returns:
        dict: A dictionary where each key is a sample name and each value is another dictionary containing mapping details for each integration site,
                in descending order of sequencing depth (count of reads mapped to the site). The mapping details include chromosome number, position,
                strand, no. of raw reads, no. of filtered reads, no. of total mapped reads, no. ofreads mapped to the locus, mean mapping quality,
                and sequence.
    """
    IS_dict = {}
    in_fqs = _find_fq(args.output_dir)
    flag_to_strand = {16: "+", 0: "-"}
    if _check_bowtie2() == True:
        idx_name = _check_genome_index(args.idx)
        counter = 1
        for clean_fq in in_fqs:
            sample_name = clean_fq.stem.split("_")[0]
            print(f"Processing file {counter} of {len(in_fqs)}.", end="\r", flush=True)
            alignments = _bowtie_map(clean_fq, args.idx, idx_name, args.output_dir)
            mappings_df = _read_sam_to_df(alignments)
            loci, total_mapped = _extract_IS(mappings_df, args.q)
            loci = loci.sort_values("depth", ascending=False).reset_index()
            IS_dict[sample_name] = {}
            if len(loci) > 0:
                for i, row in loci.iterrows():
                    IS_dict[sample_name][i] = {}
                    IS_dict[sample_name][i]["chr"] = row["chr"]
                    IS_dict[sample_name][i]["pos"] = row["pos"]
                    IS_dict[sample_name][i]["strand"] = flag_to_strand[row["flag"]]
                    IS_dict[sample_name][i]["raw_reads"] = read_no_dict[sample_name][
                        "raw_reads"
                    ]
                    IS_dict[sample_name][i]["filtered_reads"] = read_no_dict[
                        sample_name
                    ]["filtered_reads"]
                    IS_dict[sample_name][i]["total_mapped"] = total_mapped
                    IS_dict[sample_name][i]["mapped_to_locus"] = row["depth"]
                    IS_dict[sample_name][i]["mapping_quality"] = round(row["mean_q"], 2)
                    IS_dict[sample_name][i]["seq"] = row["seq"]
            else:
                IS_dict[sample_name][0] = {}
                IS_dict[sample_name][0]["chr"] = "none"
                IS_dict[sample_name][0]["pos"] = 0
                IS_dict[sample_name][0]["strand"] = "none"
                IS_dict[sample_name][0]["raw_reads"] = read_no_dict[sample_name][
                    "raw_reads"
                ]
                IS_dict[sample_name][0]["filtered_reads"] = read_no_dict[sample_name][
                    "filtered_reads"
                ]
                IS_dict[sample_name][0]["total_mapped"] = total_mapped
                IS_dict[sample_name][0]["mapped_to_locus"] = 0
                IS_dict[sample_name][0]["mapping_quality"] = 0.0
                IS_dict[sample_name][0]["seq"] = "none"
            counter += 1
    return IS_dict


def _get_gene(chr: int, pos: int) -> str:
    """
    Retrieve the gene symbol at a specific chromosome position using the Ensembl API.

    Args:
        chr (int): Chromosome number.
        pos (int): Position on the chromosome.

    Returns:
        str: Gene symbol if found, "intergenic" if no gene is found at the position, 
             or None if the API request fails.
    """
    server = ENSEMBL_SERVER
    overlap_endpoint = f"/overlap/region/human/{chr}:{pos}-{pos}?feature=gene"

    overlap_response = requests.get(
        server + overlap_endpoint, headers={"Content-Type": "application/json"}
    )
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


def format_data_frame(IS_dict: dict) -> pd.DataFrame:
    """
    Converts a dictionary of integration sites (IS) into a formatted pandas DataFrame.

    Args:
        IS_dict (dict): A dictionary where keys are integration site identifiers and values are dictionaries
                        containing integration site data.

    Returns:
        pd.DataFrame: A pandas DataFrame with the integration site data, including an additional 'gene' column
                      determined by the _get_gene function. The DataFrame columns are ordered and typed as follows:
                      - 'chr': object
                      - 'pos': int
                      - 'strand': object
                      - 'gene': object
                      - 'raw_reads': int
                      - 'filtered_reads': int
                      - 'total_mapped': int
                      - 'mapped_to_locus': int
                      - 'mapping_quality': float
                      - 'seq': object
    """
    IS_df = pd.DataFrame.from_dict(IS_dict).T
    IS_df["gene"] = IS_df.apply(lambda row: _get_gene(row["chr"], row["pos"]), axis=1)
    IS_df = IS_df[
        [
            "chr",
            "pos",
            "strand",
            "gene",
            "raw_reads",
            "filtered_reads",
            "total_mapped",
            "mapped_to_locus",
            "mapping_quality",
            "seq",
        ]
    ]
    IS_df = IS_df.astype(
        {
            "pos": "int",
            "raw_reads": "int",
            "filtered_reads": "int",
            "total_mapped": "int",
            "mapped_to_locus": "int",
            "mapping_quality": "float",
        }
    )
    return IS_df


def main(arg_list: list[str] | None = None):
    args = parse_args(arg_list)
    if pathlib.Path(args.output_dir).exists() == False:
        os.makedirs(f"{args.output_dir}/temp/")
    read_no_dict = cleanup_reads(args)
    print("\nFinished cleaning up reads.")
    print("Mapping and extracting IS.")
    IS_dict = map_IS(args, read_no_dict)
    best_loci_dict = {sample: loci[0] for sample, loci in IS_dict.items()}
    if args.save_all_loci == True:
        for sample in IS_dict.keys():
            all_IS_df = format_data_frame(IS_dict[sample])
            all_IS_df.to_csv(
                f"{args.output_dir}/{sample}_all_loci.csv", index_label="loci_no"
            )
    best_IS_df = format_data_frame(best_loci_dict)
    best_IS_df.to_csv(f"{args.output_dir}/integration_sites.csv", index_label="sample")
    print("\nDone.")


if __name__ == "__main__":
    main()
