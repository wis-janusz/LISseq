import pandas as pd
from Bio import SeqIO, SeqRecord
import LISseq

def test_find_fqgz(test_input_dir):
    names = set(filename.stem for filename in LISseq._find_fqgz(test_input_dir))
    assert names == {"0_cos.fq","1_cos.fq","2_cos.fq","3_cos.fq"}

def test_correct_fqgz(correct_fqgz):
    test_seq = LISseq._parse_fqgz(correct_fqgz)
    for read in test_seq:
        assert type(read) == SeqRecord.SeqRecord

def test_incorrect_fqgz(incorrect_fqgz):
    test_seq = LISseq._parse_fqgz(incorrect_fqgz)
    for read in test_seq:
        assert read == None

def test_cleanup(default_test_args, expected_output_path, tmp_path, expected_IS_dict_raw_reads):
    args = LISseq._parse_args(default_test_args)
    test_dict = LISseq._cleanup_reads(args)
    for i in range(1,3):
        with open(f"{tmp_path}/{i}_cos.fq_clean.fq", mode="r") as test_file, open(f"{expected_output_path}/{i}_cos.fq_clean.fq", mode="r") as exp_file:
            test_lines = test_file.readlines()
            exp_lines = exp_file.readlines()
            assert test_lines == exp_lines
    assert test_dict == expected_IS_dict_raw_reads

def test_find_fq(expected_output_path):
    names = set(filename.stem for filename in LISseq._find_fq(expected_output_path))
    assert names == {"1_cos.fq_clean", "2_cos.fq_clean"}

def test_check_bowtie_exists():
    assert LISseq._check_bowtie2() != None

def test_check_bowtie_not_exists(syspath_no_bowtie):
    assert LISseq._check_bowtie2(syspath_no_bowtie) == None

def test_check_hg38idx_default(hg38idx_default_dir):
    assert LISseq._check_hg38_index(hg38idx_default_dir) == "GRCh38_noalt_as"

def test_check_hg38idx_not_exists(hg38idx_nonexistant_dir):
    assert LISseq._check_hg38_index(hg38idx_nonexistant_dir) == ""

def test_check_hg38idx_wrong(hg38idx_wrong_dir):
    assert LISseq._check_hg38_index(hg38idx_wrong_dir) == ""

def test_read_sam(expected_output_path):
    test_df1 = LISseq._read_sam_to_df(f"{expected_output_path}/1_cos.fq_clean.sam")
    test_df2 = LISseq._read_sam_to_df(f"{expected_output_path}/2_cos.fq_clean.sam")
    assert (len(test_df1) == 5) & (len(test_df2) == 1)
    assert test_df1.columns.equals(pd.Index(["flag","chr","pos","Q"]))
    assert test_df2.columns.equals(pd.Index(["flag","chr","pos","Q"]))

def test_extract_IS(expected_output_path):
    test_df = LISseq._read_sam_to_df(f"{expected_output_path}/1_cos.fq_clean.sam")
    LISseq._extract_IS(test_df,20)
    assert False
