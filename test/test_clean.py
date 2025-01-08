import filecmp
from Bio import SeqIO, SeqRecord
import LISseq_clean

def test_fing_fqgz(test_input_dir):
    names = set(filename.stem for filename in LISseq_clean._find_fqgz(test_input_dir))
    assert names == {"0.fq","1.fq","2.fq","3.fq"}

def test_correct_fqgz(correct_fqgz):
    test_seq = LISseq_clean._parse_fqgz(correct_fqgz)
    for read in test_seq:
        assert type(read) == SeqRecord.SeqRecord

def test_incorrect_fqgz(incorrect_fqgz):
    test_seq = LISseq_clean._parse_fqgz(incorrect_fqgz)
    for read in test_seq:
        assert read == None

def test_main_function(default_test_args, expected_output_path, tmp_path):
    LISseq_clean.main(default_test_args)
    for i in range(3):
        with open(f"{tmp_path}/{i}.fq_clean.fq", mode="r") as test_file, open(f"{expected_output_path}/{i}.fq_clean.fq", mode="r") as exp_file:
            test_lines = test_file.readlines()
            exp_lines = exp_file.readlines()
            assert test_lines == exp_lines