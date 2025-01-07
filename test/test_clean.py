import pytest
from Bio import SeqIO
import LISseq_clean

def test_read_write(default_io_files, tmp_path):
    LISseq_clean.main(default_io_files)
    assert len(list(tmp_path.iterdir())) == 1

def test_default_ltr(default_io_files, default_ltr):
    assert LISseq_clean.parse_args(default_io_files).ltr == default_ltr

def test_custom_ltr(cmd_with_custom_ltr,custom_ltr):
    assert LISseq_clean.parse_args(cmd_with_custom_ltr).ltr == custom_ltr

def test_out_file_contents(default_io_files, tmp_path, correct_output_path):
    LISseq_clean.main(default_io_files)
    expected_out = SeqIO.parse(f"{correct_output_path}", format="fastq")
    expected_ids = [seq.id for seq in expected_out]
    tested_out = SeqIO.parse(f"{tmp_path}/test_out.fastq", format="fastq")
    tested_ids = [seq.id for seq in tested_out]
    assert expected_ids == tested_ids