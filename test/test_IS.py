import pytest
import pathlib
import tempfile
import os
from LISseq import _find_fq

class TestFindFq:
    def test_find_fq_happy_path(self, tmp_path):
        # Create a temporary directory structure with .fq files
        os.makedirs(f"{tmp_path}/temp", exist_ok=True)
        fq_file_1 = pathlib.Path(f"{tmp_path}/temp/sample1.fq")
        fq_file_2 = pathlib.Path(f"{tmp_path}/temp/sample2.fq")
        fq_file_1.touch()
        fq_file_2.touch()

        # Call the function and check the result
        result = _find_fq(tmp_path)
        expected = [fq_file_1, fq_file_2]
        assert sorted(result) == sorted(expected)

    def test_find_fq_no_files(self,tmp_path):
        # Create a temporary directory structure without .fq files
        os.makedirs(f"{tmp_path}/temp", exist_ok=True)

        # Call the function and check the result
        result = _find_fq(tmp_path)
        assert result == []