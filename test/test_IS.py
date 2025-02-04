import pytest
import pathlib
import os
from unittest import mock
import pandas as pd
from src.LISseq import (
    _find_fq,
    _check_bowtie2,
    _check_genome_index,
    _bowtie_map,
    _read_sam_to_df,
)


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

    def test_find_fq_no_files(self, tmp_path):
        # Create a temporary directory structure without .fq files
        os.makedirs(f"{tmp_path}/temp", exist_ok=True)

        # Call the function and check the result
        result = _find_fq(tmp_path)
        assert result == []


class TestCheckBowtie2:
    @mock.patch("shutil.which", return_value="/usr/bin/bowtie2")
    def test_check_bowtie2_available(self, mock_which):
        assert _check_bowtie2() == True

    @mock.patch("shutil.which", return_value=None)
    def test_check_bowtie2_not_found(self, mock_which):
        with pytest.raises(
            FileNotFoundError,
            match="bowtie2 not found. Please install it and add it to PATH.",
        ):
            _check_bowtie2()

    @mock.patch("shutil.which", return_value="/custom/path/bowtie2")
    def test_check_bowtie2_custom_path(self, mock_which):
        assert _check_bowtie2(path="/custom/path") == True


class TestCheckGenomeIndex:

    @mock.patch("pathlib.Path.glob")
    def test_genome_index_found(self, mock_glob):
        # Mock the return value of glob to simulate index files present
        mock_glob.return_value = [
            pathlib.Path("genome.1.bt2"),
            pathlib.Path("genome.2.bt2"),
        ]
        idx_name = _check_genome_index("dummy_dir")
        assert idx_name == "genome"

    @mock.patch("pathlib.Path.glob")
    def test_genome_index_not_found(self, mock_glob):
        # Mock the return value of glob to simulate no index files present
        mock_glob.return_value = []
        with pytest.raises(
            FileNotFoundError,
            match="No bowtie2-indexed genome found in dummy_dir, please download or build one.",
        ):
            _check_genome_index("dummy_dir")

    @mock.patch("pathlib.Path.glob")
    def test_multiple_genome_indices(self, mock_glob):
        # Mock the return value of glob to simulate multiple index files
        mock_glob.return_value = [
            pathlib.Path("genome1.1.bt2"),
            pathlib.Path("genome2.1.bt2"),
        ]
        idx_name = _check_genome_index("dummy_dir")
        assert idx_name == "genome1"


class TestBowtieMap:
    def test_bowtie_map_success(self):
        clean_fq_file = pathlib.Path("/path/to/clean.fq")
        idx_dir = "/path/to/index"
        idx_name = "index_name"
        out_dir = "/path/to/output"
        expected_out_file = f"{out_dir}/temp/{clean_fq_file.stem}.sam"

        with mock.patch("subprocess.run") as mock_run:
            mock_run.return_value.returncode = 0
            result = _bowtie_map(clean_fq_file, idx_dir, idx_name, out_dir)
            assert result == expected_out_file
            mock_run.assert_called_once_with(
                [
                    "bowtie2",
                    "-x",
                    f"{idx_dir}/{idx_name}",
                    "-U",
                    "-q",
                    clean_fq_file,
                    "-S",
                    expected_out_file,
                ],
                capture_output=True,
                text=True,
            )

    def test_bowtie_map_command_failure(self):
        clean_fq_file = pathlib.Path("/path/to/clean.fq")
        idx_dir = "/path/to/index"
        idx_name = "index_name"
        out_dir = "/path/to/output"

        with mock.patch("subprocess.run") as mock_run:
            mock_run.return_value.returncode = 1
            mock_run.return_value.stderr = "Error: bowtie2 command failed"
            with pytest.raises(Exception, match="Error: bowtie2 command failed"):
                _bowtie_map(clean_fq_file, idx_dir, idx_name, out_dir)

    def test_bowtie_map_special_characters_in_path(self):
        clean_fq_file = pathlib.Path("/path/to/clean file@#$.fq")
        idx_dir = "/path/to/index"
        idx_name = "index_name"
        out_dir = "/path/to/output"
        expected_out_file = f"{out_dir}/temp/{clean_fq_file.stem}.sam"

        with mock.patch("subprocess.run") as mock_run:
            mock_run.return_value.returncode = 0
            result = _bowtie_map(clean_fq_file, idx_dir, idx_name, out_dir)
            assert result == expected_out_file
            mock_run.assert_called_once_with(
                [
                    "bowtie2",
                    "-x",
                    f"{idx_dir}/{idx_name}",
                    "-U",
                    "-q",
                    clean_fq_file,
                    "-S",
                    expected_out_file,
                ],
                capture_output=True,
                text=True,
            )


class TestReadSamToDf:
    def test_read_sam_to_df_conversion(self, tmp_path):
        # Create a temporary SAM file with mock data
        sam_content = """@HD	VN:1.5	SO:unsorted	GO:query
@SQ	SN:chr1	LN:248956422
read1	0	chr1	158447409	42	30M	*	0	0	GTTATCCCTCAGAATTAATATTTGTCTTTC	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:30	YT:Z:UU
read2	0	chr1	158447409	42	30M	*	0	0	GTTATCCCTCAGAATTAATATTTGTCTTTC	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:30	YT:Z:UU
"""
        sam_file = tmp_path / "test.sam"
        with open(sam_file, mode="w") as in_file:
            in_file.write(sam_content)
        # Call the function
        df = _read_sam_to_df(str(sam_file))

        # Assert the DataFrame is as expected
        expected_data = {
            "flag": [0, 0],
            "chr": ["chr1", "chr1"],
            "pos": [158447409, 158447409],
            "Q": [42, 42],
            "seq": [
                "GTTATCCCTCAGAATTAATATTTGTCTTTC",
                "GTTATCCCTCAGAATTAATATTTGTCTTTC",
            ],
        }
        expected_df = pd.DataFrame(expected_data, index=["read1", "read2"])
        expected_df.index.names = ["read"]
        pd.testing.assert_frame_equal(df, expected_df)


    def test_read_sam_to_df_empty_file(self, tmp_path):
        # Create an empty SAM file
        sam_file = tmp_path / "empty.sam"
        sam_file.write_text("")

        # Call the function with empty file
        with pytest.raises(ValueError):
            _read_sam_to_df(str(sam_file))


    def test_read_sam_to_df_non_existent_file(self):
        # Call the function with a non-existent file path
        with pytest.raises(FileNotFoundError):
            _read_sam_to_df("non_existent.sam")
