import pytest
import pathlib
import os
from unittest.mock import patch, MagicMock
import pandas as pd
from src.LISseq import (
    _find_fq,
    _check_bowtie2,
    _check_genome_index,
    _bowtie_map,
    _read_sam_to_df,
    _extract_IS,
    map_IS,
    format_data_frame,
)


class TestFindFq:
    @staticmethod
    def test_find_fq_happy_path(tmp_path):
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

    @staticmethod
    def test_find_fq_no_files(tmp_path):
        # Create a temporary directory structure without .fq files
        os.makedirs(f"{tmp_path}/temp", exist_ok=True)

        # Call the function and check the result
        result = _find_fq(tmp_path)
        assert result == []


class TestCheckBowtie2:
    @staticmethod
    @patch("shutil.which", return_value="/usr/bin/bowtie2")
    def test_check_bowtie2_available(mock_which):
        assert _check_bowtie2() == True

    @staticmethod
    @patch("shutil.which", return_value=None)
    def test_check_bowtie2_not_found(mock_which):
        with pytest.raises(
            FileNotFoundError,
            match="bowtie2 not found. Please install it and add it to PATH.",
        ):
            _check_bowtie2()

    @staticmethod
    @patch("shutil.which", return_value="/custom/path/bowtie2")
    def test_check_bowtie2_custom_path(mock_which):
        assert _check_bowtie2(path="/custom/path") == True


class TestCheckGenomeIndex:
    @staticmethod
    @patch("pathlib.Path.glob")
    def test_genome_index_found(mock_glob):
        # Mock the return value of glob to simulate index files present
        mock_glob.return_value = [
            pathlib.Path("genome.1.bt2"),
            pathlib.Path("genome.2.bt2"),
        ]
        idx_name = _check_genome_index("dummy_dir")
        assert idx_name == "genome"

    @staticmethod
    @patch("pathlib.Path.glob")
    def test_genome_index_not_found(mock_glob):
        # Mock the return value of glob to simulate no index files present
        mock_glob.return_value = []
        with pytest.raises(
            FileNotFoundError,
            match="No bowtie2-indexed genome found in dummy_dir, please download or build one.",
        ):
            _check_genome_index("dummy_dir")

    @staticmethod
    @patch("pathlib.Path.glob")
    def test_multiple_genome_indices(mock_glob):
        # Mock the return value of glob to simulate multiple index files
        mock_glob.return_value = [
            pathlib.Path("genome1.1.bt2"),
            pathlib.Path("genome2.1.bt2"),
        ]
        idx_name = _check_genome_index("dummy_dir")
        assert idx_name == "genome1"


class TestBowtieMap:
    @staticmethod
    def test_bowtie_map_success():
        clean_fq_file = pathlib.Path("/path/to/clean.fq")
        idx_dir = "/path/to/index"
        idx_name = "index_name"
        out_dir = "/path/to/output"
        expected_out_file = f"{out_dir}/temp/{clean_fq_file.stem}.sam"

        with patch("subprocess.run") as mock_run:
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

    @staticmethod
    def test_bowtie_map_command_failure():
        clean_fq_file = pathlib.Path("/path/to/clean.fq")
        idx_dir = "/path/to/index"
        idx_name = "index_name"
        out_dir = "/path/to/output"

        with patch("subprocess.run") as mock_run:
            mock_run.return_value.returncode = 1
            mock_run.return_value.stderr = "Error: bowtie2 command failed"
            with pytest.raises(Exception, match="Error: bowtie2 command failed"):
                _bowtie_map(clean_fq_file, idx_dir, idx_name, out_dir)

    @staticmethod
    def test_bowtie_map_special_characters_in_path():
        clean_fq_file = pathlib.Path("/path/to/clean file@#$.fq")
        idx_dir = "/path/to/index"
        idx_name = "index_name"
        out_dir = "/path/to/output"
        expected_out_file = f"{out_dir}/temp/{clean_fq_file.stem}.sam"

        with patch("subprocess.run") as mock_run:
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
    @staticmethod
    def test_read_sam_to_df_conversion(tmp_path):
        # Create a temporary SAM file with mock data
        sam_content = """@HD	VN:1.5	SO:unsorted	GO:query
@SQ	SN:chr1	LN:248956422
read1	0	chr1	158447409	42	30M	*	0	0	GTTATCCCTCAGAATTAATATTTGTCTTTC	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:30	YT:Z:UU
read2	16	chr1	158447409	42	30M	*	0	0	GTTATCCCTCAGAATTAATATTTGTCTTTC	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:30	YT:Z:UU
"""
        sam_file = tmp_path / "test.sam"
        with open(sam_file, mode="w") as in_file:
            in_file.write(sam_content)
        # Call the function
        df = _read_sam_to_df(str(sam_file))

        # Assert the DataFrame is as expected
        expected_data = {
            "flag": [0, 16],
            "chr": ["chr1", "chr1"],
            "pos": [158447409, 158447438],
            "Q": [42, 42],
            "seq": [
                "5'LTR-GTTATCCCTCAGAATTAATATTTGTCTTTC",
                "GTTATCCCTCAGAATTAATATTTGTCTTTC-5'LTR",
            ],
        }
        expected_df = pd.DataFrame(expected_data, index=["read1", "read2"])
        expected_df.index.names = ["read"]
        pd.testing.assert_frame_equal(df, expected_df)

    @staticmethod
    def test_read_sam_to_df_empty_file(tmp_path):
        # Create an empty SAM file
        sam_file = tmp_path / "empty.sam"
        sam_file.write_text("")

        # Call the function with empty file
        with pytest.raises(ValueError):
            _read_sam_to_df(str(sam_file))

    @staticmethod
    def test_read_sam_to_df_non_existent_file():
        # Call the function with a non-existent file path
        with pytest.raises(FileNotFoundError):
            _read_sam_to_df("non_existent.sam")


class TestExtractIS:
    @staticmethod
    def test_extract_IS_with_valid_data():
        data = {
            "chr": ["chr1", "chr1", "chr1", "chr2", "chr2"],
            "pos": [100, 100, 200, 300, 300],
            "flag": [0, 0, 16, 0, 16],
            "Q": [30, 30, 25, 20, 20],
            "seq": ["ATCG", "ATCG", "GCTA", "CGTA", "TACG"],
        }
        mappings_df = pd.DataFrame(data)
        Q = 20

        loci, total_mapped = _extract_IS(mappings_df, Q)

        assert total_mapped == 5
        assert len(loci) == 0  # No loci with depth >= 1000

    @staticmethod
    def test_extract_IS_with_high_quality_threshold():
        data = {
            "chr": ["chr1", "chr1", "chr1", "chr2", "chr2"],
            "pos": [100, 100, 200, 300, 300],
            "flag": [0, 0, 16, 0, 16],
            "Q": [30, 30, 25, 20, 20],
            "seq": ["ATCG", "ATCG", "GCTA", "CGTA", "TACG"],
        }
        mappings_df = pd.DataFrame(data)
        Q = 35

        loci, total_mapped = _extract_IS(mappings_df, Q)

        assert total_mapped == 0
        assert len(loci) == 0

    @staticmethod
    def test_extract_IS_with_loci_above_threshold():
        data = {
            "chr": ["chr1"] * 1000,
            "pos": [100] * 1000,
            "flag": [0] * 1000,
            "Q": [30] * 1000,
            "seq": ["ATCG"] * 1000,
        }
        mappings_df = pd.DataFrame(data)
        Q = 20

        loci, total_mapped = _extract_IS(mappings_df, Q)

        assert total_mapped == 1000
        assert len(loci) == 1
        assert loci.iloc[0]["depth"] == 1000
        assert loci.iloc[0]["mean_q"] == 30
        assert loci.iloc[0]["seq"] == "ATCG"

    @staticmethod
    def test_extract_IS_with_no_valid_loci():
        data = {
            "chr": ["chr1", "chr1", "chr1", "chr2", "chr2"],
            "pos": [100, 100, 200, 300, 300],
            "flag": [0, 0, 16, 0, 16],
            "Q": [10, 10, 15, 5, 5],
            "seq": ["ATCG", "ATCG", "GCTA", "CGTA", "TACG"],
        }
        mappings_df = pd.DataFrame(data)
        Q = 20

        loci, total_mapped = _extract_IS(mappings_df, Q)

        assert total_mapped == 0
        assert len(loci) == 0


class TestFormatDataFrame:
    @staticmethod
    @patch("src.LISseq._get_gene", return_value="GENE1")
    def test_format_data_frame_happy_path(mock_get_gene):
        IS_dict = {
            "sample1": {
                "chr": "chr1",
                "pos": 12345,
                "strand": "+",
                "raw_reads": 1000,
                "filtered_reads": 800,
                "total_mapped": 750,
                "mapped_to_locus": 700,
                "mapping_quality": 60.5,
                "seq": "ATCG",
            },
            "sample2": {
                "chr": "chr2",
                "pos": 67890,
                "strand": "-",
                "raw_reads": 1500,
                "filtered_reads": 1200,
                "total_mapped": 1100,
                "mapped_to_locus": 1050,
                "mapping_quality": 55.3,
                "seq": "CGTA",
            },
        }

        expected_data = {
            "chr": ["chr1", "chr2"],
            "pos": [12345, 67890],
            "strand": ["+", "-"],
            "gene": ["GENE1", "GENE1"],
            "raw_reads": [1000, 1500],
            "filtered_reads": [800, 1200],
            "total_mapped": [750, 1100],
            "mapped_to_locus": [700, 1050],
            "mapping_quality": [60.5, 55.3],
            "seq": ["ATCG", "CGTA"],
        }
        expected_df = pd.DataFrame(expected_data, index=["sample1", "sample2"])

        result_df = format_data_frame(IS_dict)
        pd.testing.assert_frame_equal(result_df, expected_df)

    @staticmethod
    @patch("src.LISseq._get_gene", return_value="intergenic")
    def test_format_data_frame_intergenic(mock_get_gene):
        IS_dict = {
            "sample1": {
                "chr": "chr1",
                "pos": 12345,
                "strand": "+",
                "raw_reads": 1000,
                "filtered_reads": 800,
                "total_mapped": 750,
                "mapped_to_locus": 700,
                "mapping_quality": 60.5,
                "seq": "ATCG",
            }
        }

        expected_data = {
            "chr": ["chr1"],
            "pos": [12345],
            "strand": ["+"],
            "gene": ["intergenic"],
            "raw_reads": [1000],
            "filtered_reads": [800],
            "total_mapped": [750],
            "mapped_to_locus": [700],
            "mapping_quality": [60.5],
            "seq": ["ATCG"],
        }
        expected_df = pd.DataFrame(expected_data, index=["sample1"])

        result_df = format_data_frame(IS_dict)
        pd.testing.assert_frame_equal(result_df, expected_df)

    @staticmethod
    @patch("src.LISseq._get_gene", return_value=None)
    def test_format_data_frame_no_gene(mock_get_gene):
        IS_dict = {
            "sample1": {
                "chr": "chr1",
                "pos": 12345,
                "strand": "+",
                "raw_reads": 1000,
                "filtered_reads": 800,
                "total_mapped": 750,
                "mapped_to_locus": 700,
                "mapping_quality": 60.5,
                "seq": "ATCG",
            }
        }

        expected_data = {
            "chr": ["chr1"],
            "pos": [12345],
            "strand": ["+"],
            "gene": [None],
            "raw_reads": [1000],
            "filtered_reads": [800],
            "total_mapped": [750],
            "mapped_to_locus": [700],
            "mapping_quality": [60.5],
            "seq": ["ATCG"],
        }
        expected_df = pd.DataFrame(expected_data, index=["sample1"])

        result_df = format_data_frame(IS_dict)
        pd.testing.assert_frame_equal(result_df, expected_df)


class TestMapIS:
    @staticmethod
    @patch("src.LISseq._find_fq", return_value=[pathlib.Path("/path/to/clean1.fq"), pathlib.Path("/path/to/clean2.fq")])
    @patch("src.LISseq._check_bowtie2", return_value=True)
    @patch("src.LISseq._check_genome_index", return_value="index_name")
    @patch("src.LISseq._bowtie_map", side_effect=["/path/to/output/temp/clean1.sam", "/path/to/output/temp/clean2.sam"])
    @patch("src.LISseq._read_sam_to_df")
    @patch("src.LISseq._extract_IS")
    def test_map_IS_happy_path(mock_extract_IS, mock_read_sam_to_df, mock_bowtie_map, mock_check_genome_index, mock_check_bowtie2, mock_find_fq):
        args = MagicMock()
        args.output_dir = "/path/to/output"
        args.idx = "/path/to/index"
        args.q = 20
        read_no_dict = {
            "clean1": {"raw_reads": 1000, "filtered_reads": 800},
            "clean2": {"raw_reads": 1500, "filtered_reads": 1200},
        }

        mock_read_sam_to_df.side_effect = [
            pd.DataFrame({
                "flag": [0, 16],
                "chr": ["chr1", "chr1"],
                "pos": [100, 200],
                "Q": [42, 42],
                "seq": ["ATCG", "GCTA"],
            }),
            pd.DataFrame({
                "flag": [0, 16],
                "chr": ["chr2", "chr2"],
                "pos": [300, 400],
                "Q": [42, 42],
                "seq": ["CGTA", "TACG"],
            }),
        ]

        mock_extract_IS.side_effect = [
            (pd.DataFrame({
                "chr": ["chr1"],
                "pos": [100],
                "flag": [0],
                "depth": [1000],
                "mean_q": [42],
                "seq": ["ATCG"],
            }), 2),
            (pd.DataFrame({
                "chr": ["chr2"],
                "pos": [300],
                "flag": [0],
                "depth": [1000],
                "mean_q": [42],
                "seq": ["CGTA"],
            }), 2),
        ]

        result = map_IS(args, read_no_dict)

        expected_result = {
            "clean1": {
                0: {
                    "chr": "chr1",
                    "pos": 100,
                    "strand": "-",
                    "raw_reads": 1000,
                    "filtered_reads": 800,
                    "total_mapped": 2,
                    "mapped_to_locus": 1000,
                    "mapping_quality": 42.0,
                    "seq": "ATCG",
                }
            },
            "clean2": {
                0: {
                    "chr": "chr2",
                    "pos": 300,
                    "strand": "-",
                    "raw_reads": 1500,
                    "filtered_reads": 1200,
                    "total_mapped": 2,
                    "mapped_to_locus": 1000,
                    "mapping_quality": 42.0,
                    "seq": "CGTA",
                }
            },
        }

        assert result == expected_result

    @staticmethod
    @patch("src.LISseq._find_fq", return_value=[])
    @patch("src.LISseq._check_bowtie2", return_value=True)
    @patch("src.LISseq._check_genome_index", return_value="index_name")
    def test_map_IS_no_files(mock_check_genome_index, mock_check_bowtie2, mock_find_fq):
        args = MagicMock()
        args.output_dir = "/path/to/output"
        args.idx = "/path/to/index"
        args.q = 20
        read_no_dict = {}

        result = map_IS(args, read_no_dict)
        assert result == {}
