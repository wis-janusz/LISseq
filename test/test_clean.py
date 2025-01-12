import pytest
import pathlib
from unittest.mock import patch, MagicMock 
import gzip
from Bio import SeqIO, SeqRecord, Seq
from LISseq import _parse_args, _find_fqgz, _parse_fqgz, _clean_read, _save_clean_reads, _cleanup_reads

class TestParseArgs:
    def test_parse_args_happy_path(self):
        args = _parse_args([
            "input_dir_path", 
            "output_dir_path", 
            "-ltr", "ACTGACTGACTG", 
            "-save_all_loci", 
            "-A", "10", 
            "-q", "30", 
            "-l", "25", 
            "-idx", "CustomIndex"
        ])
        assert args.input_dir == "input_dir_path"
        assert args.output_dir == "output_dir_path"
        assert args.ltr == "ACTGACTGACTG"
        assert args.save_all_loci is True
        assert args.A == 10
        assert args.q == 30
        assert args.l == 25
        assert args.idx == "CustomIndex"

    def test_parse_args_missing_required(self):
        with pytest.raises(SystemExit):
            _parse_args([])

    def test_parse_args_optional_defaults(self):
        args = _parse_args(["input_dir_path", "output_dir_path"])
        assert args.ltr == "GGAGTGAATTAGCCCTTCCA"
        assert args.save_all_loci is False
        assert args.A == 5
        assert args.q == 20
        assert args.l == 20
        assert args.idx == "GRCh38_noalt_as"

class TestFindFqgz:
    def test_find_fqgz_happy_path(self, tmp_path):
        # Setup: Create a directory structure with .fq.gz files
        (tmp_path / "subdir").mkdir()
        (tmp_path / "file1.fq.gz").touch()
        (tmp_path / "subdir" / "file2.fq.gz").touch()

        # Execute
        result = _find_fqgz(str(tmp_path))

        # Verify
        expected = [
            tmp_path / "file1.fq.gz",
            tmp_path / "subdir" / "file2.fq.gz"
        ]
        assert sorted(result) == sorted(expected)

    def test_find_fqgz_no_files(self, tmp_path):
        # Setup: Create an empty directory
        (tmp_path / "subdir").mkdir()

        # Execute
        result = _find_fqgz(str(tmp_path))

        # Verify
        assert result == []

    def test_find_fqgz_non_existent_directory(self):
        # Setup: Define a non-existent directory path
        non_existent_path = "/non/existent/directory"

        # Execute & Verify
        with pytest.raises(FileNotFoundError):
            _find_fqgz(non_existent_path)

    def test_find_fqgz_not_a_directory(self, tmp_path):
        # Setup: Define a non-existent directory path
        (tmp_path / "file1.fq.gz").touch()
        filepath = f"{tmp_path}/file1.fq.gz"

        # Execute & Verify
        with pytest.raises(NotADirectoryError):
            _find_fqgz(filepath)

class TestParseFqgz:
    def test_parse_fqgz_valid_file(self, tmp_path):
        # Create a valid .fq.gz file
        valid_fastq_content = "@SEQ_ID\nGATTA\n+\n!!!!!\n"
        valid_gz_path = tmp_path / "valid.fq.gz"
        with gzip.open(valid_gz_path, "wt") as f:
            f.write(valid_fastq_content)

        # Test the function
        records = list(_parse_fqgz(valid_gz_path))
        assert len(records) == 1
        assert isinstance(records[0], SeqRecord.SeqRecord)
        assert str(records[0].seq) == "GATTA"

    def test_parse_fqgz_invalid_file(self, tmp_path):
        # Create an invalid .fq.gz file
        invalid_gz_path = tmp_path / "invalid.fq.gz"
        with gzip.open(invalid_gz_path, "wb") as f:
            f.write(b"not a fastq content")

        # Test the function
        records = list(_parse_fqgz(invalid_gz_path))
        assert len(records) == 0

class TestCleanRead:
    def test_clean_read_removes_ltr_and_truncates_polyA(self):
        read = SeqRecord.SeqRecord(
            seq=Seq.Seq("LTRACGTACGTAAAAA"),
            id="test",
            letter_annotations={"phred_quality": [40] * 16}
        )
        cleaned_read = _clean_read(read, "LTR", 5, 30, 5)
        assert cleaned_read.seq == Seq.Seq("ACGTACGT")

    def test_clean_read_minimum_length_requirement(self):
        read = SeqRecord.SeqRecord(
            seq=Seq.Seq("LTRACGT"),
            id="test",
            letter_annotations={"phred_quality": [40] * 7}
        )
        cleaned_read = _clean_read(read, "LTR", 5, 30, 10)
        assert cleaned_read.seq == Seq.Seq("")

    def test_clean_read_quality_and_n_filtering(self):
        read_with_n = SeqRecord.SeqRecord(
            seq=Seq.Seq("LTRACGTN"),
            id="test",
            letter_annotations={"phred_quality": [40] * 8}
        )
        cleaned_read_with_n = _clean_read(read_with_n, "LTR", 5, 30, 5)
        assert cleaned_read_with_n.seq == Seq.Seq("")

        read_low_quality = SeqRecord.SeqRecord(
            seq=Seq.Seq("LTRACGT"),
            id="test",
            letter_annotations={"phred_quality": [20] * 7}
        )
        cleaned_read_low_quality = _clean_read(read_low_quality, "LTR", 5, 30, 5)
        assert cleaned_read_low_quality.seq == Seq.Seq("")


class TestCleanupReads:
    @patch('LISseq._find_fqgz')
    @patch('LISseq._parse_fqgz')
    @patch('LISseq._clean_read')
    @patch('LISseq._save_clean_reads')
    def test_cleanup_reads_identifies_all_files(self, mock_save_clean_reads, mock_clean_read, mock_parse_fqgz, mock_find_fqgz):
        # Setup
        mock_find_fqgz.return_value = [pathlib.Path("sample1_desc.fq.gz"), pathlib.Path("sample2_desc.fq.gz")]
        mock_parse_fqgz.return_value = [MagicMock(), MagicMock()]
        mock_clean_read.side_effect = lambda read, ltr, A, q, l: read

        args = MagicMock()
        args.input_dir = "input_dir"
        args.output_dir = "output_dir"
        args.ltr = "LTR"
        args.A = 10
        args.q = 20
        args.l = 30

        # Execute
        result = _cleanup_reads(args)

        # Verify
        assert len(result) == 2
        assert "sample1" in result
        assert "sample2" in result

    @patch('LISseq._find_fqgz')
    def test_cleanup_reads_no_files(self, mock_find_fqgz):
        # Setup
        mock_find_fqgz.return_value = []

        args = MagicMock()
        args.input_dir = "input_dir"
        args.output_dir = "output_dir"
        args.ltr = "LTR"
        args.A = 10
        args.q = 20
        args.l = 30

        # Execute
        result = _cleanup_reads(args)

        # Verify
        assert result == {}

    @patch('LISseq._find_fqgz')
    @patch('LISseq._parse_fqgz')
    @patch('LISseq._clean_read')
    def test_cleanup_reads_all_filtered_out(self, mock_clean_read, mock_parse_fqgz, mock_find_fqgz):
        # Setup
        mock_find_fqgz.return_value = [pathlib.Path("sample1_desc.fq.gz")]
        mock_parse_fqgz.return_value = iter([MagicMock(), MagicMock()])
        mock_clean_read.return_value = MagicMock(id="")

        args = MagicMock()
        args.input_dir = "input_dir"
        args.output_dir = "output_dir"
        args.ltr = "LTR"
        args.A = 10
        args.q = 20
        args.l = 30

        # Execute
        result = _cleanup_reads(args)

        # Verify
        assert result == {"sample1":{"raw_reads":2,"filtered_reads":0}}
