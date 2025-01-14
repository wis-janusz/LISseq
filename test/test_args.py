import pytest
import pathlib
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
