import pytest

@pytest.fixture
def default_io_files(tmp_path):
    return ["test/test.fastq", f"{tmp_path}/test_out.fastq"]

@pytest.fixture
def cmd_with_custom_ltr(tmp_path):
    return ["test/test.fastq", f"{tmp_path}/test_out.fastq", "-ltr", "CCCCCC"]

@pytest.fixture
def custom_ltr():
    return "CCCCCC"

@pytest.fixture
def default_ltr():
    return "GGAGTGAATTAGCCCTTCCA"

@pytest.fixture
def correct_output_path():
    return "test/test_out.fastq"