import pytest

@pytest.fixture
def test_input_dir():
    return "test/test_input/"

@pytest.fixture
def correct_fqgz():
    return "test/test_input/0/0_cos.fq.gz"

@pytest.fixture
def incorrect_fqgz():
    return "test/test_input/3/3_cos.fq.gz"

@pytest.fixture
def default_test_args(tmp_path):
    return ["test/test_input", f"{tmp_path}"]

@pytest.fixture
def expected_output_path():
    return "test/expected_output/"

@pytest.fixture
def expected_IS_dict_raw_reads():
    return {"0":{"raw_reads":3, "filtered_reads":0},"1":{"raw_reads":5, "filtered_reads":3},"2":{"raw_reads":4, "filtered_reads":1}}

@pytest.fixture
def syspath_no_bowtie():
    return "/usr/bin"

@pytest.fixture
def hg38idx_default_dir():
    return "GRCh38_noalt_as"

@pytest.fixture
def hg38idx_nonexistant_dir():
    return "foo"

@pytest.fixture
def hg38idx_wrong_dir():
    return "test"





