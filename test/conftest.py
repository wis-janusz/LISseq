import pytest

@pytest.fixture
def test_input_dir():
    return "test/test_input/"

@pytest.fixture
def correct_fqgz():
    return "test/test_input/0/0.fq.gz"

@pytest.fixture
def incorrect_fqgz():
    return "test/test_input/3/3.fq.gz"

@pytest.fixture
def default_test_args(tmp_path):
    return ["test/test_input", f"{tmp_path}"]

@pytest.fixture
def expected_output_path():
    return "test/expected_output/"

@pytest.fixture
def custom_ltr():
    return "CCCCCC"

@pytest.fixture
def default_ltr():
    return "GGAGTGAATTAGCCCTTCCA"

