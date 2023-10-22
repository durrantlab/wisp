import os


def pytest_sessionstart(session):  # pytest_configure(config)
    r"""Called after the Session object has been created and
    before performing collection and entering the run test loop.
    """
    # Creates a tmp directory for writing files.
    test_tmp_path = "./tests/tmp"
    os.makedirs(test_tmp_path, exist_ok=True)
