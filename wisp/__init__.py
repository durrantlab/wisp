__version__ = "0.0.0"

from loguru import logger

logger.disable("wisp")


def enable_logging(file_path=None):
    r"""Enable logging for WISP.

    Parameters
    ----------
    file_path : :obj:`str`, optional
        Write logs to file located at this path.
    """
    logger.enable("wisp")
    if isinstance(file_path, str):
        logger.add(file_path)
