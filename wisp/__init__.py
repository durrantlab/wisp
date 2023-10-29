from loguru import logger

logger.disable("wisp")


def enable_logging():
    logger.enable("wisp")
