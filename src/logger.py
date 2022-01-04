
import logging
from datetime import datetime


def get_logger(name: str) -> logging.Logger:
    """
    returns basic logger (console and file output)

    :param name:
        name of the logger (should relate to the calling module)
    """

    # setup logger -> will write log file to the /../log directory
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    # create file handler (has current timestamp in file name)
    now = datetime.now().strftime('%Y%m%d-%H%M%S')
    fh = logging.FileHandler(f'./../log/{name}_{now}.log')
    fh.setLevel(logging.INFO)
    # create console handler
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    # add the handlers to the logger
    logger.addHandler(fh)
    logger.addHandler(ch)

    return logger
