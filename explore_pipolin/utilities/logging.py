import functools
import logging
import os
from contextlib import contextmanager
from logging import Handler, LogRecord

from typing import MutableMapping

from prefect import context

from explore_pipolin.common import Genome

_HANDLERS: MutableMapping[str, Handler] = {}
_LOG_DIR: str = 'logs'


def set_logging_dir(log_dir: str):
    global _LOG_DIR
    _LOG_DIR = log_dir


def _ensure_handler_for(genome: Genome):
    if genome.id not in _HANDLERS:
        _HANDLERS[genome.id] = _create_logging_handler(
            genome=genome, out_dir=os.path.join(_LOG_DIR, genome.id))
        logging.getLogger().addHandler(_HANDLERS[genome.id])


def _create_logging_handler(genome: Genome, out_dir: str):
    handler = logging.FileHandler(os.path.join(out_dir, f'{genome.id}.log'), mode='w')

    def my_filter(record: LogRecord):
        return 1 if (hasattr(record, 'genome_id') and record.genome_id == genome.id) else 0

    handler.addFilter(my_filter)
    handler.setFormatter(logging.Formatter('{asctime} {levelname}: {name} - {genome_id} - {message}', style='{'))
    return handler


@contextmanager
def _add_genome_id_to_logger(genome: Genome):
    logger: logging.Logger = context['logger']

    def log_record_filter(record: LogRecord):
        record.genome_id = genome.id
        return 1

    logger.addFilter(log_record_filter)
    logger.info('starting...')
    yield
    logger.info('done')
    logger.removeFilter(log_record_filter)


def genome_specific_logging(func):
    @functools.wraps(func)
    def wrapper(genome, **kwargs):
        os.makedirs(os.path.join(_LOG_DIR, genome.id), exist_ok=True)
        _ensure_handler_for(genome=genome)
        with _add_genome_id_to_logger(genome=genome):
            return func(genome, **kwargs)
    return wrapper
