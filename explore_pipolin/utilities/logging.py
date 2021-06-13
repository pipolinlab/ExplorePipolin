import functools
import logging
import os
from contextlib import contextmanager
from logging import Handler, LogRecord

from typing import MutableMapping

from prefect import context

from explore_pipolin.common import Genome
import explore_pipolin.settings as settings

_HANDLERS: MutableMapping[str, Handler] = {}


def _ensure_handler_for(genome: Genome):
    if genome.id not in _HANDLERS:
        os.makedirs(os.path.join(settings.get_instance().out_dir, genome.id), exist_ok=True)   # it's required for tests
        _HANDLERS[genome.id] = _create_logging_handler(
            genome=genome, out_dir=os.path.join(settings.get_instance().out_dir, genome.id))
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
        _ensure_handler_for(genome=genome)
        with _add_genome_id_to_logger(genome=genome):
            try:
                return func(genome, **kwargs)
            except Exception as e:
                logger: logging.Logger = context['logger']
                logger.exception(str(e))
                raise
    return wrapper
