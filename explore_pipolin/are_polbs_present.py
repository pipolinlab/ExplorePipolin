from prefect import task
from prefect import context
from prefect.engine import signals
from explore_pipolin.utilities import GQuery


@task
def are_polbs_present(gquery: GQuery):
    logger = context.get('logger')

    if len(gquery.polbs) == 0:
        logger.warning('No piPolB! => No pipolins!')
        raise signals.FAIL()
