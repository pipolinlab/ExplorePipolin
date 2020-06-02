from prefect import task
from prefect import context

from explore_pipolin.utilities import GQuery, Feature


@task(skip_on_upstream_skip=False)
def are_atts_present(gquery: GQuery):
    logger = context.get('logger')

    if len(gquery.atts) == 0 and len(gquery.denovo_atts) == 0:
        logger.warning('\n\n>>>There is piPolB, but no atts were found! Not able to define pipolin bounds!\n')
        # TODO: probably, it makes sense to output piPolB(s) alone
        # raise signals.SKIP() # let's try cutting from both sides and proceed with annotation

    elif len(gquery.atts) == 0:
        logger.warning(f'\n\n>>>No "usual" atts were found, but some atts were found by denovo search!'
                       f'For more details, check the {gquery.gquery_id}.atts file in the atts_denovo directory!\n')
        # TODO: check that it's only one repeat! Although, this shouldn't be a problem.
        atts_frames = [att.frame for att in gquery.denovo_atts]
        if len(set(atts_frames)) != 1:
            raise AssertionError('ATTs are expected to be in the same frame, as they are direct repeats!')
        if set(atts_frames).pop() == gquery.target_trnas_denovo[0].frame:
            reverse_denovo_atts = []
            for att in gquery.denovo_atts:
                reverse_denovo_atts.append(Feature(start=att.start, end=att.end, frame=-att.frame, contig=att.contig))
            gquery.atts.extend(reverse_denovo_atts)
        else:
            gquery.atts.extend(gquery.denovo_atts)
        gquery.target_trnas.extend(gquery.target_trnas_denovo)

    elif len(gquery.denovo_atts) != 0:
        logger.warning(f'\n\n>>>Some atts were found by denovo search, but we are not going to use them!'
                       f'For more details, check the {gquery.gquery_id}.atts file in the atts_denovo directory!\n')
