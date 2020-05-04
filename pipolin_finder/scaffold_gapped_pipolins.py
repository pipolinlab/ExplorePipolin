from prefect import task
from pipolin_finder.utilities import PipolinFragment, GQuery

# Useful link to check feature's qualifiers: https://www.ebi.ac.uk/ena/WebFeat/
# https://github.com/biopython/biopython/issues/1755


@task
def scaffold_pipolins(gquery: GQuery):
    if gquery.is_single_contig() or gquery.is_on_the_same_contig():
        print('>>>Scaffolding is not required!')
        if len(gquery.dict_by_contig_normalized(gquery.atts)) != 0:
            start, end = gquery.get_pipolin_bounds()
            pipolin = PipolinFragment(contig=gquery.get_contig_by_id(gquery.polbs[0].contig.contig_id),
                                      start=start, end=end)

            pipolin.atts.extend(gquery.atts)
            gquery.pipolin_fragments.append(pipolin)
        else:
            left_window, right_window = gquery.get_left_right_windows()
            pipolin = PipolinFragment(contig=gquery.get_contig_by_id(gquery.polbs[0].contig.contig_id),
                                      start=left_window[0], end=right_window[1])
            gquery.pipolin_fragments.append(pipolin)
    else:
        print('>>>Scaffolding is required!')
        gquery.try_creating_single_record()
