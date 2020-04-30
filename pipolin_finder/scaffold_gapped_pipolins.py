from prefect import task
from pipolin_finder.utilities import PipolinFragment, GQuery

# Useful link to check feature's qualifiers: https://www.ebi.ac.uk/ena/WebFeat/
# https://github.com/biopython/biopython/issues/1755


# def create_assembly_gap_record(record):
#     source_feature = SeqFeature(type='source', location=FeatureLocation(1, 100, strand=+1),
#                                 qualifiers={'mol_type': record.features[0].qualifiers['mol_type'],
#                                             'organism': record.features[0].qualifiers['organism'],
#                                             'strain': record.features[0].qualifiers['strain']})
#     assembly_gap_seq = Seq('N' * 100, alphabet=IUPACAmbiguousDNA())
#     assembly_gap_qualifiers = {'estimated_length': ['unknown'],
#                                'gap_type': ['within_scaffolds'],
#                                'linkage_evidence': ['pipolin_structure']}
#     assembly_gap_feature = SeqFeature(type='assembly_gap',
#                                       location=FeatureLocation(1, 100, strand=+1),
#                                       qualifiers=assembly_gap_qualifiers)
#     assembly_gap_record = SeqRecord(seq=assembly_gap_seq, id=record.id, name=record.name,
#                                     description=record.description, features=[source_feature, assembly_gap_feature],
#                                     annotations=record.annotations)
#
#     return assembly_gap_record


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
