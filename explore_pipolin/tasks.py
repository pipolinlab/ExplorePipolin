"""
All the tasks should be defined in this file!
"""

import os
from prefect import task
from prefect import context
from prefect.engine import signals

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from explore_pipolin.utilities.easyfig_spec import add_colours, find_and_color_amr_and_virulence
from explore_pipolin.utilities import GQuery, Contig, Feature, Orientation, PipolinFragment
from explore_pipolin.utilities.io import read_blastxml, write_repeats, write_atts_denovo
from explore_pipolin.utilities import define_gquery_id
from explore_pipolin.utilities.io import read_seqio_records
from explore_pipolin.utilities.io import read_aragorn_batch
from explore_pipolin.utilities.atts_denovo_search import find_repeats
from explore_pipolin.utilities.external_tools import blast_genome_against_seq
from explore_pipolin.utilities.external_tools import run_prokka, run_aragorn
from explore_pipolin.utilities import create_fragment_record
from explore_pipolin.utilities.io import read_gff_records
from explore_pipolin.utilities import create_att_seqfeatures
from explore_pipolin.utilities.io import write_genbank_records
from explore_pipolin.utilities.io import write_gff_records
from explore_pipolin.utilities import add_new_gb_feature
from explore_pipolin.utilities import create_new_gb_record


@task
def create_gquery(genome) -> GQuery:
    gquery = GQuery(gquery_id=define_gquery_id(genome=genome))
    genome_dict = read_seqio_records(file=genome, file_format='fasta')
    for key, value in genome_dict.items():
        contig = Contig(contig_id=key, contig_length=len(value.seq))
        gquery.contigs.append(contig)

    return gquery


@task
def run_blast_against_polb(genome, root_dir, reference):
    blast_path = os.path.join(root_dir, 'polb_blast')
    blast_genome_against_seq(genome=genome, seq=reference, seq_type='protein', output_dir=blast_path)
    return blast_path


@task
def run_blast_against_att(genome, root_dir, reference):
    blast_path = os.path.join(root_dir, 'att_blast')
    blast_genome_against_seq(genome=genome, seq=reference, seq_type='nucleotide', output_dir=blast_path)
    return blast_path


@task
def add_features_from_blast(gquery: GQuery, blast_dir, feature_type):
    entries = read_blastxml(blast_xml=os.path.join(blast_dir, f'{gquery.gquery_id}.fmt5'))
    for entry in entries:
        for hit in entry:
            feature = gquery.feature_from_blasthit(hit=hit, contig_id=entry.id)
            gquery.get_features_by_type(feature_type).append(feature)


@task
def detect_trnas_with_aragorn(genome, root_dir):
    aragorn_results = os.path.join(root_dir, 'aragorn_results')
    run_aragorn(genome, aragorn_results)
    return aragorn_results


@task
def add_features_from_aragorn(gquery: GQuery, aragorn_dir):
    entries = read_aragorn_batch(aragorn_batch=os.path.join(aragorn_dir, f'{gquery.gquery_id}.batch'))
    for contig_id, hits in entries.items():
        for hit in hits:
            feature = Feature(start=hit[0], end=hit[1], frame=hit[2], contig=gquery.get_contig_by_id(contig_id))
            gquery.trnas.append(feature)

    for att in gquery.atts:
        target_trna = gquery.find_target_trna(att)
        if target_trna is not None:
            gquery.target_trnas.append(target_trna)


@task
def are_polbs_present(gquery: GQuery):
    logger = context.get('logger')

    if len(gquery.polbs) == 0:
        logger.warning('No piPolB! => No pipolins!')
        raise signals.FAIL()


@task()
def find_atts_denovo(genome, gquery: GQuery, root_dir):
    logger = context.get('logger')

    if not gquery.is_single_contig():
        logger.warning('This step is only for complete genomes. Pass...')
        raise signals.SKIP()

    atts_denovo_dir = os.path.join(root_dir, 'atts_denovo')
    os.makedirs(atts_denovo_dir, exist_ok=True)

    left_repeats, right_repeats, sequences = find_repeats(genome, gquery, atts_denovo_dir)
    write_repeats(gquery, left_repeats, atts_denovo_dir, right_repeats, sequences)

    atts_denovo = [(i, j) for i, j in zip(left_repeats, right_repeats) if gquery.is_att_denovo(i, j)]
    write_atts_denovo(atts_denovo, gquery, atts_denovo_dir)

    return atts_denovo_dir


@task
def add_features_atts_denovo(gquery: GQuery, atts_denovo_dir):
    with open(os.path.join(atts_denovo_dir, gquery.gquery_id + '.atts')) as inf:
        _ = inf.readline()
        for line in inf:
            att_pair = [int(i) for i in line.strip().split(sep='\t')]
            gquery.denovo_atts.append(Feature(start=att_pair[0], end=att_pair[1],
                                              frame=Orientation.FORWARD, contig=gquery.contigs[0]))
            gquery.denovo_atts.append(Feature(start=att_pair[2], end=att_pair[3],
                                              frame=Orientation.FORWARD, contig=gquery.contigs[0]))

    for att in gquery.denovo_atts:
        target_trna = gquery.find_target_trna(att)
        if target_trna is not None:
            gquery.target_trnas_denovo.append(target_trna)


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


@task
def analyse_pipolin_orientation(gquery):
    gquery.is_single_target_trna_per_contig()
    for contig in gquery.contigs:
        contig.contig_orientation = gquery.get_contig_orientation(contig)


@task
def scaffold_pipolins(gquery: GQuery):
    # Useful link to check feature's qualifiers: https://www.ebi.ac.uk/ena/WebFeat/
    # https://github.com/biopython/biopython/issues/1755
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


@task
def extract_pipolin_regions(genome, gquery: GQuery, root_dir):
    if gquery.pipolin_fragments is None:
        raise AssertionError('No pipolin fragments, but should be!!!')
    genome_dict = read_seqio_records(file=genome, file_format='fasta')

    pipolins_dir = os.path.join(root_dir, 'pipolin_sequences')
    os.makedirs(pipolins_dir, exist_ok=True)

    with open(os.path.join(pipolins_dir, gquery.gquery_id + '.fa'), 'w') as ouf:
        record = SeqRecord(seq='')
        sep_record = SeqRecord(seq='N' * 100)

        for fragment in gquery.pipolin_fragments[:-1]:
            fragment_record = create_fragment_record(fragment=fragment, genome_dict=genome_dict)
            print(f'@fragment length {len(fragment_record)} from {fragment.contig.contig_id}')
            record += fragment_record
            record += sep_record

        last_record = create_fragment_record(fragment=gquery.pipolin_fragments[-1], genome_dict=genome_dict)
        print(f'@fragment length {len(last_record)} from {gquery.pipolin_fragments[-1].contig.contig_id}')
        record += last_record

        print(f'@@@ total record length {len(record)}')

        record.id = gquery.gquery_id
        record.name = gquery.gquery_id
        record.description = gquery.gquery_id
        SeqIO.write(sequences=record, handle=ouf, format='fasta')

    return pipolins_dir


@task
def annotate_pipolins(gquery, pipolins_dir, proteins, root_dir):
    prokka_dir = os.path.join(root_dir, 'prokka')
    run_prokka(gquery.gquery_id, pipolins_dir, proteins, prokka_dir)
    return prokka_dir


@task
def include_atts_into_annotation(gquery, prokka_dir, root_dir):
    gb_records = read_seqio_records(file=os.path.join(prokka_dir, gquery.gquery_id + '.gbk'), file_format='genbank')
    gff_records = read_gff_records(file=os.path.join(prokka_dir, gquery.gquery_id + '.gff'))

    att_seqfeatures = create_att_seqfeatures(record_format='gb', gquery=gquery)
    for att in att_seqfeatures:
        add_new_gb_feature(new_feature=att, record=gb_records[gquery.gquery_id])
    att_seqfeatures = create_att_seqfeatures(record_format='gff', gquery=gquery)
    for att in att_seqfeatures:
        add_new_gb_feature(new_feature=att, record=gff_records[gquery.gquery_id])

    prokka_atts_dir = os.path.join(root_dir, 'prokka_atts')
    os.makedirs(prokka_atts_dir, exist_ok=True)

    write_genbank_records(gb_records=gb_records, out_dir=prokka_atts_dir, gquery=gquery)
    write_gff_records(gff_records=gff_records, out_dir=prokka_atts_dir, gquery=gquery)

    return prokka_atts_dir


@task
def easyfig_add_colours(gquery: GQuery, in_dir, abricate_dir):
    gb_records = read_seqio_records(file=os.path.join(in_dir, gquery.gquery_id + '.gbk'), file_format='genbank')
    add_colours(gb_records[gquery.gquery_id])
    if abricate_dir is not None:
        find_and_color_amr_and_virulence(gquery, gb_records, abricate_dir)
    write_genbank_records(gb_records=gb_records, out_dir=in_dir, gquery=gquery)


@task
def set_correct_positions(gquery: GQuery, prokka_atts_dir, root_dir):
    gb_records = read_seqio_records(file=os.path.join(prokka_atts_dir, gquery.gquery_id + '.gbk'),
                                    file_format='genbank')
    # gff_records = read_gff_records(file=os.path.join(prokka_atts_dir, gquery.gquery_id + '.gff'))

    new_gb_records = {gquery.gquery_id: create_new_gb_record(gquery=gquery, gb_record=gb_records[gquery.gquery_id])}

    prokka_atts_positions_dir = os.path.join(root_dir, 'prokka_atts_positions')
    os.makedirs(prokka_atts_positions_dir, exist_ok=True)

    write_genbank_records(gb_records=new_gb_records, out_dir=prokka_atts_positions_dir, gquery=gquery)
    # write_gff_records(gff_records=gff_records, out_dir=prokka_atts_positions_dir, gquery=gquery)

    return prokka_atts_positions_dir
