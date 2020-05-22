import os
from Bio.SeqRecord import SeqRecord
from prefect import task
from Bio import SeqIO
from explore_pipolin.utilities import GQuery, Orientation
from explore_pipolin.utilities import read_seqio_records


def create_fragment_record(fragment, genome_dict):
    fragment_record = genome_dict[fragment.contig.contig_id][fragment.start:fragment.end]
    if fragment.contig.contig_orientation == Orientation.REVERSE:
        fragment_record = fragment_record.reverse_complement()
    return fragment_record


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
