import unittest
import tempfile

from explore_pipolin.download_metadata_ncbi import read_gb_ids
from explore_pipolin.download_metadata_ncbi import get_assembly_info
from explore_pipolin.download_metadata_ncbi import get_assembly_acc_and_species_name
from explore_pipolin.download_metadata_ncbi import get_and_filter_assembly_info
from explore_pipolin.download_metadata_ncbi import get_strain_and_acc_ids

from explore_pipolin.download_genomes_ncbi import download_genome_seqs
from explore_pipolin.download_genomes_ncbi import read_metadata_file
from explore_pipolin.download_genomes_ncbi import get_strain_names_and_gb_ids


class DownloadFromNCBITestCase(unittest.TestCase):
    def setUp(self) -> None:
        self.assembly = 'GCF_002860085.1'
        self.species = 'Escherichia coli'
        self.strain = 'CRE1540'
        self.gb_id = 'CP019051.1'

    def test_entrez(self):
        assembly_info = get_assembly_info(gb_id=self.gb_id)
        assembly_acc, species_name = get_assembly_acc_and_species_name(assembly_info)
        seqs_info = get_and_filter_assembly_info(assembly_info)
        strain, acc_ids = get_strain_and_acc_ids(seqs_info)

        self.assertEqual(assembly_acc, self.assembly)
        self.assertEqual(species_name, self.species)
        self.assertEqual(strain, self.strain)
        self.assertEqual(acc_ids, [self.gb_id])

        self.assertTrue(download_genome_seqs(self.gb_id))

    def test_read_gb_ids_from_file(self):
        with tempfile.NamedTemporaryFile('w') as inf:
            print('\n'.join(['LKJH12345', 'POIU_09876']), file=inf)
            inf.flush()
            ids = read_gb_ids(input_file=inf.name)
        self.assertEqual(ids, ['LKJH12345', 'POIU_09876'])

    def test_read_metadata_file(self):
        in_line = '\t'.join([self.assembly, self.species, self.strain, self.gb_id])
        with tempfile.NamedTemporaryFile('w') as inf:
            print(f'header\n{in_line}', file=inf)
            inf.flush()
            out_line = read_metadata_file(metadata_file=inf.name)

        strains, gb_ids = get_strain_names_and_gb_ids(out_line)
        self.assertEqual(strains, [self.strain])
        self.assertEqual(gb_ids, [self.gb_id])
