import os
from Bio import SeqIO


class GermlineRecord:
    def __init__(self, id, seq):
        self.id = id
        id = id.split('|')
        self.accession_number = id[0][1:]
        self.gene_name = id[1]
        self.species = id[2]
        self.functionality = id[3]
        self.label = id[4]
        # self.start_IMGT_accession_number = int(id[5].split('.')[0])
        # self.end_IMGT_accession_number = int(id[5].split('.')[-1])
        self.length = int(id[6].split('_')[0])
        self.codon_start = id[7]
        self.reverse_complement = False if id[14] == '_' else True
        self.seq = seq.upper()

    def __str__(self):
        return ">%s\n%s\n" % (self.id, self.seq)

    def __repr__(self):
        return ">%s\n%s\n" % (self.id, self.seq)


def get_germline(path):
    germline = {}
    def read_from_file(filename, constant=False):
        raw_genes = list(SeqIO.parse(os.path.join(path, filename), 'fasta'))
        func_genes, pseudo_genes = {}, {}
        for gene in raw_genes:
            name = gene.name.split('|')
            if constant:
                name = "%s_%s" % (name[1], name[4])
            else:
                name = name[1]
            record = GermlineRecord(gene.name, gene.seq)
            if record.functionality == 'F':
                func_genes[name] = record
            else: pseudo_genes[name] = record
        return func_genes, pseudo_genes


    germline['IGHV'], germline['IGHVP'] = read_from_file('IGHV.fa')
    germline['IGHD'], germline['IGHDP'] = read_from_file('IGHD.fa')
    germline['IGHJ'], germline['IGHJP'] = read_from_file('IGHJ.fa')
    germline['IGHC'], germline['IGHCP'] = read_from_file('IGHC.fa', constant=True)

    germline['IGHV-allP'] = germline['IGHV'].copy()
    germline['IGHV-allP'].update(germline['IGHVP'])

    germline['IGHC-allP'] = germline['IGHC'].copy()
    germline['IGHC-allP'].update(germline['IGHCP'])

    return germline


# TODO finish reading all organisms, cells and chains
# class GermlineOrganismCellChain:
#     def __init__(self, path, cell, chain):
#         def read_from_file(filename, constant=False):
#             raw_genes = list(SeqIO.parse(os.path.join(path, filename), 'fasta'))
#             genes = {}
#             for gene in raw_genes:
#                 name = gene.name.split('|')
#                 if constant:
#                     name = "%s_%s" % (name[1], name[4])
#                 else:
#                     name = name[1]
#                 genes[name] = GermlineRecord(gene.name, gene.seq)
#             return genes
#
#         self.V = read_from_file(cell + chain + 'V.fa')
#         self.J = read_from_file(cell + chain + 'J.fa')
#         self.C = read_from_file(cell + chain + 'C.fa')
#         if chain == 'H':
#             self.D = read_from_file(cell + chain + 'D.fa')
