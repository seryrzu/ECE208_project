import argparse
from Bio import SeqIO

from germline.germline import get_germline

def get_params():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--genome", required=True)
    parser.add_argument("-c", "--coords", required=True)
    parser.add_argument("-G", "--germline", default='../data/germline_new/Homo_sapiens/ig')
    parser.add_argument("-o", "--outfile", required=True)
    parser.add_argument("--rss-len", dest="rss_len", default=7+23+9)
    params = parser.parse_args()
    return params


def read_exact_positions(filename):
    res = {}
    with open(filename, 'r') as f:
        lines = f.readlines()
        for line in lines[1:]:
            k, v1, v2 = line.strip().split('\t')
            v1 = [x.strip('[] ') for x in v1.split(',')]
            v2 = [x.strip('[] ') for x in v2.split(',')]
            if v1[0] != '':
                res[k] = [int(x) for x in v1][0]
            elif v2[0] != '':
                res[k] = [int(x) for x in v2][0]
    return res


def extract_genes(genome, coords, germline, rss_len):
    genes = {}
    for gene, coord in coords.items():
        genes[gene] = (genome[coord-rss_len-1:coord-1], str(germline[gene].seq))
    return genes


def unflatten_dict(germline):
    res = {}
    for k in germline:
        for k2 in germline[k]:
            res[k2] = germline[k][k2]
    return res


def export_genes(genes, outfile):
    with open(outfile, 'w') as f:
        for gene_id in genes:
            gene = genes[gene_id]
            gene = gene[0] + gene[1]
            print('>'+gene_id+'\n'+gene, file=f)



def main():
    params = get_params()
    genome = str(list(SeqIO.parse(params.genome, 'fasta'))[0].seq)
    coords = read_exact_positions(params.coords)
    germline = unflatten_dict(get_germline(params.germline))
    genes = extract_genes(genome, coords, germline, params.rss_len)
    export_genes(genes, params.outfile)



if __name__ == "__main__":
    main()
