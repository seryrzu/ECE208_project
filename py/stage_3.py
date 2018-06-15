import argparse
import dendropy
import numpy as np
from Bio import SeqIO

from genes_extractor import read_exact_positions


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--tree", dest="tree_fp",
                        help="path to the branch length", required=True)
    parser.add_argument("-g", "--genome", required=True)
    parser.add_argument("-G", "--germline", default='../data/germline_new/Homo_sapiens/ig')
    parser.add_argument("-c", "--coords", required=True)
    parser.add_argument("--rss-len", dest="rss_len", default=7+23+9)
    parser.add_argument("-p", "--parameters", dest="param_fp",
                        help="path to the GTR parameter file", required=True)
    parser.add_argument("-o", "--output", dest="output",
                        help="output tree", required=True)
    params = parser.parse_args()
    return params


def extract_rss(genome, coords, rss_len):
    rss = {}

    for gene, coord in coords.items():
        rss[gene] = genome[coord-rss_len-1:coord-1]
    return rss


def read_param(param_fp):
    param_dct = dict()
    param_lst = open(param_fp,'r').readlines()
    param_dct["A"] = float(param_lst[0].strip().split()[0])
    param_dct["C"] = float(param_lst[0].strip().split()[1])
    param_dct["G"] = float(param_lst[0].strip().split()[2])
    param_dct["T"] = float(param_lst[0].strip().split()[3])
    param_dct["AC"] = float(param_lst[1].strip().split()[0])
    param_dct["AG"] = float(param_lst[1].strip().split()[1])
    param_dct["AT"] = float(param_lst[1].strip().split()[2])
    param_dct["CG"] = float(param_lst[1].strip().split()[3])
    param_dct["CT"] = float(param_lst[1].strip().split()[4])
    param_dct["GT"] = float(param_lst[1].strip().split()[5])
    return param_dct


def get_internal_nodes(tree, seq_dct, param_dct):
    '''
    :param tree: tree
    :param seq_dct: sequences dictionary, key: species name, values: sequences
    :param_dct: parameter file path
    :return: log likelihood of observing this tree
    '''
    pi = np.array([param_dct["A"], param_dct["C"], param_dct["G"], param_dct["T"]])
    nucl_index = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    R = []
    bases = list("ACGT")
    for i, n1 in enumerate(bases):
        R.append([])
        last_R = R[-1]
        for n2 in bases:
            if n2 == n1:
                last_R.append(0)
                continue
            pair = n1+n2 if n1 < n2 else n2+n1
            last_R.append(param_dct[pair] * param_dct[n2])
        last_R[i] = -sum(last_R)
    R = np.array(R)
    R /= -pi.dot(np.diag(R))

    def matrix_pow(A, threshold=1e-6):
        e_val, e_vec = np.linalg.eig(A)
        res = np.diag(np.ones_like(e_val))
        cur_e_val = np.ones_like(e_val)
        cnt = 1
        denominator = 1
        while np.abs(np.sum(cur_e_val) / denominator) > threshold:
            cur_e_val *= e_val
            res += np.diag(cur_e_val) / denominator
            cnt += 1
            denominator *= cnt
        return np.linalg.multi_dot((e_vec, res, np.linalg.inv(e_vec)))

    length = len(list(seq_dct.values())[0].strip())
    lklh_dct = {}
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            lklh_dct[node] = np.zeros((4, length))
            seq = seq_dct[node.taxon.label].strip()
            for i in range(length):
                lklh_dct[node][nucl_index[seq[i]]][i] = 1
            continue

        lklh_dct[node] = np.ones((4, length))
        for child_edge in node.child_edge_iter():
            Pt = matrix_pow(child_edge.length * R)
            lklh_dct[node] *= Pt.dot(lklh_dct[child_edge.head_node])

        argmaxes = np.argmax(lklh_dct[node], axis=0)
        seq = ''.join([bases[i] for i in argmaxes])
        node.label = seq

    return tree


def main():
    params = parse_args()
    genome = str(list(SeqIO.parse(params.genome, 'fasta'))[0].seq)
    coords = read_exact_positions(params.coords)
    rss = extract_rss(genome, coords, params.rss_len)
    tree = dendropy.Tree.get(path=params.tree_fp,
                             schema="newick",
                             rooting="default-rooted",
                             preserve_underscores="True")
    param_dct = read_param(params.param_fp)
    tree = get_internal_nodes(tree, rss, param_dct)
    tree.write(path=params.output, schema='newick')


if __name__ == "__main__":
    main()
