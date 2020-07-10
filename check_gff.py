# -*- coding: utf-8 -*-
# @Time : 2020/7/7 0:00
# @Author : Zhongyi Hua
# @FileName: check_gff.py
# @Usage: check gff error
# @Note: Please use redirect.py to redirect gene id after manual curation
# @E-mail: njbxhzy@hotmail.com

import gffutils
import portion as pt
import pandas as pd
from Bio import SeqIO
import os

get_seq = lambda seq, start, end: seq.seq[start - 1:end]


class CheckGff:
    def __init__(self, gff_path, seq_path):
        self.gff = gffutils.create_db(gff_path, ':memory:', merge_strategy='create_unique')
        self.seq = SeqIO.read(seq_path, 'fasta')
        print('Start change', os.path.basename(gff_path))

    def check_duplicate_region(self):
        region_dict = {}
        for gene in self.gff.features_of_type('gene', order_by='start'):
            region_dict[gene.id] = pt.closed(gene.start, gene.end)
        region_list = list(region_dict.values())
        id_list = list(region_dict.keys())
        for i in range(len(region_dict)-1):
            if not region_list[i] < region_list[i+1]:
                print(id_list[i], region_list[i], ' and ', id_list[i+1], region_list[i+1], ' are duplicated')
        print('check duplicated region done')

    def check_triple(self):
        for gene in self.gff.features_of_type('gene', order_by='start'):
            cds_length = 0
            for child in self.gff.children(gene, featuretype='CDS'):
                cds_length += (child.end - child.start + 1)
            if not cds_length % 3 == 0:
                print(gene.id, "'s length is not a multiple of three")

    def check_stop_codon(self):
        for gene in self.gff.features_of_type('gene', order_by='start'):
            seq_combined = ""
            strand = ""
            for cds in self.gff.children(gene, featuretype='CDS', order_by='start'):
                seq = get_seq(self.seq, cds.start, cds.end)
                strand = cds.strand
                seq_combined += seq
            if seq_combined == '':
                continue
            elif strand == '-':
                seq_combined = seq_combined.reverse_complement()
            if '*' in seq_combined.translate(table=11).rstrip('*'):
                print(gene.id, ' has stop codon')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='This script is used to check error in gff')
    parser.add_argument('-i', '--info_table', required=True,
                        help='<file_path>  information table which has two columns: gff path: '
                              'genome sequence path (in fasta format)')
    args = parser.parse_args()
    info_table = pd.read_table(args.info_table)
    for ind, row in info_table.iterrows():
        a = CheckGff(*row.to_list())
        a.check_duplicate_region()
        a.check_triple()
        a.check_stop_codon()
