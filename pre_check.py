# -*- coding: utf-8 -*-
# @Time : 2020/7/11 22:44
# @Author : Zhongyi Hua
# @FileName: pre_check.py
# @Usage: chloroplast is a circle. This script will fix the start of sequence to trnH nearby
# @Note:
# @E-mail: njbxhzy@hotmail.com

"""
这个脚本有个前提：trnH要么被注释，要么正好在头尾gap上，且头尾gap的只有这一个feature（但是不排除IR区可能跨头尾，
真是头疼， 现在只能默认IR区不跨头尾，一个可能的处理方式是在PGA注释的时候给IR设置一个非常大的值，让IR注释不出来）
"""


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
import os
import re
import warnings


class FixPosition:
    def __init__(self, gb_path_pga, out_path):
        self.path = gb_path_pga
        self.out_path = out_path

    def _fix_postiion(self):
        """
        Note: This function give tacit consent to that trnH-GUG appears in genome
        :return:
        """
        print('fix position')
        with warnings.catch_warnings(record=True) as w:
            pga_seq = SeqIO.read(self.path, 'genbank')
            t_location = None
            t_start = None
            t_end = None
            for gene in pga_seq.features:
                if gene.qualifiers.get('gene', ['not_gene'])[0] == 'trnH-GUG' and gene.type == 'gene':
                    t_location = gene.location
            if t_location is None:
                for bio_war in w:
                    try:
                        t_start, t_end = re.search(r'(\d+)\.\.(\d+)', bio_war.message.__str__()).group(1, 2)
                    except AttributeError:
                        continue
                t_strand = -1
            else:
                t_start = t_location.start
                t_end = t_location.end
                t_strand = t_location.strand
            if t_strand == 1:
                part1, part2 = pga_seq.seq[0: int(t_end)+5], pga_seq.seq[int(t_end)+5:]
                seq_fixed = part2 + part1
                seq_fixed = SeqRecord(seq=seq_fixed, id=pga_seq.id, description='')
                seq_fixed.seq = seq_fixed.seq.reverse_complement()
            else:
                part1, part2 = pga_seq.seq[0: int(t_start)-5], pga_seq.seq[int(t_start)-5:]
                seq_fixed = part2 + part1
                seq_fixed = SeqRecord(seq=seq_fixed, id=pga_seq.id, description='')
            SeqIO.write(seq_fixed, self.out_path, 'fasta')
            print('fix done')

    def check_need(self):
        pga_seq = SeqIO.read(self.path, 'genbank')
        # check housekeeping gene
        gene_name_list = []
        for gene in pga_seq.features:
            gene_name_list.append(gene.qualifiers.get('gene', ['not_gene'])[0])
        if ('matK' not in gene_name_list) and ('matk' not in gene_name_list):
            print('matK loss!')
        if ('rbcL' not in gene_name_list) and ('rbcl' not in gene_name_list):
            print('rbcL loss!')
        if 'trnH-GUG' not in gene_name_list:
            print('trnH loss!')
            return
        # check and fix position
        try:
            pga_seq.features.sort(key=lambda x: x.location.start)
        except AttributeError:
            self._fix_postiion()
        else:
            if not pga_seq.features[1].qualifiers['gene'][0] == 'trnH-GUG':
                self._fix_postiion()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description='Check and fix the start of chloroplast genome')
    parser.add_argument('-i', '--info_table', required=True,
                        help='<file_path>  information table which has two columns: seq fasta path , PGA gff path, '
                             'output fasta path')
    args = parser.parse_args()
    info_table = pd.read_table(args.info_table, names=['gb_path_pga', 'out_path'])
    for ind, row in info_table.iterrows():
        print(os.path.basename(row.to_list()[0]))
        b = FixPosition(*row.to_list())
        b.check_need()
