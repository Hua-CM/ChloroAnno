# -*- coding: utf-8 -*-
# @Time : 2020/7/11 22:44
# @Author : Zhongyi Hua
# @FileName: pre_check.py
# @Usage: chloroplast is a circle. This script will fix the start of sequence to trnH nearby
# @Note:
# @E-mail: njbxhzy@hotmail.com

"""
现在的脚本没法处理两个trnH的情况，准备补充下trnH前后有psbA的情况作为进一步判断
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import re
from tempfile import mktemp
from shutil import copy
from Bio.Blast.Applications import NcbitblastnCommandline, NcbiblastnCommandline, NcbimakeblastdbCommandline
from utilities.check import del_dir
from collections import defaultdict


class CurateSequence:
    def __init__(self, _seq_path, _out_path):
        self.genome = _seq_path
        self.query1 = os.path.abspath(os.path.join(__file__, '../ref/trnH.fa'))
        self.query2 = os.path.abspath(os.path.join(__file__, '../ref/protein.fa'))
        self.output = _out_path
        self.tmp = mktemp()

    def _make_db(self):
        os.mkdir(self.tmp)
        copy(self.genome, os.path.join(self.tmp, 'genome.fa'))
        self.genome = os.path.join(self.tmp, 'genome.fa')
        cline = NcbimakeblastdbCommandline(input_file=self.genome, dbtype='nucl')
        cline()

    def _blastn_wrapper(self):
        cline = NcbiblastnCommandline(
            query=self.query1,
            db=self.genome,
            evalue=0.001,
            out=os.path.join(self.tmp, 'blastn.res'),
            outfmt="6 qseqid sstart send sstrand evalue")
        cline()

    def _tblastn_wrapper(self):
        cline = NcbitblastnCommandline(
            query=self.query2,
            db=self.genome,
            evalue=0.001,
            out=os.path.join(self.tmp, 'tblastn.res'),
            outfmt="6 qseqid sstart send sstrand evalue")
        cline()

    def curate(self):
        self._make_db()
        self._blastn_wrapper()
        self._tblastn_wrapper()

        def _cal_dis_(gene1, gene2):
            return abs(
                (int(gene1.get('start')) + int(gene1.get('end'))) / 2 - (int(gene2.get('start')) + int(gene2.get('end'))) / 2
            )

        raw_seq = SeqIO.read(self.genome, 'fasta')
        key_lst = 'start end strand evalue'.split()
        res_dict = defaultdict()
        with open(os.path.join(self.tmp, 'blastn.res')) as f1, open(os.path.join(self.tmp, 'tblastn.res')) as f2:
            for line in f1.read().splitlines() + f2.read().splitlines():
                res_dict.setdefault(line.split()[0], []).append(dict(zip(key_lst, line.split()[1:])))
        # check housekeeping
        if ('matK' not in res_dict) and ('matk' not in res_dict):
            print('matK loss!')
        if ('rbcL' not in res_dict) and ('rbcl' not in res_dict):
            print('rbcL loss!')
        # check ambiguous nucleotide
        if re.compile('[^ATCGNatcgn]').findall(str(raw_seq.seq)):
            print('Genome contain invalid characters (not ATCGNatcgn)')
        # fix position
        ## for multiple psbA
        index_min = min(range(len(res_dict.get('psbA'))), key=lambda x: float(res_dict.get('psbA')[x].get('evalue')))
        psbA = res_dict.get('psbA')[index_min]
        ## for multiple trnH
        dis_lst = [_cal_dis_(_gene, psbA) for _gene in res_dict.get('trnH-GUG')]
        index_min = min(range(len(dis_lst)), key=dis_lst.__getitem__)
        trnh = res_dict.get('trnH-GUG')[index_min]
        if trnh.get('strand') == 'plus':
            _start = len(raw_seq) - int(trnh.get('start')) + 1 - 5
            raw_seq.seq = raw_seq.seq.reverse_complement()
        else:
            _start = int(trnh.get('start'))
        part1, part2 = raw_seq.seq[0: _start], raw_seq.seq[_start:]
        seq_fixed = part2 + part1
        seq_fixed = SeqRecord(seq=seq_fixed, id=raw_seq.id, description='')
        SeqIO.write(seq_fixed, self.output, 'fasta')
        del_dir(self.tmp)


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='Check and fix the start of chloroplast genome')
    parser.add_argument('-i', '--input', required=True,
                        help='<file_path> one sequence path per line')
    parser.add_argument('-o', '--output', required=True,
                        help='<directory>  output directory')
    args = parser.parse_args()
    with open(args.input) as f_in:
        sequence_lst = f_in.read().splitlines()
    for row in sequence_lst:
        print(os.path.basename(row))
        main_ins = CurateSequence(row, os.path.join(args.output, os.path.basename(row)))
        main_ins.curate()


if __name__ == '__main__':
    main()
