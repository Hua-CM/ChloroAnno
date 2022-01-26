# -*- coding: utf-8 -*-
# @Time : 2020/7/6 15:56
# @Author : Zhongyi Hua
# @FileName: correct.py
# @Usage: correct annotation error using two gff result.
# @Note:
# @E-mail: njbxhzy@hotmail.com

import pandas as pd
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from utilities.check import CheckCp2, add_rps12, check_cds
from utilities.gbk2gff4PGA import gbk2gff
from BCBio import GFF
from collections import defaultdict
import os
import argparse


class CorrectGff2:

    def __init__(self, gff_path_ge, gff_path_pga, seq_ins, species_pre):
        self.geseq = next(GFF.parse(gff_path_ge))
        self.pga = next(GFF.parse(gff_path_pga))
        self.seq = seq_ins
        self.prefix = species_pre
        self.curated_features = []
        self.pga_dict = self._create_dialect_(self.pga)
        self.corrected = None

    @staticmethod
    def _create_dialect_(_annotation: SeqRecord):
        _dialect = defaultdict()
        for gene in _annotation.features:
            _dialect.setdefault(gene.qualifiers.get('Name')[0], []).append(gene)
        return _dialect

    @staticmethod
    def _cal_dis_(gene1, gene2):
        return abs((gene1.location.start+gene1.location.end)/2 - (gene2.location.start+gene2.location.end)/2)

    def _query_gene_(self, gene: SeqFeature):
        _query_name = gene.qualifiers.get('Name')[0]
        _subject = self.pga_dict.get(_query_name)
        if _subject:
            if len(_subject) > 1:
                dis_lst = [self._cal_dis_(gene, _gene) for _gene in _subject]
                index_min = min(range(len(dis_lst)), key=dis_lst.__getitem__)
                return _subject[index_min]
            else:
                return _subject[0]
        else:
            return

    def correct_records(self):
        for gene in self.geseq.features:
            gene_type = gene.qualifiers['gene_biotype'][0]
            if gene.qualifiers.get('Name')[0] == 'rps12':
                continue
            if gene_type == 'protein_coding':
                try:
                    check_cds(gene, self.seq)
                    self.curated_features.append(gene)
                except:
                    fixed_gene = self._query_gene_(gene)
                    if fixed_gene:
                        self.curated_features.append(fixed_gene)
            else:  # rRNA and tRNA
                self.curated_features.append(gene)
        self.curated_features.sort(key=lambda x: x.location.start)

    def add_trna(self):
        geseq_trnas = [gene for gene in self._create_dialect_(self.geseq).keys() if gene.startswith('trn')]
        for gene in self.pga_dict.keys():
            if gene.startswith('trn') and (gene not in geseq_trnas):
                self.curated_features += self.pga_dict.get(gene)
        self.curated_features.sort(key=lambda x: x.location.start)

    def check(self):
        self.corrected = SeqRecord(id=self.geseq.id, seq='', features=self.curated_features)
        checkIns = CheckCp2(self.corrected)
        checkIns.renumber(self.prefix)
        checkIns.check_region()
        checkIns.check_cds(self.seq, add_pseudo=True)


def parseArgs():
    # command line parser
    parser = argparse.ArgumentParser(
        description='Combine GeSeq and PGA annotation results ')

    subparsers = parser.add_subparsers(help='combine for combining GeSeq and PGA results; '
                                            'rps12 for adding rps12 and renumber')
    # for step6
    parser_a = subparsers.add_parser('combine', help='combine for combining GeSeq and PGA results')
    parser_a.add_argument('-i', '--info_table', required=True,
                          help='<file_path>  information table which has four columns: Geseq gff path, '
                               'PGA genbank path, locus prefix')
    parser_a.add_argument('-o', '--output', required=True,
                          help='<directory_path>  output directory')
    parser_a.set_defaults(subcmd='combine')
    # for step8
    parser_b = subparsers.add_parser('rps12', help='adding rps12 and renumber')
    parser_b.add_argument('-i', '--info_table', required=True,
                          help='<file_path>  information table which has three columns: '
                               'Gff path, PGA genbank file, locus prefix')
    parser_b.add_argument('-o', '--output', required=True,
                          help='<directory_path>  output directory')
    parser_b.set_defaults(subcmd='rps12')
    args = parser.parse_args()
    return args


def main(args):
    if args.subcmd == 'combine':
        info_table = pd.read_table(args.info_table,
                                   names=['gff_path_ge', 'gb_path_pga', 'species_pre'])
        for ind, row in info_table.iterrows():
            # use output path as a temporary curated path
            outfile = os.path.join(args.output, os.path.split(row['gff_path_ge'])[-1])
            genome_seq = gbk2gff(row['gb_path_pga'], outfile, row['species_pre'])
            main_ins = CorrectGff2(row['gff_path_ge'], outfile, genome_seq, row['species_pre'])
            main_ins.correct_records()
            main_ins.add_trna()
            main_ins.check()
            GFF.write([main_ins.corrected], open(outfile, 'w'), include_fasta=False)

    if args.subcmd == 'rps12':
        info_table = pd.read_table(args.info_table,
                                   names=['raw_gff_path', 'pga_gb_path', 'species_pre'])
        for ind, row in info_table.iterrows():
            print(os.path.basename(row['raw_gff_path']))
            outfile = os.path.join(args.output, os.path.split(row['raw_gff_path'])[-1])
            tmp_check = CheckCp2(row['raw_gff_path'])
            tmp_check.renumber(row['species_pre'])
            # use output path as temporary file
            GFF.write([tmp_check.gff], open(outfile, 'w'), include_fasta=False)
            try:
                tmp_df = add_rps12(row['pga_gb_path'], outfile, row['species_pre'])
                tmp_df.to_csv(outfile, sep='\t', index=False, header=False)
            except:
                print('please manually add rps12')


if __name__ == '__main__':
    main(parseArgs())
