# -*- coding: utf-8 -*-
# @Time : 2020/7/6 15:56
# @Author : Zhongyi Hua
# @FileName: correct.py
# @Usage: correct annotation error using two gff result.
# @Note:
# @E-mail: njbxhzy@hotmail.com

"""
问题：整合程度不够，目前没有整合PGA中有而GeSeq中没有的
思路： 根据locus直接插（插完需要检查）
问题： attribute自动加入pseudo=True
"""

import gffutils
import pandas as pd
import numpy as np
from Bio import SeqIO
from check import CheckCp
import os

get_seq = lambda seq, start, end: seq.seq[start - 1:end]


def flatten(ls):
    def faction(lis):
        for i in lis:
            if isinstance(i, list):
                faction(i)
            else:
                flatten_list.append(i)

    flatten_list = []
    faction(ls)
    return flatten_list


class CorrectGff:
    gff_fields = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]

    def __init__(self, gff_path_ge, gff_path_pga, seq_path, gff_new_path, seq_id, species_pre):
        self.species_pre = species_pre
        self.seq_id = seq_id
        self.seq = SeqIO.read(seq_path, 'fasta')
        self.gff = gffutils.create_db(gff_path_ge, ':memory:', merge_strategy='create_unique')
        self.gff_new_path = gff_new_path
        self.gff_pga = pd.read_table(gff_path_pga, comment='#', names=CorrectGff.gff_fields)
        # avoid duplicated index in geseq_dialect and pga_dialect
        self.gff_pga.index = self.gff_pga.index + 1000
        # correction will be conducted on this object
        self.gff_corrected = pd.read_table(gff_path_ge, comment='#', names=CorrectGff.gff_fields)
        # following objects would be used to record relationship
        self.geseq_dialect = pd.DataFrame()
        self.pga_dialect = pd.DataFrame()
        print('Start correct', os.path.basename(gff_path_ge))

    @staticmethod
    def get_parent(index, attributes):
        return index, {k.split('=')[0]: k.split('=')[1] for k in attributes.split(';')}['Parent']

    def create_geseq_dialect(self):
        # gene.id - index - child_index(list)
        self.geseq_dialect['gene_id'] = self.gff_corrected.loc[self.gff_corrected['type'] == 'gene', 'attributes']. \
            apply(lambda x: x.split(';')[0].split('=')[1])
        self.geseq_dialect['index'] = self.gff_corrected.index.to_series()
        tmp_df = self.gff_corrected.loc[self.gff_corrected['type'] == 'CDS'].\
            apply(lambda x: self.get_parent(x.name, x['attributes']), axis=1, result_type="expand").\
            groupby(1).agg(list).\
            reset_index().\
            rename(columns={0: 'child_list', 1: 'gene_id'})
        self.geseq_dialect = self.geseq_dialect.merge(tmp_df, how='right')
        self.geseq_dialect['locus'] = self.gff_corrected.loc[self.gff_corrected.index.isin(self.geseq_dialect['index']), 'attributes'].\
            apply(lambda x: {k.split('=')[0]: k.split('=')[1] for k in x.split(';')}['Name']).\
            reset_index(drop=True)

    def create_pga_dialect(self):
        # gene.id - index - locus - child_index(list)
        self.pga_dialect['gene_id'] = self.gff_pga.loc[self.gff_pga['type'] == 'gene', 'attributes']. \
            apply(lambda x: x.split(';')[0].split('=')[1])
        self.pga_dialect['index'] = self.gff_pga.index.to_series()
        tmp_df = self.gff_pga.loc[self.gff_pga['type'] == 'CDS'].\
            apply(lambda x: self.get_parent(x.name, x['attributes']), axis=1, result_type="expand").\
            groupby(1).agg(list).\
            reset_index().\
            rename(columns={0: 'child_list', 1: 'gene_id'})
        self.pga_dialect = self.pga_dialect.merge(tmp_df, how='right')
        self.pga_dialect['locus'] = self.gff_pga.loc[self.gff_pga.index.isin(self.pga_dialect['index']), 'attributes'].\
            apply(lambda x: {k.split('=')[0]: k.split('=')[1] for k in x.split(';')}['Name']).\
            reset_index(drop=True)

    def _correct_record(self, gene_feature):
        """
        correct a gene feature (with children)
        :param gene_feature:
        :return:
        """
        locus = gene_feature.attributes['Name'][0]
        drop_list = self.geseq_dialect.loc[self.geseq_dialect['gene_id'] == gene_feature.id, 'child_list'].to_list()
        drop_list.append(self.geseq_dialect.loc[self.geseq_dialect['gene_id'] == gene_feature.id, 'index'].to_list())
        drop_list = [k for x in drop_list for k in x]
        if len(self.pga_dialect.loc[self.pga_dialect['locus'] == locus]) > 1:
            tmp_dict = {}
            for _ in self.pga_dialect.loc[self.pga_dialect['locus'] == locus, 'index'].to_list():
                tmp_dict.update(
                    {self.gff_pga.loc[[_]].apply(lambda x: (x['start'] + x['end'])/2, axis=1).values[0]: _}
                )
            locus_idx = tmp_dict[
                min([_ for _ in tmp_dict.keys()], key=lambda x: abs((gene_feature.start + gene_feature.end)/2-x))
            ]
            gain_list = self.pga_dialect.loc[self.pga_dialect['index'] == locus_idx, 'child_list'].to_list()
            gain_list.append(self.pga_dialect.loc[self.pga_dialect['index'] == locus_idx, 'index'].to_list())
            gain_list = [k for x in gain_list for k in x]
        elif len(self.pga_dialect.loc[self.pga_dialect['locus'] == locus]) == 0:
            print('please check', locus, gene_feature.start, gene_feature.end, '(it is not in PGA)')
            return
        else:
            gain_list = self.pga_dialect.loc[self.pga_dialect['locus'] == locus, 'child_list'].to_list()
            gain_list.append(self.pga_dialect.loc[self.pga_dialect['locus'] == locus, 'index'].to_list())
            gain_list = [k for x in gain_list for k in x]
        self.gff_corrected.drop(drop_list, inplace=True)
        self.gff_corrected = self.gff_corrected.append(self.gff_pga.loc[gain_list])

    def _add_record(self, locus_ls):
        tmp_list = np.array(self.pga_dialect.loc[self.pga_dialect['locus'].isin(locus_ls), ['index', 'child_list']]).tolist()
        gain_list = flatten(tmp_list)
        self.gff_corrected = self.gff_corrected.append(self.gff_pga.loc[gain_list])

    def _add_pseudo(self, gene):
        self.gff_corrected.loc[self.gff_corrected['attributes'].str.startswith('ID='+gene.id), 'attributes'] += ';pseudo=true'

    def correct_gff(self):
        print('Auto correct start')
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
            try:
                seq_combined.translate(table=11, cds=True)
            except:
                self._correct_record(gene)
        print('Auto correct done')

    def add_pga(self):
        add_list = []
        for locus in self.pga_dialect['locus'].to_list():
            if locus not in self.geseq_dialect['locus'].to_list():
                add_list.append(locus)
        self._add_record(add_list)

    def renumber(self):
        self.gff_corrected.to_csv(self.gff_new_path, sep='\t', index=False, header=False)
        tmp_check = CheckCp(self.gff_new_path)
        self.gff_corrected = tmp_check.renumber(self.species_pre, self.seq_id)
        self.gff_corrected.to_csv(self.gff_new_path, sep='\t', index=False, header=False)

    def check(self):
        print('Auto check start')
        gff_new = gffutils.create_db(self.gff_new_path, ':memory:', merge_strategy='create_unique')
        for gene in gff_new.features_of_type('gene', order_by='start'):
            seq_combined = ""
            strand = ""
            for cds in gff_new.children(gene, featuretype='CDS', order_by='start'):
                seq = get_seq(self.seq, cds.start, cds.end)
                strand = cds.strand
                seq_combined += seq
            if seq_combined == '':
                continue
            elif strand == '-':
                seq_combined = seq_combined.reverse_complement()
            try:
                seq_combined.translate(table=11, cds=True).rstrip('*')
            except Exception as e:
                print(gene.id, gene.attributes['Name'][0])
                print(e)
                self._add_pseudo(gene)
        self.gff_corrected.to_csv(self.gff_new_path, sep='\t', index=False, header=False)
        print('Auto check done')

    def check2(self):
        tmp_check = CheckCp(self.gff_new_path)
        tmp_check.check_region()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description='Combine GeSeq and PGA annotation results ')
    parser.add_argument('-i', '--info_table', required=True,
                        help='<file_path>  information table which has four columns: Geseq gff path, '
                             'result path, seqid, locus prefix')
    args = parser.parse_args()
    info_table = pd.read_table(args.info_table, names=['gff_path_ge', 'gff_path_pga', 'seq_path', 'gff_new_path', 'seq_id', 'species_pre'])
    for ind, row in info_table.iterrows():
        a = CorrectGff(*row.to_list())
        a.create_geseq_dialect()
        a.create_pga_dialect()
        a.correct_gff()
        a.add_pga()
        a.renumber()
        a.check()
        a.check2()
