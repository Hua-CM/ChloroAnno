# -*- coding: utf-8 -*-
# @Time : 2020/7/6 15:56
# @Author : Zhongyi Hua
# @FileName: correct.py
# @Usage: correct annotation error using two gff result.
# @Note:
# @E-mail: njbxhzy@hotmail.com

"""
思路：根据locus匹配PGA gff中对应的feature进行替换。主要解决3倍数问题和stop codon问题。采用整体替换的方法(包括children)。
考虑到gffutils为sqlite数据库，不适合修改，故考虑使用gffutils查询，pandas修改的方式
多拷贝的直接报错手改
"""

import gffutils
import pandas as pd
from Bio import SeqIO
from gff2gff4GeSeq import get_record
import os

get_seq = lambda seq, start, end: seq.seq[start - 1:end]


class CorrectGff:
    gff_fields = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]

    def __init__(self, gff_path_ge, gff_path_pga, seq_path, gff_new_path, species_pre, seq_id):
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
            print('please check', locus, gene_feature.start, gene_feature.end)
            return
        else:
            gain_list = self.pga_dialect.loc[self.pga_dialect['locus'] == locus, 'child_list'].to_list()
            gain_list.append(self.pga_dialect.loc[self.pga_dialect['locus'] == locus, 'index'].to_list())
            gain_list = [k for x in gain_list for k in x]
        self.gff_corrected.drop(drop_list, inplace=True)
        self.gff_corrected = self.gff_corrected.append(self.gff_pga.loc[gain_list])

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
            if ('*' in seq_combined.translate(table=11).rstrip('*')) or (not len(seq_combined) % 3 == 0):
                self._correct_record(gene)
        print('Auto correct done')

    def renumber(self):
        print('Renumber gene id')
        self.gff_corrected.to_csv(self.gff_new_path, sep='\t', index=False, header=False)
        gff_new = gffutils.create_db(self.gff_new_path, ':memory:', merge_strategy='create_unique')
        feature_list = []
        gene_count = 0
        for gene in gff_new.features_of_type('gene', order_by='start'):
            gene_count += 1
            gene_id = self.species_pre + '%03d' % gene_count
            gene_name = gene.attributes['Name'][0]
            gene_type = gene.attributes['gene_biotype'][0]
            gene_attributes = ['ID=' + gene_id,
                               'Name=' + gene_name,
                               'gene_biotype=' + gene_type
                               ]
            feature_list.append(get_record(gene, 'gene', gene_attributes))
            child_count = 0
            if gene_type == 'protein_coding':
                for cds in gff_new.children(gene, featuretype='CDS', order_by='start'):
                    child_count += 1
                    cds_attributes = ['ID=' + 'cds_' + gene_id + '_' + str(child_count),
                                      'Parent=' + gene_id,
                                      'product=' + cds.attributes['product'][0]]
                    cds_record = get_record(cds, 'CDS', cds_attributes)
                    cds_record.update({'phase': cds.frame})
                    feature_list.append(cds_record)
            else:
                for rna in gff_new.children(gene, featuretype=('tRNA', 'rRNA')):
                    rna_attributes = ['ID=' + 'rna_' + gene_id + '_1',
                                      'Parent=' + gene_id,
                                      'product=' + rna.attributes['product'][0]]
                    feature_list.append(get_record(gene, gene_type, rna_attributes))
                for exon in gff_new.children(gene, featuretype='exon', order_by='start'):
                    child_count += 1
                    exon_attributes = ['ID=' + 'exon_' + gene_id + '_' + str(child_count),
                                       'Parent=' 'rna_' + gene_id + '_1'
                                       ]
                    feature_list.append(get_record(exon, 'exon', exon_attributes))
        result_gff = pd.DataFrame.from_dict({index: record for index, record in enumerate(feature_list)}, 'index')
        result_gff['seqid'] = self.seq_id
        result_gff['score'] = '.'
        result_gff['source'] = 'GeSeq'
        result_gff = result_gff[["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]]
        result_gff.to_csv(self.gff_new_path, sep='\t', index=False, header=False)
        print('Renumber done')

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
                print(gene.id)
                print(e)
        print('Auto check done')
