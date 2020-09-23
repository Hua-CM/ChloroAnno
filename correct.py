# -*- coding: utf-8 -*-
# @Time : 2020/7/6 15:56
# @Author : Zhongyi Hua
# @FileName: correct.py
# @Usage: correct annotation error using two gff result.
# @Note:
# @E-mail: njbxhzy@hotmail.com

import gffutils
import pandas as pd
import numpy as np
from Bio import SeqIO
from check import CheckCp, add2_rps12
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
        # gene.id - locus - index - child_index(list)
        self.geseq_dialect['gene_id'] = self.gff_corrected.loc[self.gff_corrected['type'] == 'gene', 'attributes']. \
            apply(lambda x: x.split(';')[0].split('=')[1])
        self.geseq_dialect['index'] = self.gff_corrected.index.to_series()
        tmp_df = self.gff_corrected.loc[self.gff_corrected['type'] == 'CDS'].\
            apply(lambda x: self.get_parent(x.name, x['attributes']), axis=1, result_type="expand").\
            groupby(1).agg(list).\
            reset_index().\
            rename(columns={0: 'child_list', 1: 'gene_id'})
        try:
            tmp_df2 = self.gff_corrected.loc[self.gff_corrected['type'].isin(['tRNA', 'exon'])].\
                apply(lambda x: self.get_parent(x.name, x['attributes']), axis=1, result_type="expand")
            tmp_df2[1] = tmp_df2.apply(lambda x: str.replace(x[1], 'rna_', ''), axis=1)
            tmp_df2 = tmp_df2.\
                groupby(1).\
                agg(list).\
                reset_index().\
                rename(columns={0: 'child_list', 1: 'gene_id'})
            tmp_df = tmp_df.append(tmp_df2, ignore_index=True)
        except:
            pass
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
        tmp_df2 = self.gff_pga.loc[self.gff_pga['type'].isin(['tRNA', 'exon'])].\
            apply(lambda x: self.get_parent(x.name, x['attributes']), axis=1, result_type="expand")
        tmp_df2[1] = tmp_df2.apply(lambda x: str.replace(x[1], 'rna_', ''), axis=1)
        tmp_df2 = tmp_df2.\
            groupby(1).\
            agg(list).\
            reset_index().\
            rename(columns={0: 'child_list', 1: 'gene_id'})
        tmp_df =tmp_df.append(tmp_df2, ignore_index=True)
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

    def _find_start_codon(self, gene, _seq_combined):

        def _find_start(_seq, _count):
            _seq = _seq[3:]
            _count += 1
            if _seq == '' or _count > 9:
                return None
            try:
                _seq.translate(cds=True, table=11)
                return _count*3
            except:
                _find_start(_seq, _count)

        shift_num = _find_start(_seq_combined, 0)
        if shift_num is None:
            self._add_pseudo(gene)
            return print(gene.id, 'Maybe RNA-editing, but no alternative start codon find')
        else:
            print(gene.id, 'Mayebe RNA-editing, and a alternative start codon found')
            if gene.strand == '-':
                self.gff_corrected.loc[self.gff_corrected['attributes'].str.startswith('ID=' + gene.id), 'end'] \
                 = int(self.gff_corrected.loc[self.gff_corrected['attributes'].str.startswith('ID=' + gene.id),
                                              'end']) - 3*shift_num
            else:
                self.gff_corrected.loc[self.gff_corrected['attributes'].str.startswith('ID=' + gene.id), 'start'] \
                 = int(self.gff_corrected.loc[self.gff_corrected['attributes'].str.startswith('ID=' + gene.id),
                                              'start']) + 3*shift_num

    def correct_gff(self):
        print('Auto correct start')
        for gene in self.gff.features_of_type('gene', order_by='start'):
            if gene.attributes['gene_biotype'][0] == 'protein_coding':
                seq_combined = ""
                for cds in self.gff.children(gene,
                                             featuretype='CDS',
                                             order_by='start',
                                             reverse=False if gene.strand == '+' else True):
                    seq = get_seq(self.seq, cds.start, cds.end)
                    if cds.strand == '-':
                        seq_combined += seq.reverse_complement()
                    else:
                        seq_combined += seq
                try:
                    seq_combined.translate(cds=True, table=11)
                except:
                    self._correct_record(gene)
            elif gene.attributes['gene_biotype'][0] == 'tRNA':
                ge_st, ge_ed = gene.start, gene.end
                idx_list = self.pga_dialect.loc[self.pga_dialect['locus'] == gene.attributes['Name'][0], 'index'].to_list()
                if idx_list:
                    pga_st_list = self.gff_pga.loc[idx_list, 'start'].to_list()
                    pga_ed_list = self.gff_pga.loc[idx_list, 'end'].to_list()
                    if (ge_st not in pga_st_list) or (ge_ed not in pga_ed_list):
                        self._correct_record(gene)
        print('Auto correct done')

    def add_trna(self):
        add_list = []
        for locus in [x for x in self.pga_dialect['locus'].to_list() if x.startswith('trn')]:
            if locus not in self.geseq_dialect['locus'].to_list():
                add_list.append(locus)
        self._add_record(add_list)

    def renumber(self):
        self.gff_corrected.to_csv(self.gff_new_path, sep='\t', index=False, header=False)
        tmp_check = CheckCp(self.gff_new_path)
        self.gff_corrected = tmp_check.renumber(self.seq_id, self.species_pre)
        self.gff_corrected.to_csv(self.gff_new_path, sep='\t', index=False, header=False)

    def check(self):
        print('Auto check start')
        gff_new = gffutils.create_db(self.gff_new_path, ':memory:', merge_strategy='create_unique')
        for gene in gff_new.features_of_type('gene', order_by='start'):
            if gene.attributes['Name'] == ['rps12']:
                continue
            seq_combined = ""
            for cds in gff_new.children(gene,
                                        featuretype='CDS',
                                        order_by='start',
                                        reverse=False if gene.strand == '+' else True):
                seq = get_seq(self.seq, cds.start, cds.end)
                if cds.strand == '-':
                    seq_combined += seq.reverse_complement()
                else:
                    seq_combined += seq
            if seq_combined == '':
                continue
            elif seq_combined.__len__() <= 33:
                print('The CDS length of', gene.id, 'is less than 33 bp')
            try:
                seq_combined.translate(table=11, cds=True)
            except Exception as e:
                if e.__str__() == "First codon 'ACG' is not a start codon":
                    self._find_start_codon(gene, seq_combined)
                else:
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

    # wrapper for sub-command line
    def combine(_args):
        info_table = pd.read_table(_args.info_table,
                                   names=['gff_path_ge', 'gff_path_pga', 'seq_path', 'gff_new_path', 'seq_id',
                                          'species_pre'])
        for ind, row in info_table.iterrows():
            a = CorrectGff(*row.to_list())
            a.create_geseq_dialect()
            a.create_pga_dialect()
            a.correct_gff()
            a.add_trna()
            a.renumber()
            a.check()
            a.check2()

    def rps12(_args):
        info_table = pd.read_table(_args.info_table,
                                   names=['raw_gff_path', 'new_gff_path', 'pga_gb_path',
                                          'seq_id', 'species_pre'])
        for ind, row in info_table.iterrows():
            print(os.path.basename(row['raw_gff_path']))
            tmp_check = CheckCp(row['raw_gff_path'])
            tmp_df = tmp_check.renumber(*row[['seq_id', 'species_pre']])
            try:
                tmp_df = add2_rps12(row['pga_gb_path'], tmp_df, row['species_pre'])
            except:
                print('please manually add rps12')
            finally:
                tmp_df.to_csv(row['new_gff_path'], sep='\t', index=False, header=False)

    # command line parser
    parser = argparse.ArgumentParser(
        description='Combine GeSeq and PGA annotation results ')

    subparsers = parser.add_subparsers(help='combine for combining GeSeq and PGA results; '
                                            'rps12 for adding rps12 and renumber')
    # for step6
    parser_a = subparsers.add_parser('combine', help='add help')
    parser_a.add_argument('-i', '--info_table', required=True,
                          help='<file_path>  information table which has four columns: Geseq gff path, '
                               'result path, seqid, locus prefix')
    parser_a.set_defaults(func=combine)
    # for step8
    parser_b = subparsers.add_parser('rps12', help='add help')
    parser_b.add_argument('-i', '--info_table', required=True,
                          help='<file_path>  information table which has four columns: Geseq gff path, '
                               'result path, seqid, locus prefix')
    parser_b.set_defaults(func=rps12)
    args = parser.parse_args()
    args.func(args)
