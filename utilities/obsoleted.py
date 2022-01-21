# -*- coding: utf-8 -*-
# @Time    : 2022/1/20 20:59
# @Author  : Zhongyi Hua
# @File    : obsoleted.py
# @Usage   : Just store some obsoleted codes
# @Note    : 
# @E-mail  : njbxhzy@hotmail.com

import gffutils
import portion as pt
import numpy as np
from check import standard_name_list
import pandas as pd
from Bio import SeqIO

get_seq = lambda seq, start, end: seq.seq[start - 1:end]


def get_record(record_feature, record_type, attributes):
    """
    Transform a gff feature of to a dict
    :param record_feature: a record feature (Class gffutils.Feature)
    :param record_type: gene/CDS/tRNA/rRNA...
    :param attributes: a list of character that positioned in gff3 column9
    :return: gff_df
    """
    feature_record = {'type': record_type,
                      'start': record_feature.start,
                      'end': record_feature.end,
                      'strand': record_feature.strand,
                      'phase': '.',
                      'attributes': ";".join(attributes)}
    return feature_record


class CheckCp:
    standard_name_list = standard_name_list

    def __init__(self, gff_path):
        self.gff = gffutils.create_db(gff_path, ':memory:', merge_strategy='create_unique')

    def check_hs(self):
        """
        Check housekeeping gene: matK, rbcL
        :return: a message
        """
        gene_name_list = []
        for gene in self.gff.features_of_type('gene', order_by='start'):
            gene_name_list.append(gene.attributes['Name'][0])
        if ('matK' not in gene_name_list) and ('matk' not in gene_name_list):
            print('matK loss!')
        if ('rbcL' not in gene_name_list) and ('rbcl' not in gene_name_list):
            print('rbcL loss!')

    def check_name(self):
        """
        Check unstandard name
        :return: messages
        """
        gene_name_list = []
        for gene in self.gff.features_of_type('gene', order_by='start'):
            gene_name_list.append(gene.attributes['Name'][0])
        for gene_name in gene_name_list:
            if not ((gene_name in CheckCp.standard_name_list) or gene_name.startswith('orf')):
                print('check ' + gene_name)

    def check_region(self):
        """
        check duplicated gene region
        :return: messages
        """
        region_list = []
        locus_list = []
        for gene in self.gff.features_of_type('gene', order_by='start'):
            region_list.append(pt.closed(gene.start, gene.end))
            locus_list.append([gene.attributes['Name'][0]])
        for i in range(len(region_list) - 1):
            if not region_list[i] < region_list[i + 1]:
                print(locus_list[i], region_list[i], ' and ', locus_list[i + 1], region_list[i + 1],
                      ' are duplicated')
        print('check duplicated region done')

    def renumber(self, seq_id, species_pre):
        print('Renumber gene id')
        feature_list = []
        gene_count = 0
        for gene in self.gff.features_of_type('gene', order_by='start'):
            gene_count += 1
            gene_id = species_pre + '%03d' % gene_count
            gene_attributes = ['ID=' + gene_id]
            gene_type = gene.attributes['gene_biotype'][0]
            gene_attributes += [_[0] + '=' + _[1][0] for _ in gene.attributes.items() if not _[0] == 'ID']
            feature_list.append(get_record(gene, 'gene', gene_attributes))
            child_count = 0
            if gene_type == 'protein_coding':
                for cds in self.gff.children(gene, featuretype='CDS', order_by='start'):
                    child_count += 1
                    cds_attributes = ['ID=' + 'cds_' + gene_id + '_' + str(child_count),
                                      'Parent=' + gene_id,
                                      'product=' + cds.attributes['product'][0]]
                    cds_record = get_record(cds, 'CDS', cds_attributes)
                    cds_record.update({'phase': cds.frame})
                    feature_list.append(cds_record)
            else:
                for rna in self.gff.children(gene, featuretype=('tRNA', 'rRNA')):
                    rna_attributes = ['ID=' + 'rna_' + gene_id + '_1',
                                      'Parent=' + gene_id,
                                      'product=' + rna.attributes['product'][0]]
                    feature_list.append(get_record(gene, gene_type, rna_attributes))
                for exon in self.gff.children(gene, featuretype='exon', order_by='start'):
                    child_count += 1
                    exon_attributes = ['ID=' + 'exon_' + gene_id + '_' + str(child_count),
                                       'Parent=' 'rna_' + gene_id + '_1'
                                       ]
                    feature_list.append(get_record(exon, 'exon', exon_attributes))

        _result_gff = pd.DataFrame(feature_list)
        _result_gff['seqid'] = seq_id
        _result_gff['score'] = '.'
        _result_gff['source'] = 'GeSeq'
        _result_gff = _result_gff[["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]]
        print('Renumber done')
        return _result_gff

    def check_cds(self, seq_path):
        geo_seq = SeqIO.read(seq_path, 'fasta')
        print('Auto check start')
        for gene in self.gff.features_of_type('gene', order_by='start'):
            if gene.attributes['gene_biotype'] == ['protein_coding']:
                if gene.attributes['Name'] == ['rps12']:
                    continue
                seq_combined = ""
                cds_count = 0
                for cds in self.gff.children(gene,
                                             featuretype='CDS',
                                             order_by='start',
                                             reverse=False if gene.strand == '+' else True):

                    if (cds_count == 0) and (not cds.frame == '0'):
                        print('check ', gene.id, ' frame')
                    seq = get_seq(geo_seq, cds.start, cds.end)
                    if cds.strand == '-':
                        seq_combined += seq.reverse_complement()
                    else:
                        seq_combined += seq
                    cds_count += 1
                if seq_combined == '':
                    continue
                elif seq_combined.__len__() <= 33:
                    print('The CDS length of', gene.id, 'is less than 33 bp')
                try:
                    seq_combined.translate(table=11, cds=True)
                except Exception as e:
                    print(gene.id)
                    print(e)
        print('Auto check done')

def add_rps12(geseq_gff, new_gff, seq_path, species_pre):
    """

    :param geseq_gff: raw geseq gff file
    :param new_gff: renumbered gff
    :param species_pre:
    :return:
    """
    if type(new_gff) == str:
        new_gff = pd.read_table(new_gff,
                                comment='#',
                                names=["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])
    raw_gff = pd.read_table(geseq_gff,
                            comment='#',
                            names=["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])
    seq = SeqIO.read(seq_path, 'fasta')
    gene_count = sum(new_gff['type'] == 'gene')
    # get rps12
    features_list = []
    gene_features = raw_gff[(raw_gff['attributes'].str.contains('rps12')) & (raw_gff['type'] == 'gene')]
    part1 = gene_features[gene_features.duplicated(subset=['start', 'end'])].iloc[0]
    part2 = gene_features.drop_duplicates(keep=False)
    err_list = []
    for idx, part in part2.iterrows():
        gene_count += 1
        gene_id = species_pre + '%03d' % gene_count
        part_attributes = ['ID=' + gene_id,
                           'Name=rps12',
                           'gene_biotype=protein_coding'
                           ]
        cds1 = raw_gff[(raw_gff['start'] == part.start) & (raw_gff['type'] == 'exon')].iloc[0]
        cds2 = raw_gff[(raw_gff['end'] == part.end) & (raw_gff['type'] == 'exon')].iloc[0]
        cds1_attributes = ['ID=' + 'cds_' + gene_id + '_1',
                           'Parent=' + gene_id,
                           'product=30S ribosomal protein S12']
        cds2_attributes = ['ID=' + 'cds_' + gene_id + '_2',
                           'Parent=' + gene_id,
                           'product=30S ribosomal protein S12']
        cds3_attributes = ['ID=' + 'cds_' + gene_id + '_3',
                           'Parent=' + gene_id,
                           'product=30S ribosomal protein S12']
        features_list.append(get_record(part1, 'gene', part_attributes + ['part=1/2']))
        features_list.append(get_record(part, 'gene', part_attributes + ['part=2/2']))
        features_list.append(get_record(part1, 'CDS', cds1_attributes))
        features_list.append(get_record(cds1, 'CDS', cds2_attributes))
        features_list.append(get_record(cds2, 'CDS', cds3_attributes))
        # check translation
        seq_part1 = get_seq(seq, part1.start, part1.end)
        if part1.strand == '-':
            seq_part1 = seq_part1.reverse_complement()
        seq_combined = ''
        for feature in [cds1, cds2]:
            seq_combined += get_seq(seq, feature.start, feature.end)
        if cds1.strand == '-':
            seq_combined.reverse_complement()
        seq_combined = seq_part1 + seq_combined
        try:
            seq_combined.translate(table=11, cds=True)
        except Exception as e:
            print(gene_id)
            print(e)
            err_list.append(gene_id)
    rps12_df = pd.DataFrame.from_dict({idx: feature for idx, feature in enumerate(features_list)}, 'index')
    rps12_df.loc[rps12_df['type'] == 'CDS', 'phase'] = 0
    rps12_df['seqid'] = new_gff.seqid.to_list()[0]
    rps12_df['source'] = 'GeSeq'
    rps12_df['score'] = '.'
    for err_id in err_list:
        rps12_df.loc[rps12_df['attributes'].str.startswith('ID='+err_id), 'attributes'] += ';pseudo=true'
    new_gff = new_gff.append(rps12_df, sort=False)
    return new_gff


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
        self.geseq_dialect = self.gff_corrected.loc[self.gff_corrected.index.isin(self.geseq_dialect['index']), 'attributes'].\
            apply(lambda x: {k.split('=')[0]: k.split('=')[1] for k in x.split(';')}['Name']).\
            reset_index().\
            rename(columns={'attributes': 'locus'}).\
            merge(self.geseq_dialect)


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
        self.pga_dialect = self.gff_pga.loc[self.gff_pga.index.isin(self.pga_dialect['index']), 'attributes'].\
            apply(lambda x: {k.split('=')[0]: k.split('=')[1] for k in x.split(';')}['Name']).\
            reset_index().\
            rename(columns={'attributes': 'locus'}).\
            merge(self.pga_dialect)

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