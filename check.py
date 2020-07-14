# -*- coding: utf-8 -*-
# @Time : 2020/7/13 0:38
# @Author : Zhongyi Hua
# @FileName: check.py
# @Usage: 
# @Note:
# @E-mail: njbxhzy@hotmail.com

"""
把各个部分的check任务全部拆出来，这样方便反复使用
1. 保守基因情况
2. un normal name
3. 重复区域
4. renumber
5. cds check
"""

import gffutils
import portion as pt
import pandas as pd
from Bio import SeqIO

get_seq = lambda seq, start, end: seq.seq[start - 1:end]

tb_cds = pd.read_table('cds.txt')
tb_rna = pd.read_table('rna.txt')

standard_name_list = [record['name'] for record in tb_cds.to_dict('records')]
standard_name_list += [record['name'] for record in tb_rna.to_dict('records')]


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
        for gene in self.gff.features_of_type('Name', order_by='start'):
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

    def renumber(self, species_pre, seq_id):
        print('Renumber gene id')
        feature_list = []
        gene_count = 0
        for gene in self.gff.features_of_type('gene', order_by='start'):
            gene_count += 1
            gene_id = species_pre + '%03d' % gene_count
            gene_name = gene.attributes['Name'][0]
            gene_type = gene.attributes['gene_biotype'][0]
            gene_attributes = ['ID=' + gene_id,
                               'Name=' + gene_name,
                               'gene_biotype=' + gene_type
                               ]
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
        result_gff = pd.DataFrame.from_dict({index: record for index, record in enumerate(feature_list)}, 'index')
        result_gff['seqid'] = seq_id
        result_gff['score'] = '.'
        result_gff['source'] = 'GeSeq'
        result_gff = result_gff[["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]]
        print('Renumber done')
        return result_gff

    def check_cds(self, seq_path):
        geo_seq = SeqIO.read(seq_path, 'fasta')
        print('Auto check start')
        for gene in self.gff.features_of_type('gene', order_by='start'):
            seq_combined = ""
            strand = ""
            for cds in self.gff.children(gene, featuretype='CDS', order_by='start'):
                seq = get_seq(geo_seq, cds.start, cds.end)
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
