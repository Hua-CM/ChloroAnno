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
import os

get_seq = lambda seq, start, end: seq.seq[start - 1:end]

location = os.path.abspath(os.path.join(__file__, '..'))
tb_cds = pd.read_table(os.path.join(location, 'ref/cds.txt'))
tb_rna = pd.read_table(os.path.join(location, 'ref/rna.txt'))

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


def add2_rps12(pga_gb, new_gff, species_pre):
    def _get_record(location, fea_type, attributes):
        return {'type': fea_type,
                'start': location.start+1,
                'end': location.end,
                'strand': '-' if location.strand == -1 else '+',
                'phase': '.',
                'attributes': ";".join(attributes)}
    if type(new_gff) == str:
        new_gff = pd.read_table(new_gff,
                                comment='#',
                                names=["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])
    genome = SeqIO.read(pga_gb, 'genbank')
    features_list = []
    err_list = []
    gene_count = sum(new_gff['type'] == 'gene')
    rps12_list = [ele for ele in genome.features if ele.type == 'CDS' and ele.qualifiers.get('gene') == ['rps12']]
    part1 = [part for part in rps12_list if len(part.location.parts) == 1][0]
    part2_list = [part for part in rps12_list if len(part.location.parts) > 1]
    for part in part2_list:
        gene_count += 1
        gene_id = species_pre + '%03d' % gene_count
        part_attributes = ['ID=' + gene_id,
                           'Name=rps12',
                           'gene_biotype=protein_coding'
                           ]
        features_list.append(_get_record(part1.location, 'gene', part_attributes + [
            'Exception=trans-splicing;part=1/2']))
        features_list.append(_get_record(part.location, 'gene', part_attributes + [
            'Exception=trans-splicing;part=2/2']))
        cds_count = 1
        cds_attributes = ['ID=' + 'cds_' + gene_id + '_' + str(cds_count),
                          'Parent=' + gene_id,
                          'product=30S ribosomal protein S12']
        features_list.append(_get_record(part1.location, 'CDS', cds_attributes))
        seq_part1 = get_seq(genome, part1.location.start+1, part1.location.end)
        if part1.location.strand == -1:
            seq_part1 = seq_part1.reverse_complement()
        seq_combined = ''
        for cds in part.location.parts:
            cds_count += 1
            cds_attributes = ['ID=' + 'cds_' + gene_id + '_' + str(cds_count),
                              'Parent=' + gene_id,
                              'product=30S ribosomal protein S12']
            features_list.append(_get_record(cds, 'CDS', cds_attributes))
            if cds.strand == -1:
                seq_combined += get_seq(genome, cds.start+1, cds.end).reverse_complement()
            else:
                seq_combined += get_seq(genome, cds.start+1, cds.end)
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


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description='Check your chloroplast genome gff file')

    parser.add_argument('-i', '--info_table', required=True,
                        help='<file_path>  meta information table which has five columns: Geseq gff path, seq path, '
                             'seqid, locus prefix, renumber result path. If you do not need renumber, just provide gff '
                             'and seq')

    parser.add_argument('-c', '--cds', action='store_true', default=False,
                        help='check cds legitimacy')

    parser.add_argument('-s', '--hs', action='store_true', default=False,
                        help='check house-keeping gene (matK, rbcL)')

    parser.add_argument('-n', '--name', action='store_true', default=False,
                        help='check whether gene names are legal name')

    parser.add_argument('-r', '--region', action='store_true', default=False,
                        help='check whether gene region duplicated')

    parser.add_argument('-e', '--renumber', action='store_true', default=False,
                        help='renumber gene locus suffix')

    args = parser.parse_args()

    info_table = pd.read_table(args.info_table, names=['raw_gff_path', 'seq_path', 'seq_id', 'prefix', 'result'])
    for ind, row in info_table.iterrows():
        print(os.path.basename(row['raw_gff_path']))
        tmp_check = CheckCp(row['raw_gff_path'])
        if args.cds:
            tmp_check.check_cds(row['seq_path'])
        if args.renumber:
            result_gff = tmp_check.renumber(row['seq_id'], row['prefix'])
            result_gff.to_csv(row['result'], sep='\t', index=False)
        if args.hs:
            tmp_check.check_hs()
        if args.name:
            tmp_check.check_name()
        if args.region:
            tmp_check.check_region()
