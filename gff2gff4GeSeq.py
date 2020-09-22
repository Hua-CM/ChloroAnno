# -*- coding: utf-8 -*-
# @Time : 2020/7/7 14:26
# @Author : Zhongyi Hua
# @FileName: gff2gff4GeSeq.py
# @Usage: Parse a GeSeq gff result to a submission gff file (e.g. remove useless note, add parent/child relationship)
# @Note: Please use checkgff.py in HuaSmallTools to check biological error after change.
# @E-mail: njbxhzy@hotmail.com

import pandas as pd
import gffutils
import os
from check import CheckCp

location = os.path.abspath(os.path.join(__file__, '..'))
tb_cds = pd.read_table(os.path.join(location, 'ref/cds.txt'))
tb_rna = pd.read_table(os.path.join(location, 'ref/rna.txt'))
product_dict = {record['name']: record['product'] for record in tb_cds.to_dict('records')}
rna_dict = {record['name']: record['product'] for record in tb_rna.to_dict('records')}
replace_dict = {'petE':	 'petG',
                'psbG':	 'ndhK',
                'lhbA':	 'psbZ'}


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


def change_gff(raw_gff_path, new_gff_path, seqid, species_pre):
    print('Start change', os.path.basename(new_gff_path))
    gff_file = gffutils.create_db(raw_gff_path, ':memory:', merge_strategy='create_unique')
    feature_list = []
    gene_count = 0
    for gene in gff_file.features_of_type('gene', order_by='start'):
        gene_count += 1
        gene_id = species_pre + '%03d' % gene_count
        gene_name = gene.attributes['gene'][0]
        if (gene_name == 'rps12') or ('fragment' in gene_name):
            continue
        elif gene_name in replace_dict.keys():
            gene_name = replace_dict[gene_name]
        gene_type = gene.attributes['gene_biotype'][0]
        gene_attributes = ['ID=' + gene_id,
                           'Name=' + gene_name,
                           'gene_biotype=' + gene_type
                           ]
        feature_list.append(get_record(gene, 'gene', gene_attributes))
        # parse child feature
        child_count = 0
        # protein coding gene
        if gene_type == 'protein_coding':
            next_phase = 0
            for exon in gff_file.children(gene,
                                          featuretype='exon',
                                          order_by='start',
                                          reverse=False if gene.strand == '+' else True):
                child_count += 1
                cds_attributes = ['ID=' + 'cds_' + gene_id + '_' + str(child_count),
                                  'Parent=' + gene_id,
                                  'product=' + product_dict.get(gene_name, 'hypothetical protein')]
                cds_record = get_record(exon, 'CDS', cds_attributes)
                cds_record.update({'phase': next_phase})
                next_phase = (3 - ((exon.end - exon.start+1-next_phase) % 3)) % 3
                feature_list.append(cds_record)
            if child_count == 0:
                cds_attributes = ['ID=' + 'cds_' + gene_id + '_1',
                                  'Parent=' + gene_id,
                                  'product=' + product_dict.get(gene_name, 'hypothetical protein')]
                # change phase
                cds_record = get_record(gene, 'CDS', cds_attributes)
                cds_record.update({'phase': next_phase})
                feature_list.append(cds_record)
        # RNA
        else:
            rna_attributes = ['ID=' + 'rna_' + gene_id,
                              'Parent=' + gene_id,
                              'product=' + rna_dict.get(gene_name, 'hypothetical protein')]
            feature_list.append(get_record(gene, gene_type, rna_attributes))
            for exon in gff_file.children(gene,
                                          featuretype='exon',
                                          order_by='start'):
                child_count += 1
                rna_attributes = ['ID=' + 'exon_' + gene_id + '_' + str(child_count),
                                  'Parent=' 'rna_' + gene_id
                                  ]
                feature_list.append(get_record(exon, 'exon', rna_attributes))
    result_gff = pd.DataFrame.from_dict({index: record for index, record in enumerate(feature_list)}, 'index')
    result_gff['seqid'] = seqid
    result_gff['score'] = '.'
    result_gff['source'] = 'GeSeq'
    result_gff = result_gff[["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]]
    result_gff.to_csv(new_gff_path, sep='\t', index=False, header=False)
    print('change done')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description='Change GeSeq result gff to a version that meet submission requirement')
    parser.add_argument('-i', '--info_table', required=True,
                        help='<file_path>  information table which has four columns: Geseq gff path, '
                             'result path, seqid, locus prefix')
    args = parser.parse_args()
    info_table = pd.read_table(args.info_table, names=['raw_gff_path', 'new_gff_path', 'seq_id', 'species_id'])
    for ind, row in info_table.iterrows():
        change_gff(*row.to_list())
        tmp_check = CheckCp(row['new_gff_path'])
        tmp_check.check_name()
        tmp_check.check_region()

