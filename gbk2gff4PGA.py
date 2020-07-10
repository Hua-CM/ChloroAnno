# -*- coding: utf-8 -*-
# @Time : 2020/1/1 17:12
# @Author : Zhongyi Hua
# @FileName: gbk2gff4PGA.py
# @Usage: Parse PGA GenBank file to a clean gff3 profile
# @Note:
# @E-mail: njbxhzy@hotmail.com

# load module
import pandas as pd
import os
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

tb_cds = pd.read_table('cds.txt')
tb_rna = pd.read_table('rna.txt')
product_dict = {record['name']: record['product'] for record in tb_cds.to_dict('records')}
rna_dict = {record['name']: record['product'] for record in tb_rna.to_dict('records')}


def get_phase(genbank_feature):
    """
    Get the phase of a record
    :param genbank_feature: a record feature (Class SeqFeature)
    :return: phase ("."/0/1/2)
    """
    try:
        phase = str(int(genbank_feature.qualifiers['codon_start'][0])-1)
    except KeyError:
        phase = '.'
    return phase


def fix_location(record_feature, plus_num):
    """
    Fix a subseq loaction problem. If you use position to subseq a child feature from a SeqRecord instance, the start
    location of the child feature would return to zero. So we should plus the real start position to it.
    :param record_feature: a record feature (Class SeqFeature)
    :param plus_num: the real start position
    :return: none
    """
    if len(record_feature.location.parts) == 1:
        record_feature.location = FeatureLocation(record_feature.location.start + plus_num,
                                                  record_feature.location.end + plus_num,
                                                  strand=record_feature.location.strand)
    else:
        location_list = []
        for _ in record_feature.location.parts:
            location_list.append(FeatureLocation(_.start + plus_num,
                                                 _.end + plus_num,
                                                 strand=_.strand)
                                 )
        record_feature.location = CompoundLocation(location_list)


def remove_none_location(genome_record):
    """
    Since chloroplast is a circle, some annotation may cross the end-begin gap, remove them
    :param genome_record: a Class Bio.SeqRecord.SeqRecord object
    :return: a Class Bio.SeqRecord.SeqRecord object
    """
    unwanted_list = []
    for index, x in enumerate(genome_record.features):
        if x.location is None:
            unwanted_list.append(index)
    for unwanted_index in sorted(unwanted_list, reverse=True):
        del genome_record.features[unwanted_index]


def get_record(record_feature, attributes):
    """
    Transform the SeqFeature Class to a line and add it to a gff3 DataFrame
    :param record_feature: a record feature (Class SeqFeature)
    :param attributes: a list of character that positioned in gff3 column9
    :return: a dict have six columns
    """
    feature_record = {'type': record_feature.type,
                      'start': record_feature.location.start+1,
                      'end': record_feature.location.end,
                      'strand': '+' if record_feature.strand == 1 else '-',
                      'phase': get_phase(record_feature) if record_feature.type == 'CDS' else '.',
                      'attributes': ';'.join(attributes)}
    return feature_record


def main(genbank_path, new_gff_path, seq_id, species_id):
    print('Start change', os.path.basename(new_gff_path))
    records_list = []
    genome = SeqIO.read(genbank_path, "genbank")
    remove_none_location(genome)
    genome.features.sort(key=lambda x: x.location.start)
    gene_count = 0
    IR_count = 0
    for ele in genome.features:
        if ele.type == 'gene':
            if ele.qualifiers['gene'][0] == 'rps12':
                continue
            gene_count += 1
            ele.id = species_id + '%03d' % gene_count
            for child_feature in genome[ele.location.start:ele.location.end].features:
                fix_location(child_feature, ele.location.start)
                if child_feature.type != 'gene' and \
                        child_feature.location.start == ele.location.start and \
                        child_feature.location.end == ele.location.end:
                    child_feature.type = 'mRNA' if child_feature.type == 'CDS' else child_feature.type
                    if child_feature.qualifiers['gene'][0] == ele.qualifiers['gene'][0]:
                        # This module for protein coding gene CDS region
                        if child_feature.type == 'mRNA':
                            gene_attributes = ['ID=' + ele.id,
                                               'Name=' + ele.qualifiers['gene'][0],
                                               'gene_biotype=protein_coding']
                            records_list.append(get_record(ele, gene_attributes))
                            cds_count = 0
                            for cds in reversed(child_feature.location.parts):
                                cds_count += 1
                                cds_feature = SeqFeature(cds,
                                                         type='CDS',
                                                         qualifiers={'codon_start': child_feature.qualifiers['codon_start'][0]})
                                cds_feature.id = 'cds_' + species_id + '%03d' % gene_count + '_' + '%d' % cds_count
                                cds_attributes = ['ID=' + cds_feature.id,
                                                  'Parent=' + ele.id,
                                                  'product=' + child_feature.qualifiers['product'][0]
                                                  ]
                                records_list.append(get_record(cds_feature, cds_attributes))
                        # This module for rRNA and tRNA exon
                        else:
                            # gene
                            gene_attributes = ['ID=' + ele.id,
                                               'Name=' + ele.qualifiers['gene'][0],
                                               'gene_biotype=' + child_feature.type]
                            records_list.append(get_record(ele, gene_attributes))
                            # rna
                            child_feature.id = 'rna_' + species_id + '%03d' % gene_count
                            child_attributes = ['ID=' + child_feature.id,
                                                'Name=' + child_feature.qualifiers['gene'][0],
                                                'Parent=' + ele.id,
                                                'product=' + child_feature.qualifiers['product'][0]
                                                ]
                            records_list.append(get_record(child_feature, child_attributes))
                            exon_count = 0
                            # exon
                            for exon in reversed(child_feature.location.parts):
                                exon_count += 1
                                exon_feature = SeqFeature(exon, type='exon')
                                exon_feature.id = 'exon_' + species_id + '%03d' % gene_count + '_' + '%d' % exon_count
                                exon_attributes = ['ID=' + exon_feature.id, 'Parent=' + child_feature.id]
                                records_list.append(get_record(exon_feature, exon_attributes))
        elif ele.type == 'repeat_region':
            IR_count += 1
            gene_attributes = ['ID=IR' + str(IR_count), 'note=Inverted repeats']
            records_list.append(get_record(ele, gene_attributes))
    records_dict = {index: record for index, record in enumerate(records_list)}
    result_gff = pd.DataFrame.from_dict(records_dict, 'index')
    result_gff['seqid'] = seq_id
    result_gff['score'] = '.'
    result_gff['source'] = 'PGA'
    result_gff = result_gff[["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]]
    result_gff.to_csv(new_gff_path, sep='\t', header=False, index=False, encoding='utf8')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='This is the script for change PGA annotation result GenBank file  \
                                                  to its gff3 version')
    parser.add_argument('-i', '--info_table', required=True,
                        help='<file_path>  The sample information table')
    args = parser.parse_args()
    info_table = pd.read_table(args.info_table, names=['genbank_path', 'new_gff_path', 'seq_id', 'species_id'])
    for ind, row in info_table.iterrows():
        main(*row.to_list())
