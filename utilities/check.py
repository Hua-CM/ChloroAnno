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

import portion as pt
import pandas as pd
from Bio import SeqIO
from BCBio import GFF
from Bio.SeqRecord import SeqRecord
import os, json


file_location = os.path.abspath(os.path.join(__file__, '../..'))
with open(os.path.join(file_location, 'ref/cds.json')) as f_in:
    product_dict = json.load(f_in)
with open(os.path.join(file_location, 'ref/rna.json')) as f_in:
    rna_dict = json.load(f_in)

standard_name_list = list(product_dict.keys())
standard_name_list += list(rna_dict.keys())


def check_cds(gene, genome_seq):
    """
    :param gene: SeqRecord with sub_features attribute
    :return: Bool
    """
    check_seq = ''
    gene.sub_features.sort(key=lambda x: x.location.start, reverse=True if gene.strand == -1 else False)
    for cds in gene.sub_features:
        check_seq += cds.location.extract(genome_seq)
    check_seq.translate(table=11, cds=True)


# Based on SeqRecord and SeqFearure
class CheckCp2:
    standard_name_list = standard_name_list

    def __init__(self, gff_path):
        self.gff = gff_path if isinstance(gff_path, SeqRecord) else next(GFF.parse(gff_path))

    def check_hs(self):
        """
        Check housekeeping gene: matK, rbcL
        :return: a message
        """
        gene_name_list = []
        for gene in self.gff.features:
            gene_name_list.append(gene.qualifiers.get('Name')[0])
        if ('matK' not in gene_name_list) and ('matk' not in gene_name_list):
            print('matK loss!')
        if ('rbcL' not in gene_name_list) and ('rbcl' not in gene_name_list):
            print('rbcL loss!')

    def check_name(self):
        for gene in self.gff.features:
            gene_name = gene.qualifiers.get('Name')[0]
            if not ((gene_name in CheckCp2.standard_name_list) or gene_name.startswith('orf')):
                print('check ' + gene_name)

    def check_region(self):
        """
        check duplicated gene region
        :return: messages
        """
        region_list = []
        locus_list = []
        for gene in self.gff.features:
            region_list.append(pt.closed(gene.location.start.__int__(), gene.location.end.__int__()))
            locus_list.append(gene.qualifiers.get('Name')[0])
        for i in range(len(region_list) - 1):
            if not region_list[i] < region_list[i + 1]:
                print(locus_list[i], region_list[i], ' and ', locus_list[i + 1], region_list[i + 1],
                      ' are duplicated')
        print('check duplicated region done')

    def renumber(self, prefix):
        print('Renumber gene id')
        gene_count = 0
        for gene in self.gff.features:
            gene.id = '%s%03d' % (prefix, gene_count)
            gene.qualifiers['ID'] = [gene.id]
            sub_count = 0
            gene_count += 1
            gene_type = gene.qualifiers['gene_biotype'][0]
            for subfeature in gene.sub_features:
                sub_count += 1
                subfeature.id = '%s_%s_%1d' % ('rna' if gene_type in ['tRNA', 'rRNA'] else 'cds', gene.id, sub_count)
                subfeature.qualifiers['ID'] = [subfeature.id]
                subfeature.qualifiers['Parent'] = [gene.id]
                if gene_type in ['tRNA', 'rRNA']:
                    exon_count = 0
                    for exon in subfeature.sub_features:
                        exon_count += 1
                        exon.id = 'exon_%s_%1d' % (gene.id, exon_count)
                        exon.qualifiers['ID'] = [exon.id]
                        exon.qualifiers['Parent'] = [subfeature.id]

    def check_cds(self, seq_ins, add_pseudo=False):
        self.gff.seq = seq_ins
        for gene in self.gff.features:
            gene_type = gene.qualifiers['gene_biotype'][0]
            if gene.qualifiers.get('Name')[0] == 'rps12':
                continue
            if gene_type == 'protein_coding':
                try:
                    check_cds(gene, self.gff.seq)
                except Exception as e:
                    if add_pseudo:
                        gene.qualifiers['pseudo'] = ['true']
                    print(gene.id)
                    print(e)


def add_rps12(pga_gb, new_gff, species_pre):
    def _get_record(location, fea_type, attributes):
        return {'type': fea_type,
                'start': location.start+1,
                'end': location.end,
                'strand': '-' if location.strand == -1 else '+',
                'phase': '.',
                'attributes': ";".join(attributes)}

    def get_seq(seq, start, end):
        return seq.seq[start - 1:end]

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
            'exception=trans-splicing;part=1/2']))
        features_list.append(_get_record(part.location, 'gene', part_attributes + [
            'exception=trans-splicing;part=2/2']))
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


def parseArgs():
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
    return args


def main(args):
    info_table = pd.read_table(args.info_table, names=['raw_gff_path', 'seq_path', 'prefix', 'result'])
    for ind, row in info_table.iterrows():
        print(os.path.basename(row['raw_gff_path']))
        tmp_check = CheckCp2(row['raw_gff_path'])
        if args.cds:
            tmp_check.check_cds(row['seq_path'])
        if args.renumber:
            tmp_check.renumber(row['prefix'])
            GFF.write(tmp_check.gff)
        if args.hs:
            tmp_check.check_hs()
        if args.name:
            tmp_check.check_name()
        if args.region:
            tmp_check.check_region()


if __name__ == '__main__':
    main(parseArgs())
