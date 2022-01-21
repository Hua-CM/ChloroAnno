# -*- coding: utf-8 -*-
# @Time : 2022/1/18 22:26
# @Author : Zhongyi Hua
# @File : curate_geseq.py.py
# @Usage :
# @Note :
# @E-mail : njbxhzy@hotmail.com

from BCBio import GFF
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from utilities.check import CheckCp2, product_dict, rna_dict


class MySeqFeature(SeqFeature):
    def inherit(self, _parent):
        self.__dict__ = _parent.__dict__

    def update_subfeature(self, _subfeature=None):
        if not hasattr(self, 'sub_features'):
            self.sub_features = []
        if _subfeature is not None:
            self.sub_features.append(_subfeature)


def tidy_gff(gff_path, _prefix):
    gff_file = next(GFF.parse(gff_path))
    gene_count = 0
    feature_lst = []
    gff_file.features.sort(key=lambda x: x.location.start)
    for feature in gff_file.features:
        if not feature.type == 'gene':
            continue
        gene_type = feature.qualifiers.get('gene_biotype')[0]
        gene_name = feature.qualifiers.get('gene')[0]
        gene_count += 1
        gene_id = _prefix + '%03d' % gene_count
        _gene = MySeqFeature()
        _gene.update_subfeature()
        gene_qualifiers = {'ID':  [gene_id],
                           'gene_biotype': ['tRNA' if gene_type == 'tRNA' else 'protein_coding'],
                           'Name': [gene_name]
                           }
        _gene.id, _gene.type, _gene.location, _gene.qualifiers = gene_id, 'gene', feature.location, gene_qualifiers
        if gene_type in ['tRNA', 'rRNA']:  # single tRNA and rRNA
            _rna = MySeqFeature()
            _rna.update_subfeature()
            rna_qualifiers = {'ID': ['rna_%s' % gene_id],
                              'Parent': [gene_id],
                              'product': [rna_dict.get(gene_name, 'hypothetical protein')]
                              }
            _rna.id, _rna.type, _rna.location, _rna.qualifiers = \
                'rna_%s' % gene_id, gene_type, feature.location, rna_qualifiers
            _gene.update_subfeature(_rna)
        # multiple exons
        if feature.sub_features:
            if gene_type == 'tRNA':
                parent_feature = _rna
                sub_prefix = 'exon'
            else:  # protein coding
                parent_feature = _gene
                sub_prefix = 'cds'
            sub_count = 0
            for subfeature in feature.sub_features:
                if subfeature.type == 'exon':
                    sub_count += 0
                    sub_feature = SeqFeature(id='%s_%s_%s' % (sub_prefix, gene_id, sub_count),
                                             type='exon' if gene_type == 'tRNA' else 'CDS',
                                             location=subfeature.location,
                                             qualifiers={'ID': ['exon_%s_%s' % (gene_id, sub_count)],
                                                         'Parent': [parent_feature.id]})
                    if gene_type == 'protein_coding':
                        sub_feature.qualifiers['product'] = [product_dict.get(gene_name, 'hypothetical protein')]
                    parent_feature.update_subfeature(sub_feature)
        elif gene_type == 'protein_coding':  # single CDS
            _cds = SeqFeature(id='cds_%s' % gene_id,
                              type='CDS',
                              location=feature.location,
                              qualifiers={'ID': ['cds_%s' % gene_id],
                                          'Parent': [gene_id],
                                          'product': [product_dict.get(gene_name, 'hypothetical protein')]})
            _gene.update_subfeature(_cds)
        feature_lst.append(_gene)
    return SeqRecord(id=gff_file.id, seq=gff_file.seq, features=feature_lst)


def parseArgs():
    import argparse
    parser = argparse.ArgumentParser(
        description='Change GeSeq result gff to a version that meet submission requirement')
    parser.add_argument('-i', '--info_table', required=True,
                        help='<file_path>  information table which has three columns: Geseq gff path, '
                             'result path, locus prefix')
    args = parser.parse_args()
    return args


def main(args):
    with open(args.info_table) as f_in:
        info_table = f_in.read().splitlines()
    for row in info_table:
        raw_gff_path, new_gff_path, species_id = row.split()
        curated_gff = tidy_gff(raw_gff_path, species_id)
        tmp_check = CheckCp2(curated_gff)
        tmp_check.check_name()
        tmp_check.check_region()
        with open(new_gff_path, 'w') as f_out:
            GFF.write([curated_gff], f_out, include_fasta=False)


if __name__ == '__main__':
    main(parseArgs())
