# -*- coding: utf-8 -*-
# @Time    : 2021/9/2 21:12
# @Author  : Zhongyi Hua
# @File    : gbk2gff.py
# @Usage   : change Genbank to gff
# @Note    : mainly for NCBI gff
# @E-mail  : njbxhzy@hotmail.com

from collections import defaultdict
from pathlib import Path
import argparse


from BCBio import GFF
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, CompoundLocation

from ABC import MySeqFeature

FEATURE_ORDER = {'gene': 0, 'CDS': 1, 'tRNA': 2, 'rRNA': 3}
BIOTYPE = {'CDS': 'protein_coding', 'rRNA': 'rRNA', 'tRNA': 'tRNA'}


def parse_rps12(_cds: SeqFeature):
    def judge_length(gene_length: int):
        _dis = {'114': abs(gene_length - 114),
                '232': abs(gene_length - 232),
                '26': abs(gene_length - 26)}
        return min(_dis, key=_dis.get)

    _gene1 = SeqFeature(type='gene')
    _gene2 = SeqFeature(type='gene')
    _gene2.location = []
    _gene2.sub_features = []
    gene_qualifiers = defaultdict(list)
    # gene common
    for _key in ['gene', 'locus_tag']:
        gene_qualifiers[_key] = _cds.qualifiers.get(_key)
    gene_qualifiers['exception'] = ['trans-splicing']
    gene_qualifiers['gene_biotype'] = ['protein_coding']
    gene_qualifiers['gbkey'] = ['gene']
    gene_qualifiers['ID'] = gene_qualifiers['locus_tag']
    # cds common
    cds_qualifiers = defaultdict(list)
    for _key in ['gene', 'locus_tag', 'product', 'protein_id', 'transl_table', 'db_xref']:
        cds_qualifiers[_key] = _cds.qualifiers.get(_key)
    cds_qualifiers['gbkey'] = ['CDS']
    cds_qualifiers['Parent'] = cds_qualifiers['locus_tag']
    cds_count = 1
    for location in _cds.location.parts:
        if judge_length(location.__len__()) == '114':
            _gene1.location = location
            _gene1.strand = location.strand
            _gene1.qualifiers = gene_qualifiers.copy()
            _gene1.qualifiers['part'] = ['1/2']
            _cds1 = SeqFeature(location=location, type='CDS', strand=location.strand, qualifiers=cds_qualifiers.copy())
            _cds1.qualifiers['ID'] = ['cds-' + cds_qualifiers['locus_tag'][0] + '-1']
            _gene1.sub_features = [_cds1]
        else:
            cds_count += 1
            _gene2.qualifiers = gene_qualifiers.copy()
            _gene2.qualifiers['part'] = ['2/2']
            _gene2.location.append(location)
            _cds23 = SeqFeature(location=location, type='CDS', strand=location.strand, qualifiers=cds_qualifiers.copy())
            _cds23.qualifiers['ID'] = ['cds-' + cds_qualifiers['locus_tag'][0] + '-%d' % cds_count]
            _gene2.sub_features += [_cds23]
    _gene2.location = CompoundLocation(_gene2.location)
    _gene2.strand = _gene2.location.strand
    return [_gene1, _gene2]


def gbk2gff(seq_in):
    """Parse GenBank with locus tag

    Args:
        seq_in (SeqRecord): The import GenBank file in SeqRecord class

    Returns:
        top_feature: The BCBio.GFF readable file 
    """
    features = [_ for _ in seq_in.features if _.type in ['gene', 'CDS', 'tRNA', 'rRNA']]
    features.sort(key=lambda x: (x.location.start, FEATURE_ORDER[x.type]))
    out_features = defaultdict()
    rps12_features = []
    for feature in features:
        _feature = MySeqFeature()
        _feature.inherit(feature)
        if _feature.qualifiers.get('gene') == ['rps12']:
            if _feature.type == 'CDS':
                rps12_features.append(_feature)
            continue
        if _feature.type == 'gene':
            _feature.qualifiers['ID'] = _feature.qualifiers['locus_tag']
            _feature.qualifiers['gbkey'] = ['gene']
            out_features.setdefault(_feature.qualifiers['locus_tag'][0], _feature).update_subfeature()
        elif _feature.type == 'CDS':
            _count = 0
            for location in _feature.location.parts:
                _count += 1
                new_qualifiers = defaultdict(list)
                for _key in ['gene', 'locus_tag', 'product', 'protein_id', 'transl_table', 'db_xref']:
                    new_qualifiers.setdefault(_key, _feature.qualifiers.get(_key))
                new_qualifiers['ID'] = ['cds-' + new_qualifiers['locus_tag'][0] + '-%d' % _count]
                new_qualifiers['Parent'] = new_qualifiers['locus_tag']
                new_qualifiers['gbkey'] = ['CDS']
                out_features.setdefault(_feature.qualifiers['locus_tag'][0], _feature).update_subfeature(
                    SeqFeature(location, 'CDS', strand=location.strand, id=new_qualifiers['ID'], qualifiers=new_qualifiers)
                )
        elif _feature.type in ['rRNA', 'tRNA']:
            _count = 0
            _feature.qualifiers['ID'] = ['rna-' + _feature.qualifiers['locus_tag'][0]]
            _feature.qualifiers['Parent'] = _feature.qualifiers['locus_tag']
            _feature.qualifiers['gbkey'] = [_feature.type]
            for location in _feature.location.parts:
                _count += 1
                new_qualifiers = defaultdict(list)
                for _key in ['gene', 'locus_tag', 'product', 'db_xref']:
                    new_qualifiers.setdefault(_key, _feature.qualifiers.get(_key))
                new_qualifiers['ID'] = ['exon-' + new_qualifiers['locus_tag'][0] + '-%d' % _count]
                new_qualifiers['Parent'] = _feature.qualifiers['ID']
                new_qualifiers['gbkey'] = ['exon']
                _feature.update_subfeature(
                    SeqFeature(location, 'exon', strand=location.strand, id=new_qualifiers['ID'], qualifiers=new_qualifiers)
                )
            out_features.setdefault(_feature.qualifiers['locus_tag'][0], _feature).update_subfeature(_feature)
        else:
            continue
    features = list(out_features.values())
    for rps12 in rps12_features:
        features += parse_rps12(rps12)
    # add gene_biotype in qualifiers
    for feature in features:
        if feature.type == 'gene':
            try:
                feature.qualifiers['gene_biotype'] = [BIOTYPE[feature.sub_features[0].type]]
            except:
                continue
    # add CDS phase for protein_coding gene
    for feature in features:
        if feature.qualifiers.get('gene_biotype') == ['protein_coding']:
            set_phase(feature)
    top_feature = SeqRecord(id=seq_in.id, seq=seq_in.seq, features=features)
    return top_feature


def parse_args():
    """Parse arguments
    """
    parser = argparse.ArgumentParser(
        description='Change gff to genbank format')
    parser.add_argument('-i', '--info_table', required=True,
                        help='<file_path>  information table (tab-separate) which has two columns:genbank path, '
                             'gff path')
    _args = parser.parse_args()
    return _args


def main(args):
    """main interface
    """
    with open(args.info_table) as f_in:
        for _line in f_in.read().strip().split('\n'):
            seq_in = SeqIO.read(Path(_line.split('\t')[0]), 'fasta')
            genome = gbk2gff(seq_in)
            GFF.write([genome], open(_line.split('\t')[1], 'w'))


if __name__ == '__main__':
    main(parse_args())
