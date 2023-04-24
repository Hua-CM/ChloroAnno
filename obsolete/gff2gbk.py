# -*- coding: utf-8 -*-
# @Time : 2020/6/18 21:49
# @Author : Zhongyi Hua
# @FileName: gff2genbank.py
# @Usage: change gff to Genbank
# @Note:  only for gffs in GWH!!!
# @E-mail: njbxhzy@hotmail.com

from collections import defaultdict
import argparse
import time

from BCBio import GFF
from Bio import SeqIO
from Bio import SeqRecord
from Bio.SeqFeature import ExactPosition, FeatureLocation

from ABC import MySeqFeature


def gff2genbank(gff_path, seq_path, organism):
    gff_file = [_ for _ in GFF.parse(gff_path)]
    asm_seq = SeqIO.read(seq_path, 'fasta')
    genbank_header = defaultdict(str, {'organism': organism,
                                       'organelle': 'plastid:chloroplast',
                                       'molecule_type': 'DNA',
                                       'topology': 'circular',
                                       'data_file_division': 'PLN',
                                       'date': time.strftime("%d-%b-%Y", time.localtime()).upper(),
                                       'source': 'CGIR'})

    seqfeature_dict = defaultdict(MySeqFeature)
    # source feature
    seqfeature_dict['source'] = MySeqFeature(FeatureLocation(ExactPosition(1), ExactPosition(len(asm_seq))),
                                             strand=1,
                                             qualifiers={'organism': organism,
                                                         'organelle': 'plastid:chloroplast',
                                                         'mol_type': "genomic DNA"},
                                             type='source')
    seqfeature_dict['source'].location = [seqfeature_dict['source'].location]
    # gene
    for feature in gff_file[0].features:
        gene = MySeqFeature()
        gene.inherit(feature)
        seqfeature_dict.setdefault(gene.qualifiers['ID'][0], gene).update_location(gene.location)
        # CDS
        for subfeature in gene.sub_features:
            child = MySeqFeature()
            child.inherit(subfeature)
            if gene.qualifiers['gene_biotype'][0] == 'protein_coding':
                _prefix = 'cds_'
            elif gene.qualifiers['gene_biotype'][0] == 'rRNA':
                _prefix = 'rrna_'
            elif gene.qualifiers['gene_biotype'][0] == 'tRNA':
                _prefix = 'trna_'
            else:
                _prefix = 'other_'
            seqfeature_dict.setdefault(_prefix + gene.qualifiers['ID'][0], child).update_location(child.location)
    # reset information like NCBI
    for _key, _feature in seqfeature_dict.items():
        # qualifiers
        _old = _feature.qualifiers
        if _feature.type == 'gene':
            new_qualitier = {'gene': _old.get('gene') if 'gene' in _old else _old.get('Name'),
                             'locus_tag': _old['ID'][0],
                             'db_xref': 'GeneID:' + _old['Accession'][0]}
        elif _feature.type in ['CDS', 'rRNA', 'tRNA']:
            new_qualitier = {'gene': seqfeature_dict[_old['Parent'][0]].qualifiers.get('gene'),
                             'locus_tag': _old['Parent'][0],
                             'product': _old.get('product')[0] if _old.get('product') else '',
                             'db_xref': seqfeature_dict[_old['Parent'][0]].qualifiers.get('db_xref'),
                             'transl_table': '11'}
            if _feature.type == 'CDS':
                new_qualitier['protein_id'] = _old.get('Protein_Accession')
        else:
            new_qualitier = _old
        _feature.qualifiers = new_qualitier
    # output
    seqfeature_list = list(seqfeature_dict.values())
    return SeqRecord.SeqRecord(id=asm_seq.id,
                               seq=asm_seq.seq,
                               features=seqfeature_list,
                               annotations=genbank_header)


def getArgs():
    parser = argparse.ArgumentParser(
        description='Change gff to genbank format')
    parser.add_argument('-i', '--info_table', required=True,
                        help='<file_path>  information table (tab-separate) which has four columns: gff path, '
                             'fasta path, organism name, output path')
    _args = parser.parse_args()
    return _args


def main(args):
    with open(args.info_table) as f_in:
        for _line in f_in.read().strip().split('\n'):
            try:
                genome = gff2genbank(*_line.split('\t')[:3])
                SeqIO.write(genome, _line.split('\t')[3], 'genbank')
            except:
                print(_line.split('\t')[1])


if __name__ == '__main__':
    main(getArgs())
