# -*- coding: utf-8 -*-
# @Time : 2020/10/10 15:05
# @Author : Zhongyi Hua
# @FileName: gff2tbl.py
# @Usage: 
# @Note:
# @E-mail: njbxhzy@hotmail.com

import pandas as pd
import gffutils


def gff2gbl(_gff_path, _tbl_path):
    gff_file = gffutils.create_db(_gff_path, ':memory:', merge_strategy='create_unique')
    _df_list = []
    seq_accession = None
    for _gene in gff_file.features_of_type('gene', order_by='start'):
        seq_accession = _gene.seqid
        if _gene.attributes['Name'][0] == 'rps12':
            continue
        _feature_list = []
        _feature_list.append({'c1': _gene.start if _gene.strand == '+' else _gene.end,
                              'c2': _gene.end if _gene.strand == '+' else _gene.end,
                              'c3': 'gene'})
        _feature_list.append({'c4': 'gene', 'c5': _gene.attributes['Name'][0]})
        _feature_list.append({'c4': 'locus_tag', 'c5': _gene.attributes['ID'][0]})
        if _gene.attributes.get('pseudo'):
            _feature_list.append({'c4': 'pseudo'})
            _df_list.append(pd.DataFrame(_feature_list))
            continue
        _df_list.append(pd.DataFrame(_feature_list))
        # child
        _feature_list = []
        _product = None
        if _gene.attributes['gene_biotype'][0] == 'protein_coding':
            for _cds in gff_file.children(_gene,
                                          order_by='start',
                                          reverse=False if _gene.strand == '+' else True):
                _feature_list.append({'c1': _cds.start if _gene.strand == '+' else _cds.end,
                                      'c2': _cds.end if _gene.strand == '+' else _cds.start})
                _product = _cds.attributes['product'][0]
            _feature_list[0].update({'c3': 'CDS'})
            _feature_list.append({'c4': 'codon_start', 'c5': '1'})
            _feature_list.append({'c4': 'product', 'c5': _product})
            _feature_list.append({'c4': 'transl_table', 'c5': '11'})
            _df_list.append(pd.DataFrame(_feature_list))
        else:
            for _rna in gff_file.children(_gene,
                                          order_by='start',
                                          reverse=False if _gene.strand == '+' else True):
                if _rna.featuretype in ['tRNA', 'rRNA']:
                    _product = _rna.attributes['product'][0]
                    continue
                else:
                    if _gene.strand == '+':
                        _feature_list.append({'c1': _rna.start, 'c2': _rna.end})
                    else:
                        _feature_list.append({'c1': _rna.end, 'c2': _rna.start})
            if len(_feature_list) == 0:
                _feature_list.append({'c1': _gene.start if _gene.strand == '+' else _gene.end,
                                      'c2': _gene.start if _gene.strand == '-' else _gene.end,
                                      'c3': _gene.attributes['gene_biotype'][0]})
                _feature_list.append({'c4': 'product', 'c5': _product})
                _feature_list.append({'c4': 'gene', 'c5': _gene.attributes['Name'][0]})
                _df_list.append(pd.DataFrame(_feature_list))
            else:
                _feature_list[0].update({'c3': _gene.attributes['gene_biotype'][0]})
                _feature_list.append({'c4': 'product', 'c5': _product})
                _feature_list.append({'c4': 'gene', 'c5': _gene.attributes['Name'][0]})
                _df_list.append(pd.DataFrame(_feature_list))
    _result_df = pd.concat(_df_list)
    _result_df.to_csv(_tbl_path, sep='\t', index=False, header=False)
    # take as file and character
    _correct_list = ['>Feature ' + seq_accession]
    with open(_tbl_path, 'r') as f:
        _tbl = f.readlines()
        for _line in _tbl:
            _correct_list.append(_line.rstrip().replace('.0', ''))
    with open(_tbl_path, 'w') as f:
        f.write('\n'.join(_correct_list))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description='Convert gff file to Feature table for submitting to NCBI')

    parser.add_argument('-i', '--input', required=True,
                        help='<file_path>  gff file path OR a two columns meta file(gff path / tbl path)')
    parser.add_argument('-o', '--output',
                        help='<file_path> Feature table(.tbl) path')
    parser.add_argument('-b', '--batch', action='store_true', default=False,
                        help='Use a meta file for Batch conversion')
    args = parser.parse_args()
    if args.batch:
        _meta = pd.read_table(args.input, names=['v1', 'v2'])
        for idx, row in _meta.iterrows():
            gff2gbl(row['v1'], row['v2'])
    else:
        gff2gbl(args.input, args.output)
