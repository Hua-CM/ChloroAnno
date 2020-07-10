# -*- coding: utf-8 -*-
# @Time : 2020/6/18 21:49
# @Author : Zhongyi Hua
# @FileName: gff2genbank.py
# @Usage: change gff to Genbank
# @Note:  only for our gff
# @E-mail: njbxhzy@hotmail.com

import gffutils
import pandas as pd
from Bio import SeqIO
from Bio import SeqRecord
from Bio import SeqFeature
from Bio.SeqFeature import ExactPosition, CompoundLocation, FeatureLocation
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio.Seq import Seq
from collections import OrderedDict

tb_map = pd.read_table(r'D:\ToxicDatabase\序列提交\NCBI\cds.txt', names=['gene_name', 'product'])
tb_trna = pd.read_table(r'D:\ToxicDatabase\序列提交\NCBI\tRNA.txt')
product_dict = {record['gene_name']: record['product'] for record in tb_map.to_dict('records')}
trna_dict = {record['trna']: record['product'] for record in tb_trna.to_dict('records')}
getstrand = lambda x: 1 if x.strand == '+' else -1


def gff2genbank(gff_path, seq_path, organism):
    gff_file = gffutils.create_db(gff_path, ':memory:', merge_strategy='create_unique')
    asm_seq = SeqIO.read(seq_path, 'fasta')
    seqfeature_list = list()
    # source feature
    seqfeature_list.append(
        SeqFeature.SeqFeature(FeatureLocation(ExactPosition(0), ExactPosition(len(asm_seq))),
                              strand=1,
                              qualifiers=OrderedDict([('organism', organism),
                                                      ('organelle', 'plastid:chloroplast'),
                                                      ('mol_type', 'genomic DNA')
                                                      ]),
                              type='source'
                              )
                            )
    # gene
    for gene in gff_file.features_of_type('gene', order_by='start'):
        seqfeature_list.append(
            SeqFeature.SeqFeature(FeatureLocation(ExactPosition(gene.start-1),
                                                  ExactPosition(gene.end),
                                                  strand=getstrand(gene)),
                                  strand=getstrand(gene),
                                  type='gene',
                                  qualifiers=OrderedDict([('gene', gene.attributes['Name'][0]),
                                                          ('locus_tag', gene.id)
                                                          ])
                                  )
                                )
        # CDS
        cds_list = []
        for cds in gff_file.children(gene, featuretype='CDS', order_by='start'):
            cds_list.append(FeatureLocation(ExactPosition(cds.start - 1),
                                            ExactPosition(cds.end),
                                            strand=getstrand(cds)))
        if len(cds_list) > 1:
            if cds_list[0].strand == -1:
                cds_list.reverse()
            location = CompoundLocation(cds_list)
        elif len(cds_list) == 1:
            location = cds_list[0]
        else:
            location = None
        if location is not None:
            seqfeature_list.append(
                SeqFeature.SeqFeature(location,
                                      strand=getstrand(gene),
                                      type='CDS',
                                      qualifiers=OrderedDict([('gene', gene.attributes['Name'][0]),
                                                              ('codon_start', 1),
                                                              ('transl_table', 11),
                                                              ('product', product_dict[gene.attributes['Name'][0]]),
                                                              ('translation', location.extract(asm_seq.seq).translate(table=11).__str__().strip('*'))
                                                              ])
                                      )
            )
        # exon
        for exon in gff_file.children(gene, featuretype='exon', order_by='start'):
            seqfeature_list.append(
                SeqFeature.SeqFeature(FeatureLocation(ExactPosition(exon.start - 1),
                                                      ExactPosition(exon.end),
                                                      strand=getstrand(exon)),
                                      strand=getstrand(exon),
                                      type='exon',
                                      qualifiers=OrderedDict([('gene', gene.attributes['Name'][0]),
                                                              ])
                                      )
            )
        # rRNA
        for rrna in gff_file.children(gene, featuretype='rRNA', order_by='start'):
            seqfeature_list.append(
                SeqFeature.SeqFeature(FeatureLocation(ExactPosition(rrna.start - 1),
                                                      ExactPosition(rrna.end),
                                                      strand=getstrand(rrna)),
                                      strand=getstrand(rrna),
                                      type='rRNA',
                                      qualifiers=OrderedDict([('gene', gene.attributes['Name'][0]),
                                                              ('product', product_dict[gene.attributes['Name'][0]])
                                                              ])
                                      )
            )
        # tRNA
        trna_list = []
        for trna in gff_file.children(gene, featuretype='tRNA', order_by='start'):
            trna_list.append(FeatureLocation(ExactPosition(trna.start - 1),
                                             ExactPosition(trna.end),
                                             strand=getstrand(trna)))
        if len(trna_list) > 1:
            if trna_list[0].strand == -1:
                trna_list.reverse()
            location = CompoundLocation(trna_list)
        elif len(trna_list) == 1:
            location = trna_list[0]
        else:
            location = None
        if location is not None:
            seqfeature_list.append(
                SeqFeature.SeqFeature(location,
                                      strand=getstrand(gene),
                                      type='tRNA',
                                      qualifiers=OrderedDict([('gene', gene.attributes['Name'][0]),
                                                              ('product',  trna_dict[gene.attributes['Name'][0]])
                                                              ])
                                      )
            )
    # annotation
    asm_annotation = {'molecule_type': 'DNA',
                      'topology': 'circular',
                      'date': '11-JUN-2020',
                      'organism': organism}
    return SeqRecord.SeqRecord(id=asm_seq.id,
                               seq=Seq(str(asm_seq.seq), IUPACAmbiguousDNA()),
                               features=seqfeature_list,
                               annotations=asm_annotation)

