# -*- coding: utf-8 -*-
# @Time : 2023/4/16 20:54
# @Author : Zhongyi Hua
# @FileName: input.py
# @Usage: Convert result from various databases and softwares to List[CoreSeqFeature]
# @Note:
# @E-mail: njbxhzy@hotmail.com
from pathlib import Path
from datetime import date
import warnings

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from BCBio import GFF
from .ABC import CoreSeqFeature, fix_location, feature2tbl
from .functions import Check


def read_geseq(file_path: Path):
    """Since geseq gff file is illegal, and only directly handled by
    Tidy function, a function was written for it directly.

    1. Remove "Note" and "gbkey" in GeSeq gff file qualifiers.
    2. Add CDS and tRNA to correct gene parent

    Args:
        file_path (Path): The GeSeq raw gff path

    Returns:
        top_feature (CoreSeqFeature): The SeqRecord contains CoreSeqFeature
    """
    gff_obj = next(GFF.parse(file_path))
    gene_feature_lst = []
    gene_feature = None
    rna_feature = None
    for feature in gff_obj.features:
        # Since the CoreSeqFeature initialize the RNA feature when CoreSeqFeature
        # update_subfeature the RNA subfeature, I use this method to update the RNA
        # subfeature before the next gene
        if feature.type == 'gene':
            if rna_feature:
                gene_feature.update_subfeature(rna_feature)
                rna_feature = None
            feature.qualifiers.pop('Note')
            feature.qualifiers.pop('gbkey')
            gene_feature = CoreSeqFeature(feature)
            gene_feature_lst.append(gene_feature)
            f_type = gene_feature.qualifiers['gene_biotype'][0]
            if f_type in ['tRNA', 'rRNA']:
                # To avoid pointer problem, use SeqFeature instead of MySeqFeature
                gene_feature.sub_features = []
                rna_feature = SeqFeature(location=gene_feature.location,
                                         type=f_type,
                                         qualifiers=gene_feature.qualifiers)
                rna_feature.sub_features = []
        elif feature.type  == 'CDS':
            gene_feature.update_subfeature(feature)
        elif feature.type in ['tRNA', 'rRNA']:
            # The tRNA and rRNA in GeSeq gff is the exon feature actually
            feature.type = 'exon'
            rna_feature.sub_features.append(feature)
    # Remove rps12
    gene_feature_lst = [_ for _ in gene_feature_lst if _.qualifiers['gene'] != ['rps12'] ]
    top_record = SeqRecord(id=gff_obj.id,
                           seq=gff_obj.seq,
                           features=gene_feature_lst)
    check_ins = Check(top_record)
    check_ins.renumber('GeSeq')
    top_record = check_ins.record
    return top_record

# The input IO
def read_anno(file_path: Path, file_type: str):
    """Read Chloroplast annotation from file and convert it 
    to the CoreSeqFeature

    Args:
        file_path (Path): A chloroplast gennome annotation file
        file_type (str): The file file format
    
    Return:
        top_feature (CoreSeqFeature): The SeqRecord contains CoreSeqFeature
    """
    if file_type == 'gb':
        # The GenBank header may raise an warning
        warnings.filterwarnings("ignore")
        return _gbk2core(SeqIO.read(file_path, 'gb'))
    if file_type == 'gff':
        return _gff2core(next(GFF.parse(file_path)))


def _gbk2core(genbank: SeqRecord):
    """Parse GenBank to CoreSeqFeature

    Args:
        genbank (SeqRecord): The import GenBank file in SeqRecord class
    
    **NOTICE:** The `SeqRecord.annotations` for GenBank file is empty, and thus, I keep it.

    Returns:
        top_record (List[CoreSeqFeature]): The SeqRecord contains CoreSeqFeature
    """
    FEATURE_ORDER = {'gene': 0, 'CDS': 1, 'tRNA': 2, 'rRNA': 3}
    features = [_ for _ in genbank.features if _.type in ['gene', 'CDS', 'tRNA', 'rRNA']]
    features.sort(key=lambda x: (x.location.start, FEATURE_ORDER[x.type]))
    gene_feature_lst = []
    for feature in features:
        if feature.type == 'gene':
            if feature.qualifiers.get('gene') == ['rps12'] or feature.qualifiers.get('Name') == ['rps12']:
                continue
            _feature = CoreSeqFeature(feature)
            for child_feature in genbank[_feature.location.start: _feature.location.end].features:
                if child_feature.type != 'gene':
                    # Avoid treating CDS inserted bewteen two exons from a tRNA as a subfeature
                    if _feature.qualifiers.get('gene_biotype') == ['tRNA'] and child_feature.type == 'CDS':
                        continue
                    fix_location(child_feature, _feature.location.start)
                    _feature.update_subfeature(child_feature)
            # Check children
            ## Avoid errors (a lot of CDS for a gene) caused by long distance between two exons.
            ## No gene on chloroplast has more than two exons, except for clpP and rps12.
            ## Therefore, for gene with more than two exons, only keep the first one and last one
            if len(_feature.sub_features) > 2 and _feature.qualifiers.get('gene') != ['clpP']:
                _tmp_lst = sorted(_feature.sub_features, key=lambda x: x.location.start)
                _tmp_lst = [_tmp_lst [0], _tmp_lst[-1]]
                _feature.sub_features = _tmp_lst
                del _tmp_lst
            gene_feature_lst.append(_feature)
    # Check duplicated features and remove
    position_set = set()
    gene_lst = []
    for gene in gene_feature_lst:
        if (gene.location.start, gene.location.end) in position_set:
            continue
        position_set.add((gene.location.start, gene.location.end))
        gene_lst.append(gene)
    # For non-Gene features , such as IR,Source.
    # Some GenBank files also contain mRNA and exon, remember to remove .
    gene_lst +=  [_ for _ in genbank.features if _.type not in ['gene', 'CDS', 'tRNA', 'rRNA', 'exon', 'intron','mRNA']]
    top_record = SeqRecord(id=genbank.id,
                           seq=genbank.seq,
                           features=gene_lst,
                           annotations=genbank.annotations)
    return top_record


def _gff2core(gff: SeqRecord):
    """Parse GFF to CoreSeqFeature

    Args:
        gff (SeqRecord): The import GFF file in SeqRecord class

    Returns:
        top_record (List[CoreSeqFeature]): The SeqRecord contains CoreSeqFeature
    """
    gene_feature_lst = []
    for feature in gff.features:
        if feature.type == 'gene':
            if feature.qualifiers.get('gene') == ['rps12'] or feature.qualifiers.get('Name') == ['rps12']:
                continue
            _feature = CoreSeqFeature(feature)
            _feature.sub_features = []
            for child_feature in feature.sub_features:
                if child_feature.type in ['tRNA', 'rRNA', 'CDS']:
                    # The tRNA and rRNA always have "exon" subfeatures, CoreSeqFeature would handle it when output
                    _feature.update_subfeature(child_feature)
                elif child_feature.type == 'mRNA':
                    # Some software with mRNA in their results. But even these softwares,
                    # the CDSs still belong to gene directly.
                    for child2_feature in child_feature.sub_features:
                        if child2_feature.type == 'CDS':
                            _feature.update_subfeature(child2_feature)
            gene_feature_lst.append(_feature)
    gene_feature_lst +=  [_ for _ in gff.features if _.type not in ['gene', 'CDS', 'tRNA', 'rRNA', 'exon', 'mRNA']]
    top_record = SeqRecord(id=gff.id, seq=gff.seq, features=gene_feature_lst)
    return top_record

# The output IO
def write_anno(top_record, file_path: Path, file_type: str):
    """Read Chloroplast annotation from file and convert it
     to the CoreSeqFeature

    Args:
        file_path (Path): _description_
        file_type (str): _description_
    
    Return:
        top_feature (CoreSeqFeature): The SeqRecord contains CoreSeqFeature
    """
    if file_type == 'gb':
        out_record = _core2gbk(top_record)
        SeqIO.write(out_record, file_path, 'gb')
    elif file_type == 'gff':
        out_record = _core2gff(top_record)
        GFF.write([out_record], open(file_path, 'w'), include_fasta=False)
    elif file_type == 'tbl':
        # A sting! Just write!
        out_record = _core2tbl(top_record)
        file_path.write_text(out_record)

def _keep_keys(feature, key_lst):
    new_qualifiers = {}
    for _key in key_lst:
        if _value := feature.qualifiers.get(_key):
            new_qualifiers[_key] = _value
    feature.qualifiers = new_qualifiers


def _core2gbk(top_record: SeqRecord):
    """_summary_
    The qualifiers should be output in GenBank file:
        gene: gene, locus_tag
        tRNA: gene, locus_tag, product
        rRNA: gene, locus_tag, product
        CDS:  gene, locus_tag, codon_start, transl_table, product
    
    **NOTICE:**
    1. The seq in SeqRecord should not be null for output GenBank file.
    2. The CDSs in CoreSeqFeature.sub_features are seperate and need to combine.
    3. The feature list for GenBank output is not *nested*. If nested, the SeqIO.write only
     write the top level features.
    4. The `SeqRecord.annotations` for GenBank File should not be empty.
      
    Args:
        top_record (SeqRecord): The SeqRecord with CoreSeqFeature features.
    
    Return:
        output_record (SeqRecord): The SeqRecord with unnested (Core)SeqFeature features.
    """
    feature_lst = []
    keep_gb_keys = ['locus_tag', 'gene', 'codon_start', 'product', 'transl_table', 'pseudo', 'trans_splicing']
    for gene in top_record.features:
        # non-gene features
        feature_lst.append(gene)
        if gene.type == 'gene':
            if 'pseudo' in gene.qualifiers.keys():
                gene.qualifiers['pseudo'] = [None]
                _keep_keys(gene, keep_gb_keys)
            # In some extremely rare cases, there are pseudogenes who don't possess CDS subfeatures
            if gene.qualifiers.get('gene_biotype') == ['protein_coding'] and len(gene.sub_features) > 0:
                _tmp_cds_feature = gene.gen_gb_cds()
                _keep_keys(_tmp_cds_feature, keep_gb_keys)
                feature_lst.append(_tmp_cds_feature)
            else:
                # For RNAs
                for _subfeature in gene.sub_features:
                    _keep_keys(_subfeature, keep_gb_keys)
                feature_lst +=  gene.sub_features
        # In some extremely rare cases, there are pseudogenes whose type is pseudo
        if gene.type == 'pseudogene':
            gene.qualifiers['pseudo'] = [None]
    # Only for GenBank File
    if top_record.annotations:
        annotations = top_record.annotations
    else:
        annotations = {
            'molecule_type': 'DNA',
            'topology': 'circular',
            'data_file_division': 'PLN',
            'date': date.today().strftime('%d-%b-%Y'),
            'accessions': [top_record.id],
            'sequence_version': 1,
            'keywords': [''],
        }
    output_record = SeqRecord(id=top_record.id,
                              seq=top_record.seq,
                              features=feature_lst,
                              annotations=annotations)
    return output_record


def _core2gff(top_record):
    """Technically, the SeqRecord with a list of CoreSeqFeature genes could
     be output directly. However, (1) some values in qualifers need to remove.
      especially 'translation', and 'codon_start' in CDS subfeatures;(3) exon
      sub_features under protein_coding gene need to be removed.
    
    **NOTICE:** This funciton directly modify the input object. Beacuse after this 
    function, object will be written to file.

    Args:
        top_record (_type_): _description_
    """
    keep_cds_keys = ['ID', 'Parent', 'Name', 'gene', 'phase', 'product', 'transl_table']
    for gene in top_record.features:
        # non-gene features
        if gene.type == 'gene' and gene.qualifiers.get('gene_biotype') == ['protein_coding']:
            exclude_exon_sub = []
            for sub_feature in gene.sub_features:
                if sub_feature.type == 'CDS':
                    new_qualifiers = {}
                    for _key in keep_cds_keys:
                        _value = sub_feature.qualifiers.get(_key)
                        if _value:
                            new_qualifiers[_key] = _value
                    sub_feature.qualifiers = new_qualifiers
                    exclude_exon_sub.append(sub_feature)
            gene.sub_features = exclude_exon_sub
        if 'pseudo' in gene.qualifiers:
            gene.qualifiers['pseudo'] = ['true']
        if 'trans_splicing' in gene.qualifiers:
            gene.qualifiers['trans_splicing'] = ['true']
    return top_record


def _core2tbl(top_record: SeqRecord):
    """Convert the to a tbl format. It should be used, and only used, if you
    want submit your annotaion to public database.

    !!NOTICE!!: Unlike other methods, the keep keys for feature table was stored in ABC.feature2tbl function.

    Args:
        top_record (): _description_

    Returns:
        : _description_
    """
    string_lst = [f'>Feature {top_record.id}']
    for gene in top_record.features:
        # non-gene features not output
        # In some extremely rare cases, there are pseudogenes whose type is pseudogene.
        # And since feature are generated from characters directly, the pseudogene and
        # trans_splicing need to set to [''], instead of [None] in gbk formats.
        if gene.type == 'pseudogene':
            gene.qualifiers['pseudo'] = ['']
        if gene.type == 'gene':
            if 'pseudo' in gene.qualifiers.keys():
                gene.qualifiers['pseudo'] = ['']
            if 'trans_splicing' in gene.qualifiers.keys():
                gene.qualifiers['trans_splicing'] = ['']
                for _subfeature in gene.sub_features:
                    _subfeature.qualifiers['trans_splicing'] = ['']
            string_lst.append(feature2tbl(gene))
            # Under extremely rare cases, pseudogenes don't have sub_features
            if gene.qualifiers.get('gene_biotype') == ['protein_coding'] and len(gene.sub_features) > 0:
                string_lst.append(feature2tbl(gene.gen_gb_cds()))
            else:
                # For tRNA and rRNA, no need to output exon for tbl format
                for rna in gene.sub_features:
                    string_lst.append(feature2tbl(rna))
    return '\n'.join(string_lst) + '\n'
