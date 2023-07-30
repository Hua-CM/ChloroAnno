# -*- coding: utf-8 -*-
# @Time : 2023/4/17 15:45
# @Author : Zhongyi Hua
# @FileName: ABC
# @Usage: The Basic class for this pipeline
# @Note:
# @E-mail: njbxhzy@hotmail.com
import json
import string
import random
from pathlib import Path
from copy import deepcopy

from Bio.SeqFeature import SeqFeature, CompoundLocation, FeatureLocation

file_location = (Path(__file__) /  '../..').resolve()
with open(file_location / 'ref/cds.json') as f_in:
    cds_dict = json.load(f_in)
with open(file_location / 'ref/rna.json') as f_in:
    rna_dict = json.load(f_in)
with open(file_location / 'ref/region.json') as f_in:
    region_set_list = [set(_) for _ in json.load(f_in).items()]
product_dict = {**cds_dict, **rna_dict}
standard_name_list = list(product_dict.keys()) + list(rna_dict.keys())


def fix_location(record_feature, plus_num):
    """Fix a subseq loaction problem.
    If you use position to subseq a child feature from a SeqRecord instance, the start
    location of the child feature would return to zero. So we should plus the real start position to it.

    Args:
        record_feature (SeqFeature): The sub_feature SeqFeature whose location need to fix.
        plus_num (int): The real start position
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


def id_generate(size=6, chars=string.ascii_uppercase):
    """The possible duplication between the old and new IDs can lead to conflicts 
    between the deletion of the old ID and the addition of the new ID. We solve 
    this problem by specifying a temporary ID.

    Args:
        size (int, optional): The tmp id length.
            Defaults to 6.
        chars (List[chr], optional): A list of characters used to generate tmp id. 
            Defaults to string.ascii_uppercase.

    Returns:
        fa_dict (dict): {new seq id: sequences[SeqRcord]}
        pd_gff  (pd.DataFrame): The ordered gff
    """
    return ''.join(random.choice(chars) for _ in range(size))


class MySeqFeature(SeqFeature):
    """This class was write for all features, while the CoreSeqFeature only 
    for 'gene'.
    """
    def __init__(self, _parent):
        self.__dict__ = deepcopy(_parent.__dict__)
        if not hasattr(self, 'sub_features'):
            self.sub_features = []
   
    def update_subfeature(self, _subfeature=None):
        """_summary_

        Args:
            _subfeature (_type_, optional): _description_. Defaults to None.
        """
        if _subfeature is not None:
            self.sub_features.append(_subfeature)
 
"""All possible qualifiers
GenBank
    gene
    codon_start(cds)
    transl_table(cds)
    product 
    locus_tag
GFF
    ID
    Name
    Parent
    gene_biotype: protein_coding, tRNA, rRNA (Gene)
    phase (CDS)
    product
Equivalent:
    locus_tag = ID(gene) | Parent(CDS, tRNA, rRNA)
    gene = Name
    codon_start = phase (not exactly identical)

**Notice**:
    1. exon and mRNA not necessary for plasmid genome annotation file. 
    2. Except the 'type', other attributes belong to SeqFeature is not critical.
    The values in qualifiers are critical
    3. These attributes do not need modification: seq
    4. {'pseudo': ['true']} in GFF format and {'pseudo': ['']} in GenBank format.
    Moreover, the pseudogene may don't have CDS sub features, and its type may be 'pseudogene'.
    The CoreSeqFeature don't hanlde this difference, which means you need handle it when output

Don't worry ID and product. The ID just be unique is OK, because the renumber method
in correct module will renumber the result. For product, correct name method in 
correct module will correct it as well.
"""

class CoreSeqFeature(MySeqFeature):
    """The CoreSeqFeature should always be the GENE. CDS, rRNA tRNA should be integrated
    into the CoreSeqFeature by `update_subfeatures`

    Args:
        MySeqFeature (_type_): _description_
    """
    ALL_KEY = {
        'gene': {'locus_tag', 'ID', 'gene', 'Name', 'gene_biotype'},
        'tRNA': {'locus_tag', 'ID', 'Parent', 'gene', 'Name', 'product'},
        'rRNA': {'locus_tag', 'ID', 'Parent', 'gene', 'Name', 'product'},
        'CDS':  {'locus_tag', 'ID', 'Parent', 'gene', 'Name', 'codon_start',
                 'phase', 'product', 'transl_table'}
    }

    def __init__(self, _parent):
        self.__dict__ = deepcopy(_parent.__dict__)
        if not hasattr(self, 'sub_features'):
            self.sub_features = []
        if self.type == 'gene':
            self._init_gene()

    def _init_gene(self) -> None:
        add_dict = {'locus_tag': 'ID',
                    'ID': 'locus_tag',
                    'gene': 'Name',
                    'Name': 'gene'}

        lost_keys = list(self.ALL_KEY.get('gene') - set(self.qualifiers.keys()))
        # The gene feature dose not have id
        for _key in lost_keys:
            try:
                # To get name
                if self.qualifiers.get('gene'):
                    gene_name = self.qualifiers.get('gene')[0]
                elif self.qualifiers.get('Name'):
                    gene_name = self.qualifiers.get('Name')[0]
                else:
                    gene_name = self.qualifiers.get('locus_tag')[0]

                if _key == 'gene_biotype':
                    if len(self.sub_features) != 0:
                        if self.sub_features[0].type in ['mRNA',' CDS']:
                            self.qualifiers['gene_biotype'] = ['protein_coding']
                        elif self.sub_features[0].type in ['tRNA', 'rRNA']:
                            self.qualifiers['gene_biotype'] = [self.sub_features[0].type]
                    else:
                        # To guess type
                        if gene_name in cds_dict:
                            gene_type = 'protein_coding'
                        elif gene_name.startswith('trn'):
                            gene_type = 'tRNA'
                        elif gene_name.startswith('rrn'):
                            gene_type = 'rRNA'
                        else:
                            print('Please check', gene_name, 'gene type. Just set to "protein_coding"')
                            gene_type = 'protein_coding'
                        self.qualifiers['gene_biotype'] = [gene_type]
                else:
                    self.qualifiers[_key] = self.qualifiers.get(add_dict.get(_key))
                # I've had enough of these SBs, could they tidy their files before upload to a public database?
                if ('Name' in lost_keys) and ('gene' in lost_keys):
                    self.qualifiers['Name'] = gene_name
                    self.qualifiers['gene'] = gene_name
                # If assign this at first, the last line would set both values to None
                if ('locus_tag' in lost_keys) and ('ID' in lost_keys):
                    self.qualifiers['locus_tag'] = [id_generate()]
                    self.qualifiers['ID']        = self.qualifiers['locus_tag']
            except KeyError:
                print(_key)
            except TypeError:
                print(_key)

    def _set_phase(self, ignore_strand=False):
        """set CDS phase in sub_feature of genes

        Args:
            _gene (SeqFeature) : A gene type SeqFeature, with CDS type sub_features
            gnore_strand (bool): For rps12, whose second part have two different directions CDS.

        Return:
            None
        """
        cds_length = 0
        # Note: ignore_strand should at the first beacause rps12 strand was None and will raise a TypeError
        if ignore_strand or self.strand > 0:
            _reverse = False
        else:
            _reverse = True
        self.sub_features.sort(key=lambda x: x.location.start, reverse=_reverse )
        for _cds in self.sub_features:
            _cds.qualifiers['phase'] = ['%d' % ((3-cds_length % 3) % 3)]
            cds_length += len(_cds.location)

    def _init_cds(self, _subfeature: SeqFeature) -> None:
        """Fill CDS qualifiers
        Args:
            _subfeature (SeqFeature): The CDS sub_feature for a gene.

        Return:
            cds_lst (List[SeqFeature]) or None: If there are multiple CDSs, return a list, otherwise, None.
        """
        lost_keys = list(self.ALL_KEY.get('CDS') - set(_subfeature.qualifiers.keys()))
        _subfeature.qualifiers['transl_table'] = ['11']
        # 'gene', 'Name', 'Parent', 'ID', and 'locus_tag' will be added when subfeatures add to subfeatures 
        # For GFF
        if 'codon_start' in lost_keys:
            _subfeature.qualifiers['codon_start'] = [str(int(_subfeature.qualifiers.get('phase')[0])+1)]
        # For GenBank
        if 'phase' in lost_keys:
            _subfeature.qualifiers['phase'] = [str(int(_subfeature.qualifiers.get('codon_start')[0])-1)]
        if isinstance(_subfeature.location, CompoundLocation):
            # It means that CDS subfeature is from the GenBank file and need to seperate
            cds_lst = []
            for _idx, cds_ele_location in enumerate(_subfeature.location.parts):
                _subcds = MySeqFeature(_subfeature)
                _subcds.qualifiers['ID'] = ['cds_' + self.qualifiers['ID'][0] + '_' + str(_idx+1)]
                _subcds.location = cds_ele_location
                cds_lst.append(_subcds)
            return cds_lst
        return None

    def _init_rna(self, _subfeature: SeqFeature) -> None:
        """Some tRNA subfeature has two exons. This function was written to tidy
         these tRNA subfeatures' exon subfeatures.

        Args:
            _subfeature (SeqFeature): A tRNA sub_feature for a gene.
        """
        # For GenBank input, the location of whose tRNA subfeatures with exon was CompoundLocation
        exon_count = 1
        # Just use gene prefix, not rRNA prefix
        exon_id = self.qualifiers.get('ID')[0]
        exon_lst = []
        # **IMPORTANCE:** Do not use MySeqFeature. The deepcopy was hard to control
        if isinstance(_subfeature.location, CompoundLocation):
            for exon_location in _subfeature.location.parts:
                exon_feature = SeqFeature(location = exon_location,
                                          type = 'exon',
                                          qualifiers=deepcopy(_subfeature.qualifiers))
                exon_feature.qualifiers['ID'] = [f'exon_{exon_id}_{exon_count:1d}']
                exon_feature.qualifiers['Parent'] = _subfeature.qualifiers.get('ID')
                exon_count += 1
                exon_lst.append(exon_feature)
            _subfeature.sub_features += exon_lst
        # For GFF with exon  subfeature input.
        elif len(_subfeature.sub_features) > 0:
            for exon_feature in _subfeature.sub_features:
                exon_feature.qualifiers['ID'] = [f'exon_{exon_id}_{exon_count:1d}']
                exon_feature.qualifiers['Parent'] = _subfeature.qualifiers.get('ID')
                exon_count += 1
        # For GenBank or GFF single-exon tRNAs that do not have sub_features. Add a exon feature
        else:
            exon_feature = SeqFeature(location = _subfeature.location,
                                      type = 'exon',
                                      qualifiers=deepcopy(_subfeature.qualifiers))
            exon_feature.qualifiers['ID'] = [f'exon_{exon_id}_{exon_count:1d}']
            exon_feature.qualifiers['Parent'] = _subfeature.qualifiers.get('ID')
            exon_lst.append(exon_feature)
            _subfeature.sub_features += exon_lst

    def _init_subfeature(self, _subfeature: SeqFeature) -> None:
        _subfeature.qualifiers['ID'] = [self.qualifiers['ID'][0] + '_' + id_generate()]
        _subfeature.qualifiers['Parent'] = self.qualifiers['ID']
        _subfeature.qualifiers['locus_tag'] = self.qualifiers['locus_tag']
        _subfeature.qualifiers['gene'] = self.qualifiers['gene']
        _subfeature.qualifiers['Name'] = self.qualifiers['Name']
        if 'product' not in self.qualifiers.keys():
            _subfeature.qualifiers['Name'] = [product_dict.get(_subfeature.qualifiers['gene'][0])]
    
    def gen_gb_cds(self):
        """The CDS stored in CoreSeqFeature was as a list of SeqFeature[SimpleLocation], 
        and thus, it need to convert to a SeqFeature[CompondLocation] for write the GenBank
        file

        Return:
            cds_feature: A combined, single CDS feature for GenBank output
        """
        if len(self.sub_features) == 0:
            return None
        if len(self.sub_features) == 1:
            return self.sub_features[0]
        simple_locations = [_.location for _ in self.sub_features]
        cds_feature = SeqFeature(location = CompoundLocation(simple_locations),
                                 type='CDS',
                                 qualifiers = self.sub_features[0].qualifiers)
        cds_feature.location = CompoundLocation(simple_locations)
        cds_feature.qualifiers['ID'] = ['cds_' + self.qualifiers['ID'][0]]
        return cds_feature

    def update_subfeature(self, _subfeature=None):
        """Need update ID | Parent before add it to the subfeature

        Args:
            _subfeature (_type_, optional): _description_. Defaults to None.
        """
        # For safe, all subfeatures change to MySeqFeature to ensure it has sub_feature attributes
        if not hasattr(_subfeature, 'sub_features'):
            _subfeature.sub_features = []
        if _subfeature is not None:
            self._init_subfeature(_subfeature)
            if _subfeature.type == 'CDS':
                other_cds = self._init_cds(_subfeature)
                if other_cds:
                    self.sub_features += other_cds
                    self._set_phase()
                else:
                    self.sub_features.append(_subfeature)
            else:
                self._init_rna(_subfeature)
                self.sub_features.append(_subfeature)

def feature2tbl(record: SeqFeature):
    """Convert a SeqFeature to tbl output record
    
    **NOTICE**
        1. "codon_start", "phase", 'translation' for CDS need to remove from the CDS seqfure
        2. !!IMPORTANCE!! The start position need to plus one to convert to 1-based index
        3. Not applicable for rps12.
    
    Format: A string
        6102\t6063\tCDS\t\t\n
        5164\t4938\t\t\t\n
		\t\t\tproduct\tribosomal protein S16
		\t\t\ttransl_table\t11
		\t\t\tprotein_id\tgb|QQO80509.1|

    Args:
        record (SeqFeature) : A GenBank-lik SeqFeature (CompoundLocation)

    Returns:
        record_string (str) : A string for tbl output
    """
    # Part1: Location
    part1_line_lst = []
    for _location in record.location.parts:
        if _location.strand > 0:
            part1_line_lst.append([str(_location.start+1), str(_location.end), '', '' ,''])
        else:
            part1_line_lst.append([str(_location.end), str(_location.start+1), '', '' ,''])
    part1_line_lst[0][2] = record.type

    # Part2: Qualifier features
    keep_keys = set(['transl_table', 'product', 'locus_tag', 'gene', 'codon_start', 'pseudo', 'trans_splicing'])
    if record.type == 'gene':
        keep_keys.remove('product')
    part2_lint_lst = [['' ,'', '', _key, _value[0]]  for _key, _value in record.qualifiers.items() if _key in keep_keys]

    # Combine
    record_string = '\n'.join(['\t'.join(line).rstrip() for line in part1_line_lst + part2_lint_lst])
    return record_string
