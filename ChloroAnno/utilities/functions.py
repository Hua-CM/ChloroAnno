# -*- coding: utf-8 -*-
# @File    :   functions.py
# @Time    :   2023/04/21 15:35:27
# @Author  :   Zhongyi Hua
# @Usage   :   The functions related to chloroplast annotation directly.
# @Note    :   
# @E-mail  :   njbxhzy@hotmail.com



"""
把各个部分的check任务全部拆出来, 这样方便反复使用
1. 保守基因情况
2. un normal name
3. 重复区域
4. renumber
5. cds check
"""
import re
from pathlib import Path
from collections import defaultdict
import subprocess as sp
from shutil import copy

from difflib import SequenceMatcher
import portion as pt
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, CompoundLocation, SimpleLocation
from Bio.SeqRecord import SeqRecord

from .ABC import standard_name_list, CoreSeqFeature


def check_cds(gene, genome_seq):
    """Check whether the CDS sequence is legal

    Args:
        gene (SeqFeature): A gene type SeqFeature with CDS sub_features
        genome_seq (Bio.Seq.Seq): The genomic sequence
    """
    check_seq = ''
    gene.sub_features.sort(key=lambda x: x.location.start, reverse=bool(gene.location.strand == -1))
    for cds in gene.sub_features:
        check_seq += cds.location.extract(genome_seq)
    check_seq.translate(table=11, cds=True)


def correct_name(gene_name):
    """Use a distance method to correct gene name.

    Args:
        gene_name (str): The error gene name, like psbA1

    Returns:
        gene_name (str) or 'incorrect': The corrected gene name.
    """
    _similarities = [SequenceMatcher(None, gene_name, _name2).quick_ratio() for _name2 in standard_name_list]
    index_max = max(range(len(_similarities)), key=_similarities.__getitem__)
    if _similarities[index_max] > 0.8:
        return standard_name_list[index_max]
    return 'incorrect'


class Check:
    """The class provide all check methods

    **NOTICE:** Check class ignore non-gene features

    """
    standard_name_list = standard_name_list

    def __init__(self, top_record):
        """Initiation

        Args:
            top_record (SeqRecord): A SeqRecord object has a list of CoreSeqFeature sub features in
             'features' attribute.
        """
        self.record = top_record

    def check_hs(self):
        """
        Check housekeeping gene: matK, rbcL
        :return: a message
        """
        gene_name_list = [gene.qualifiers.get('Name')[0] for gene in self.record.features if gene.type == 'gene']
        if ('matK' not in gene_name_list) and ('matk' not in gene_name_list):
            print('matK loss!')
        if ('rbcL' not in gene_name_list) and ('rbcl' not in gene_name_list):
            print('rbcL loss!')

    def check_name(self):
        """_summary_

        Returns:
            _type_: _description_
        """
        return_dict = {}
        for gene in self.record.features:
            if gene.type == 'gene':
                gene_name = gene.qualifiers.get('Name')[0]
                if not ((gene_name in Check.standard_name_list) or gene_name.startswith('orf')):
                    curated_name = correct_name(gene_name)
                    return_dict.setdefault(gene_name, curated_name)
                    print('check ' + gene_name + ' . Is it ' + curated_name + "?")
        return return_dict

    def check_region(self):
        """
        check duplicated gene region
        :return: messages
        """
        print('check overlap region start')
        region_list = []
        locus_list = []
        for gene in self.record.features:
            if gene.type == 'gene':
                region_list.append(pt.closed(int(gene.location.start), int(gene.location.end)))
                locus_list.append(gene.qualifiers.get('Name')[0])
        for i in range(len(region_list) - 1):
            if region_list[i].overlaps(region_list[i + 1]):
            #if set(locus_list[i:i+2]) not in region_set_list: # I forget why I write this line, just mask it for now.
                print(locus_list[i], region_list[i], ' and ', locus_list[i + 1], region_list[i + 1],
                      ' are overlap')
        print('check overlap region done')

    def renumber(self, prefix: str):
        """Renumber gene id from one to the last

        Args:
            prefix (str): The gene prefix. e.g. "AAA" in "AAA001"
        
        Return:
            gene_count (int): The total gene number
        """
        print('Renumber gene id')
        gene_count = 1
        for gene in self.record.features:
            if gene.type == 'gene':
                gene.id = f'{prefix}{gene_count:03d}'
                gene.qualifiers['ID'] = [gene.id]
                gene.qualifiers['locus_tag'] = gene.qualifiers['ID']
                sub_count = 0
                gene_count += 1
                gene_type = gene.qualifiers['gene_biotype'][0]
                for subfeature in gene.sub_features:
                    sub_count += 1
                    _sub_prefix = 'rna' if gene_type in ['tRNA', 'rRNA'] else 'cds'
                    subfeature.id = f'{_sub_prefix}_{gene.id}_{sub_count:1d}'
                    subfeature.qualifiers['ID'] = [subfeature.id]
                    subfeature.qualifiers['Parent'] = [gene.id]
                    subfeature.qualifiers['locus_tag'] = gene.qualifiers['ID']
                    if gene_type in ['tRNA', 'rRNA']:
                        exon_count = 0
                        for exon in subfeature.sub_features:
                            exon_count += 1
                            exon.id = f'exon_{gene.id}_{exon_count:1d}'
                            exon.qualifiers['ID'] = [exon.id]
                            exon.qualifiers['Parent'] = [subfeature.id]
        print('Renumber gene id done')
        return gene_count - 1

    def check_cds(self, seq_ins: Bio.Seq.Seq, add_pseudo=False):
        """check if whether the CDSs is illegal

        Args:
            seq_ins (Bio.Seq.Seq): The chloroplast genome sequence
            add_pseudo (bool, optional): Whether add a "pseudo" attribute to the qualifier. Defaults to False.
        """
        self.record.seq = seq_ins
        for gene in self.record.features:
            if gene.type == 'gene':
                # Only useful for gene with gene_biotype
                gene_type = gene.qualifiers['gene_biotype'][0]
                gene.id = gene.qualifiers['ID'][0]
                if gene.qualifiers.get('Name')[0] == 'rps12':
                    continue
                if gene_type == 'protein_coding':
                    try:
                        check_cds(gene, self.record.seq)
                    except Exception as e:
                        if add_pseudo:
                            gene.qualifiers['pseudo'] = ['true']
                        print(gene.id)
                        print(e)


class Correct:
    """Correct annotation results
    """
    def __init__(self, top_record1: SeqRecord, top_record2: SeqRecord, geo_seq: Bio.Seq.Seq):
        """Combine result from two sources, typically two softwares (e.g. PGA and CPGAVAS2),
        to validate the final results.
        
        **NOTICE**:
        1. The more reliable result should be record1 (e.g. PGA), while the more informative
           result should be record2 (e.g. GeSeq2 or CPGAVAS2).
        2. Only non-gene features (e.g. source, IR, misc_feature) in record1 will be retained.

        Args:
            top_record1 (SeqRecord): _description_
            top_record2 (SeqRecord): _description_
            seq (Bio.Seq.Seq): 
        """
        self.record1 = top_record1
        self.record2 = top_record2
        self.seq = geo_seq
        self.curated = []
        self.gene_dict = defaultdict()
 

    def _create_link(self):
        """I use the gene name (e.g. psbA) to link different file (e.g. GeSeq and PGA).
        Thus, I need create a dialect.

        Returns:
            dialect (dict): {'name': [gene1[SeqFeature], gene2, ....]}
        """
        for gene in self.record1.features:
            if gene.type == 'gene':
                self.gene_dict.setdefault(gene.qualifiers.get('Name')[0], []).append(gene)
            else:
                # Some features may be not genes (e.g. IR, source). I don't curate these features.
                self.curated.append(gene)
        # Drop out non-gene features
        self.record2.features = [ gene for gene in self.record2.features if gene.type == 'gene']


    @staticmethod
    def _cal_dis_(gene1, gene2):
        """For gene with multiple copies, particullary tRNA, I use the gene distance 
        to decide which gene is the one I want to replace.

        Args:
            gene1 (SeqFeature): _description_
            gene2 (SeqFeature): _description_

        Returns:
            Distance (float) : The distance between two genes.
        """
        return abs((gene1.location.start+gene1.location.end)/2 - (gene2.location.start+gene2.location.end)/2)

    def _query_gene_(self, gene: CoreSeqFeature):
        """Since I use gene name to make the link, multiple gene copies would cause 
        a problem to map. So, this function is to fix this problem.

        Args:
            gene (CoreSeqFeature): A gene feature from record2

        Returns:
            gene (CoreSeqFeature) or None: A gene feature from record1. It there is 
             no corresponding gene feature in record1, then return None.
        """
        _query_name = gene.qualifiers.get('Name')[0]
        _subject = self.gene_dict.get(_query_name)
        if _subject:
            if len(_subject) > 1:
                dis_lst = [self._cal_dis_(gene, _gene) for _gene in _subject]
                index_min = min(range(len(dis_lst)), key=dis_lst.__getitem__)
                return _subject[index_min]
            return _subject[0]
        return None

    def correct_records(self):
        """I take records from the more informative software results.
        If the record from more informative one was incorrect, then, use
        its coresponding record from record1 to replace it.
        """
        _tmp_curated = []
        for gene in self.record2.features:
            gene_type = gene.qualifiers['gene_biotype'][0]
            if gene.qualifiers.get('Name')[0] == 'rps12':
                continue
            if gene_type == 'protein_coding':
                try:
                    check_cds(gene, self.seq)
                    _tmp_curated.append(gene)
                except:
                    fixed_gene = self._query_gene_(gene)
                    if fixed_gene:
                        _tmp_curated.append(fixed_gene)
            else:  # rRNA and tRNA
                _tmp_curated.append(gene)
        # Sort and remove duplicated records in the final list.
        ## Duplicated records are always casued by two pseudo genes exsit in the record2.
        ## For example, two pseudo petB genes were annotated in record2 at 3383..4024 and 
        ##  91021..91662, while there is only one in record1 (3382..4023). Thus, after correct, there 
        ## will two 3382..4023 petB genes in the result list.
        _tmp_curated.sort(key=lambda x: x.location.start)
        _tmp_gene = None
        for _gene in _tmp_curated:
            if _gene == _tmp_gene:
                continue
            _tmp_gene = _gene
            self.curated.append(_gene)

    def add_trna(self):
        """Some tRNAs were not annotated by PGA, and thus I use record2 to make up it.
        I also want to check duplicated RNAs and only keep duplicated records in record1.
        """
        record1_trna_locations = pt.empty()
        for gene_name, gene_lst in self.gene_dict.items():
            if gene_name.startswith('trn'):
                for gene in gene_lst:
                    record1_trna_locations = record1_trna_locations.union(pt.closed(int(gene.location.start), int(gene.location.end)))
        record2_trnas = [gene for gene in self.record2.features if gene.qualifiers['Name'][0].startswith('trn')]
        for gene in self.gene_dict.keys():
            if gene.startswith('trn') and (gene not in record2_trnas):
                # Check whether region dulication
                for trna_gene in self.gene_dict.get(gene):
                    trna_gene_location = pt.closed(int(trna_gene.location.start), int(trna_gene.location.end))
                    if not trna_gene_location.overlaps(record1_trna_locations):
                        self.curated.append(trna_gene)
        self.curated.sort(key=lambda x: x.location.start)
    
    def correct(self):
        """The main interface

        Returns:
            top_record (SeqRecord): The SeqRecord contains CoreSeqFeature
        """
        self._create_link()
        self.correct_records()
        self.add_trna()
        return SeqRecord(id=self.record1.id, seq=Seq(''), features=self.curated)

"""rps12 feature
1. SeqFeature numbers
for one rps12 gene. One gene feature and one CDS feature in GenBank,
but two gene feature and three CDS feature in gff

2. trans_splicing in qualifiers
# trans_splicing in GenBank SeqFeature qualifiers
 {'trans_splicing': ['']}
# trans_splicing in GFF SeqFeature qualifiers
 {'exception': ['trans-splicing']}

3. part in qualifiers
# GenBank Gene SeqFeature do not have 'part'
 {'part': ['1']} or  {'part': ['2']}

**NOTICE**: I use a mixture type for store: one gene SeqFeature with multiple CDS features
, and change to gff format when output 
"""

def add_rps12(pga_gb: Path, prefix: str, num_start: int):
    """Add rps12 infomation to gff file using PGA GenBank

    Args:
        pga_gb    (Path): The path of PGA output GenBank file
        prefix    (str) : The gene prefix
        num_start (int) : From which number to index
    
    Return:
        feature_lst (List[CoreSeqFeature]): A list of two rps12 genes.
    """
    genome = SeqIO.read(pga_gb, 'gb')
    features_lst = []
    rps12_list = [ele for ele in genome.features if ele.type == 'CDS' and ele.qualifiers.get('gene') == ['rps12']]
    part1 = [part for part in rps12_list if len(part.location.parts) == 1][0]
    part2_lst = [part for part in rps12_list if len(part.location.parts) > 1]
    gene_count = num_start
    for part2 in part2_lst:
        gene_id = f'{prefix}{gene_count:03d}'
        gene_qualifiers = {
            'locus_tag': [gene_id],
            'ID': [gene_id],
            'gene': ['rps12'],
            'Name': ['rps12'],
            'gene_biotype': ['protein_coding'],
            'product': ['ribosomal protein S12'],
            'trans_splicing': None
        }

        gene_location = CompoundLocation([part1.location,
                                          SimpleLocation(part2.location.start,
                                                         part2.location.end,
                                                         part2.location.strand)])
        gene_feature = SeqFeature(location= gene_location,
                                  type='gene',
                                  qualifiers=gene_qualifiers,
        )
        gene_feature = CoreSeqFeature(gene_feature)

        # Part1 at first
        tmp_locations = part1.location + part2.location
        cds_id = 1
        for cds_location in tmp_locations.parts:
            cds_qualifiers = {
                'ID': [f'cds_{gene_id}_{cds_id:1d}'],
                'locus_tag': [gene_id],
                'Parent': [gene_id],
                'codon_start': ['1'],
                'transl_table':['11'],
                'gene': ['rps12'],
                'Name': ['rps12'],
                'gene_biotype': ['protein_coding'],
                'product': ['ribosomal protein S12'],
                'trans_splicing': None
            }
            cds_feaure = SeqFeature(location=cds_location,
                                    type='CDS',
                                    qualifiers=cds_qualifiers)
            gene_feature.update_subfeature(cds_feaure)
            cds_id += 1
        
        features_lst.append(gene_feature)
        gene_count += 1
    return features_lst


class CurateSequence:
    """This Class was written for curate the raw sequence before annotation
    1. The strand direction. Although none of public database specicfies a specific
     strand, it always assumes the minus strand for trnH as the default direction.
    2. The start position. Since the chloroplast is the circular genome. The start position
     in fasta is defined by people. It always assumes the trnH gene should at the head to 
     avoid any gene locates on the gap.
    3. Whether there were amubiguous nucleotides in sequences.
    """
    def __init__(self, _seq_path: Path, tmp_dir: Path):
        self.genome = _seq_path
        self.query1 = (Path(__file__) / '../../ref/trnH.fa').resolve()
        self.query2 = (Path(__file__) / '../../ref/protein.fa').resolve()
        self.tmp = tmp_dir

    def _blast_wrapper(self, genome_seq: Path):
        """The blast wrapper to do blast

        Args:
            genome_seq (Path): The fasta path
        """
        # makeblastdb
        copy(genome_seq, self.tmp / 'genome.fasta')
        mkdbcline = [
            "makeblastdb",
            "-in", str(self.tmp / 'genome.fasta'),
            "-dbtype", "nucl"]
        sp.run(mkdbcline, check=True)
        # trnH
        # blastn command (trnH)
        ncline = [
            "blastn",
            "-query", self.query1,
            "-db", str(self.tmp / 'genome.fasta'),
            "-evalue", "0.001",
            "-word_size", "11",
            "-outfmt", "6 qseqid sstart send sstrand evalue",
            "-out", str(self.tmp / 'blastn.res')
        ]
        # Execute blastn
        sp.run(ncline, check=True)

        # tblastn command (matK, rbcL and psbA)
        tncline = [
            "tblastn",
            "-query", self.query2,
            "-db", str(self.tmp / 'genome.fasta'),
            "-evalue", "0.001",
            "-outfmt", "6 qseqid sstart send sstrand evalue",
            "-out", str(self.tmp / 'tblastn.res')
        ]
        # Execute tblastn
        sp.run(tncline, check=True)
    
    def _curate(self, genome_seq):
        self._blast_wrapper(genome_seq)
        raw_seq = SeqIO.read(genome_seq, 'fasta')
        key_lst = ['start', 'end', 'strand', 'evalue']
        res_dict = defaultdict()
        with open(self.tmp / 'blastn.res', encoding='utf-8') as f1, open(self.tmp / 'tblastn.res', encoding='utf-8') as f2:
            for line in f1.read().splitlines() + f2.read().splitlines():
                res_dict.setdefault(line.split()[0], []).append(dict(zip(key_lst, line.split()[1:])))
        if 'trnH-GUG' not in res_dict or 'psbA' not in res_dict:
            return raw_seq, None
        return raw_seq, res_dict


    def curate(self):
        """Main curate function

        Returns:
            seq_fixed (SeqRecord): The fixed sequence
        """

        def _cal_dis_(gene1, gene2):
            return abs(
                (int(gene1.get('start')) + int(gene1.get('end'))) / 2 - (int(gene2.get('start')) + int(gene2.get('end'))) / 2
            )
        ## for trnH on gap
        ## # rotate the chloroplast circle genome on a certain length to
        ## #  move the trnH out off the gap and rerun this method
        raw_seq, res_dict = self._curate(self.genome)
        loop_count = 0
        while not res_dict:
            # random rotate
            part1, part2 = raw_seq.seq[0: 1000], raw_seq.seq[1000:]
            tmp_seq = part2 + part1
            tmp_seq  = SeqRecord(seq=tmp_seq, id=raw_seq.id, description='')
            SeqIO.write(tmp_seq, self.tmp / 'rotate.fasta','fasta')
            raw_seq, res_dict = self._curate(self.tmp / 'rotate.fasta')
            loop_count += 1
            if loop_count > 3:
                print('Please check this chloroplast genome manually. Could not find trnH or psbA on it.')
                return None
        # check housekeeping
        if ('matK' not in res_dict) and ('matk' not in res_dict):
            print('matK loss!')
        if ('rbcL' not in res_dict) and ('rbcl' not in res_dict):
            print('rbcL loss!')
        # check ambiguous nucleotide
        if re.compile('[^ATCGNatcgn]').findall(str(raw_seq.seq)):
            print('Genome contain invalid characters (not ATCGNatcgn)')
        # fix position
        ## for multiple psbA
        index_min = min(range(len(res_dict.get('psbA'))), key=lambda x: float(res_dict.get('psbA')[x].get('evalue')))
        psbA = res_dict.get('psbA')[index_min]
        ## for multiple trnH
        dis_lst = [_cal_dis_(_gene, psbA) for _gene in res_dict.get('trnH-GUG')]
        index_min = min(range(len(dis_lst)), key=dis_lst.__getitem__)
        trnh = res_dict.get('trnH-GUG')[index_min]
        if trnh.get('strand') == 'plus':
            _start = len(raw_seq) - int(trnh.get('start')) + 1 - 5
            raw_seq.seq = raw_seq.seq.reverse_complement()
        else:
            _start = int(trnh.get('start'))
        part1, part2 = raw_seq.seq[0: _start], raw_seq.seq[_start:]
        seq_fixed = part2 + part1
        seq_fixed = SeqRecord(seq=seq_fixed, id=raw_seq.id, description='')
        return seq_fixed


    def _iden_quartered_structure(self):
        """Identify LSC, SSR, and IRs.

        Returns:
            seq_LSC, seq_IRa, seq_SSC, seq_IRb (Bio.Seq.Seq):  
        """
        raw_seq = SeqIO.read(self.genome, 'fasta')
        # To avoid someting on gap
        double_seq = raw_seq + raw_seq
        SeqIO.write(double_seq, self.tmp / 'doule_seq.fasta', 'fasta')
        
        mkdbcline = [
            "makeblastdb",
            "-in", str(self.tmp / 'doule_seq.fasta'),
            "-dbtype", "nucl"
        ]
        sp.run(mkdbcline, check=True)

        # 2. blastn command
        ncline = [
            "blastn",
            "-query", str(self.tmp / 'doule_seq.fasta'),
            "-db", str(self.tmp / 'doule_seq.fasta'),
            "-strand", "both",
            "-perc_identity", "70",
            "-outfmt", "6 qlen length qstart qend sstart send",
            "-out", str(self.tmp / 'double_seq.res')
        ]
        sp.run(ncline, check=True)

        max_aln_n = 0
        min_ir = 1000
        origin_len = len(raw_seq)
        locations = set()

        with open(self.tmp / 'double_seq.res', encoding='utf-8') as f1:
            for _line in f1.read().splitlines():
                qlen, length, qstart, qend, sstart, send = [int(_) for _  in _line.split('\t')]
                if qlen < min_ir:
                    continue
                if length >= origin_len:
                    continue
                location = tuple(sorted([qstart, qend, sstart, send]))
                # hit across origin and repeat
                if location[-1] - location[0] > origin_len:
                    continue
                # self to self or self to repeat self
                if len(set(location)) != 4:
                    continue
                if length < max_aln_n:
                    continue
                else:
                    max_aln_n = length
                locations.add(location)
            locations = list(locations)
            locations.sort(key=lambda x: x[1]-x[0], reverse=True)
            locations.sort(key=lambda x: x[0])
            locations = locations[:2]
            # ira_start, ira_end, irb_start, irb_end
            a = slice(locations[0][0]-1, locations[0][1])
            b = slice(locations[0][1], locations[0][2]-1)
            c = slice(locations[0][2]-1, locations[0][3])
            d = slice(locations[1][1], locations[1][2]-1)
            if (locations[0][2]-locations[0][1]) > (
                    locations[1][2]-locations[1][1]):
                region_LSC = b
                region_IRa = c
                region_SSC = d
                region_IRb = a
            else:
                region_LSC = d
                region_IRa = a
                region_SSC = b
                region_IRb = c

        seq_LSC = double_seq[region_LSC]
        seq_IRa = double_seq[region_IRa]
        seq_SSC = double_seq[region_SSC]
        seq_IRb = double_seq[region_IRb]

        return seq_LSC, seq_IRa, seq_SSC, seq_IRb

    def iso(self, ref_path: Path):
        """Identify isomers 

        Args:
            ref_path (Path): reference genome path

        Returns:
            out_seq (SeqRecord): The selected path
        """
        seq_LSC, seq_IRa, seq_SSC, seq_IRb = self._iden_quartered_structure()

        SeqIO.write(seq_SSC, self.tmp / 'SSC.fa' , 'fasta' )
        
        copy(ref_path, self.tmp / 'ref_genome.fasta')
        #makeblastdb
        # 1. makeblastdb command (创建参考基因组数据库)
        mkdbcline = [
            "makeblastdb",
            "-in", str(self.tmp / 'ref_genome.fasta'),
            "-dbtype", "nucl"
        ]
        sp.run(mkdbcline, check=True)

        # 2. blastn command (执行SSC序列比对)
        ncline=[
            "blastn",
            "-query", str(self.tmp / 'SSC.fa'),
            "-db", str(self.tmp / 'ref_genome.fasta'),
            "-strand", "both",
            "-max_hsps", "1",
            "-outfmt", "6 qseqid sseqid sstrand evalue",
            "-out", str(self.tmp / 'SSC.res')
        ]
        sp.run(ncline, check=True)
        
        strand = 0
        for _line in Path(self.tmp/ 'SSC.res').read_text(encoding='utf-8').strip().split('\n'):
            _tmp_eles = _line.split('\t')
            if _tmp_eles[1] == 'plus':
                strand += 1
            else:
                strand -= 1
        #seqid1strand, seqid2strand, refstrand
        if strand > 0:
            out_seq = seq_LSC + seq_IRa + seq_SSC + seq_IRb
            out_seq.id = seq_LSC.id
            out_seq.description = ''
            return out_seq
        if strand < 0:
            out_seq = seq_LSC + seq_IRa + seq_SSC.reverse_complement() + seq_IRb
            out_seq.id = seq_LSC.id
            out_seq.description = ''
            return out_seq
        print('Please check', self.genome.stem, 'manually')
        return None
