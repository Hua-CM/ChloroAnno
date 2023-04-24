# -*- coding: utf-8 -*-
# @Time : 2020/4/16 20:59
# @Author : Zhongyi Hua
# @FileName: correct2.py
# @Usage: correct annotation error using CPGAVAS2 and PGA result.
# @Note:
# @E-mail: njbxhzy@hotmail.com

from pathlib import Path
from collections import defaultdict
from datetime import date

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline
from utilities import Check, Correct, CurateSequence, add_rps12, read_geseq
from utilities import read_anno
#from utilities import continue_on_error

# !!NOTICE!! All methods in this scirpt can only be used in a loop because they are under continue decorator.

#@continue_on_error
# This is the old version !!
def isomerism_old(asm_path1:Path, asm_path2:Path, ref_path: Path, tmp_dir:Path):
    """_summary_

    Args:
        asm_path1 (Path): _description_
        asm_path2 (Path): _description_
        ref_path (Path): _description_
        tmpdir (Path): A temporary directory

    Returns:
        isomer (Path): The selected isomer path
    """
    path1_seq = SeqIO.read(asm_path1, 'fasta')
    path2_seq = SeqIO.read(asm_path2, 'fasta')
    ref_seq = SeqIO.read(ref_path, 'fasta')
    asm_dict = {path1_seq.id: asm_path1,
                path2_seq.id: asm_path2,
                ref_seq.id:   ref_path}
    query_path = (Path(__file__) / '../ref/ndh.fa').resolve()
    subject_path = tmp_dir / 'subject.fasta'

    SeqIO.write([path1_seq, path2_seq, ref_seq], subject_path, 'fasta')

    dbcmd = NcbimakeblastdbCommandline(input_file=subject_path, dbtype='nucl')
    qcmd =  NcbiblastnCommandline(
        query=query_path,
        db=subject_path,
        evalue=0.001,
        out=tmp_dir / 'blastn.res',
        outfmt="6 qseqid sseqid sstrand evalue")

    dbcmd()
    qcmd()

    #{seqid1: 0, seqid2: 0, refid: 0}
    # if one gene on plus strand, add 1, otherwise minus 1
    tmpstrand_dict = defaultdict(int)
    for _line in Path(tmp_dir / 'blastn.res').read_text(encoding='utf-8').strip().split('\n'):
        _tmp_eles = _line.split('\t')
        if _tmp_eles[1] == 'plus':
            tmpstrand_dict[_tmp_eles[1]] += 1
        tmpstrand_dict[_tmp_eles[1]] -= 1
    #seqid1strand, seqid2strand, refstrand
    if tmpstrand_dict.get(path1_seq.id) == tmpstrand_dict.get(ref_seq.id):
        print(asm_path1)
        return asm_dict.get(path1_seq.id)
    if tmpstrand_dict.get(path2_seq.id) == tmpstrand_dict.get(ref_seq.id):
        print(asm_path2)
        return asm_dict.get(path2_seq.id)
    print('Please check', asm_path1, 'manually')
    return None

def isomerism(fasta_path:Path, ref_path: Path, tmp_dir:Path):
    """Select the right path from the chloroplast paths

    Args:
        fasta_path (Path): The input fasta path
        ref_path (Path): A 
        tmp_dir (Path): A temporary directory path

    Returns:
         selected_path (Bio.SeqRecord): The selected sequence
    """
    main_ins = CurateSequence(fasta_path, tmp_dir)
    selected_path = main_ins.iso(ref_path)
    return selected_path

#@continue_on_error
def curate(fasta_path: Path, tmp_dir: Path):
    """ Check and curate sequence before annotation

    Args:
        fasta_path (Path): The input fasta path
        tmp_dir (Path): A temporary directory path

    Returns:
        seq_fixed (Bio.Seq.Seq): The fixed sequence
    """
    main_ins = CurateSequence(fasta_path, tmp_dir)
    seq_fixed = main_ins.curate()
    return seq_fixed

#@continue_on_error
def tidy_geseq(geseq_gff_path: Path):
    """Tidy gff output

    1. Remove "note" and "gbkey" in GeSeq gff file qualifiers.
    2. Add CDS and tRNA to correct gene parent
    
    Args:
        geseq_gff_path (Path): _description_
    """
    raw_record = read_geseq(geseq_gff_path)
    return raw_record
 
#@continue_on_error
def correct(pga_path:Path, cpgavas2_path: Path, prefix: str):
    """Correct chloroplast genome annotation using results from two softwares

    Args:
        pga_path (Path): _description_
        cpgavas2_path (Path): _description_
    
    Return:
        renumbered_record (SeqRecord): A SeqRecord object with CoreSeqFeature subfeatures
    """
    record1 = read_anno(pga_path, 'gb')
    record2 = read_anno(cpgavas2_path, 'gb')
    correct_ins = Correct(record1, record2, record1.seq)
    corrected_record = correct_ins.correct()
    # Set id for check report
    corrected_record.id = record1.id
    # Check corrected gene
    check_ins = Check(corrected_record)
    # Renumber first to ensure the gene id in check report is in line with the result
    total_gene_num = check_ins.renumber(prefix)
    check_ins.check_name()
    check_ins.check_region()
    check_ins.check_cds(record1.seq, add_pseudo=True)
    renumbered_record = check_ins.record
    # Add rps12
    rps12_lst = add_rps12(pga_path, prefix, total_gene_num+1)
    renumbered_record.features += rps12_lst
    # For safe, return a SeqRecord with annotations
    renumbered_record.annotations = {
            'molecule_type': 'DNA',
            'topology': 'circular',
            'data_file_division': 'PLN',
            'date': date.today().strftime('%d-%b-%Y'),
            'accessions': [renumbered_record.id],
            'sequence_version': 1,
            'keywords': [''],
        }
    return renumbered_record

#@continue_on_error
def check(input_path: Path, seq_path=None):
    """Check annotatiosn validity 

    Args:
        input_path (Path): The GenBank / GFF file inputs.
        seq_path (path, optional): If check the gff file, a fasta sequence is necessary
    """
    record1 = read_anno(input_path, 'gb')
    if seq_ins:
        seq_ins = SeqIO.read(seq_path, 'fasta').seq
    else:
        seq_ins = record1.seq
    tmp_check = Check(record1)
    tmp_check.check_cds(seq_ins)
    tmp_check.check_name()
    tmp_check.check_region()

#@continue_on_error
def convert(input_path: Path, ftype):
    """Convert file input to a SeqRecord with CoreSeqFeature

    Args:
        input_path (Path): Input file path
        ftype (str): File type

    Returns:
        record (SeqRecord):  A SeqRecord with CoreSeqFeature
    """
    record = read_anno(input_path, ftype)
    return record
