# -*- coding: utf-8 -*-
# @Time : 2023/4/16 18:32
# @Author : Zhongyi Hua
# @FileName: isomerism.py
# @Usage: Check the isomerism
# @Note:
# @E-mail: njbxhzy@hotmail.com

"""总的思想是通过SSC区域进行BLAST,所以核心问题是如何提取SSC区域。
解题思路:提取SSC区域的基因进行BLAST
进一步思路: 在数据库里存一个基因，分别去blast query sequence和ref sequence，从而解决问题
大部分ndh基因位于SSC区域，准备使用这部分基因进行
"""

from pathlib import Path
from tempfile import mktemp
from shutil import copy, rmtree
from collections import defaultdict
import argparse

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline

class CurateSequence:
    """_summary_
    """
    def __init__(self, _seq_path, _out_path):
        self.genome = _seq_path
        self.query1 = (Path(__file__) / '../ref/trnH.fa').resolve()
        self.query2 = (Path(__file__) / '../ref/protein.fa').resolve()
        self.output = _out_path
        self.tmp = Path(mktemp())

    def _make_db(self):
        self.tmp.mkdir()
        copy(self.genome, self.tmp / 'genome.fa')
        self.genome = self.tmp / 'genome.fa'
        cline = NcbimakeblastdbCommandline(input_file=self.genome, dbtype='nucl')
        cline()




def main():
    """Main interface
    """
    parser = argparse.ArgumentParser(
        description='Check and fix the start of chloroplast genome')
    parser.add_argument('-i', '--input', required=True, type=Path, metavar='<File path>',
                        help='Three column per line: path1 fasta, path2 fasta, and reference fasta')
    parser.add_argument('-o', '--output', required=True, type=Path, metavar='<Dir path>',
                        help='Output directory')
    args = parser.parse_args()

    if not args.output.exists():
        args.output.mkdir()
    sequence_lst = args.input.read_text().strip().split('\n')
    tmpdir = Path(mktemp())
    tmpdir.mkdir()
    for _row in sequence_lst:
        asm_path1, asm_path2, ref_path = _row.split('\t')
        choosed_path = check_iso(Path(asm_path1), Path(asm_path2), Path(ref_path), tmpdir)
        copy(choosed_path, args.output / choosed_path.name)


if __name__ == '__main__':
    main()
