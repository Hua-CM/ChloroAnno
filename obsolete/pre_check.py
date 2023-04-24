# -*- coding: utf-8 -*-
# @Time : 2020/7/11 22:44
# @Author : Zhongyi Hua
# @FileName: pre_check.py
# @Usage: chloroplast is a circle. This script will fix the start of sequence to trnH nearby
# @Note:
# @E-mail: njbxhzy@hotmail.com

"""
现在的脚本没法处理两个trnH的情况，准备补充下trnH前后有psbA的情况作为进一步判断
"""

import re
import argparse
from tempfile import mktemp
from shutil import copy, rmtree
from pathlib import Path
from collections import defaultdict

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbitblastnCommandline, NcbiblastnCommandline, NcbimakeblastdbCommandline





def main():
    """Main interface
    """
    parser = argparse.ArgumentParser(
        description='Check and fix the start of chloroplast genome')
    parser.add_argument('-i', '--input', required=True, type=Path, metavar='<File path>',
                        help='one sequence path per line')
    parser.add_argument('-o', '--output', required=True, type=Path, metavar='<Dir path>',
                        help='Output directory')
    args = parser.parse_args()
    if not args.output.exists():
        args.output.mkdir()
    with open(args.input) as f_in:
        sequence_lst = f_in.read_text().strip().split('\n')
    for row in sequence_lst:
        print(Path(row).name)
        main_ins = CurateSequence(row, args.output / Path(row).name)
        main_ins.curate()


if __name__ == '__main__':
    main()
