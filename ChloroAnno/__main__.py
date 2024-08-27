# -*- coding: utf-8 -*-
# @File    :   __main__.py
# @Time    :   2023/04/19 02:02:40
# @Author  :   Zhongyi Hua 
# @Version :   0.0.1
# @E-mail  :   njbxhzy@hotmail.com
# @Note    :   None
# @Usage   :   main interface

import argparse
from pathlib import Path
from tempfile import mktemp
from shutil import rmtree

from Bio import SeqIO
from tqdm import tqdm

from task import correct, curate, isomerism, check, tidy_geseq, convert
from utilities import parse_meta, dft, write_anno

def parse_args():
    """Parse arguments
    """
    # command line parser
    parser = argparse.ArgumentParser(
        description='A toolkit for Chloroplast Annotation')
    parser.add_argument('-t', '--task', type=str, metavar='<choice>', required=True,
                        choices=['check', 'correct', 'convert', 'curate', 'iso', 'tidy'],
                        help='Task type. Refer to README for more details')
    parser.add_argument('-i', '--info', required=True, type=Path, metavar='<File path>',
                        help='Meta info table. Refer to README for more details')
    parser.add_argument('--outfmt', type=str, metavar='<choice>', default='gb',
                        choices=['gb', 'gff', 'tbl'],
                        help='Output format. Default: gb')
    parser.add_argument('-o', '--output', type=Path, metavar='<Dir path>', default=Path('ChloroAnnoOutput'),
                        help='Output directory')
    parser.add_argument('--tmp', type=Path, metavar='<Dir path>', default=Path(mktemp()),
                        help='Specify the tempoary directory manually. Default: The system temporary directory')
    args = parser.parse_args()
    return args


def main(args):
    """Main interface
    """
    sample_lst = parse_meta(args.info)
    tmp = args.tmp
    if not args.output.exists():
        args.output.mkdir()
    if args.task == 'iso':
        if not tmp.exists():
            tmp.mkdir()
        for _sample in tqdm(sample_lst):
            print(f'start check {_sample["inpath1"].stem} isomerism')
            selected_path = isomerism(_sample['inpath1'], _sample['refpath'], tmp)
            if selected_path:
                SeqIO.write(selected_path, Path(args.output) / Path(_sample['inpath1']).name, 'fasta')
        rmtree(tmp)
    elif args.task == 'curate':
        if not tmp.exists():
            tmp.mkdir()
        for _sample in tqdm(sample_lst):
            print(f'start curate {_sample["inpath1"].stem} sequence')
            seq_fixed = curate(_sample['inpath1'], tmp)
            if seq_fixed:
                SeqIO.write(seq_fixed, Path(args.output) / Path(_sample['inpath1']).name, 'fasta')
            else:
                continue
        rmtree(tmp)
    elif args.task == 'correct':
        for _sample in tqdm(sample_lst):
            try:
                print(f'start correct {_sample["inpath1"].stem} annotation')
                out_record = correct(_sample['inpath1'], _sample['inpath2'], _sample['prefix'])
                if _sample.get('organism'):
                    out_record.annotations['organism'] = _sample.get('organism')
                write_anno(out_record,
                        Path(args.output) / (Path(_sample['inpath1']).stem + '.' + args.outfmt),
                        args.outfmt)
            except:
                print(f'Correct {_sample["inpath1"].stem} annotation fail')
                continue
    elif args.task == 'check':
        for _sample in tqdm(sample_lst):
            if 'inpath2' in _sample:
                check(_sample['inpath1'], _sample['inpath2'])
            else:
                check(_sample['inpath1'])
    elif args.task == 'tidy':
        for _sample in tqdm(sample_lst):
            out_record = tidy_geseq(_sample['inpath1'])
            write_anno(out_record,
                       Path(args.output) / Path(_sample['inpath1']).name,
                       'gff')
    elif args.task == 'convert':
        for _sample in tqdm(sample_lst):
            if 'informat1' in _sample:
                f_type = _sample['informat1']
            else:
                f_type = dft( _sample['inpath1'])
                if not f_type:
                    raise TypeError('Please specify the input file format using "informat1"')
            record = convert(_sample['inpath1'], f_type)
            write_anno(record,
                       Path(args.output) / (Path(_sample['inpath1']).stem + '.' + args.outfmt),
                       args.outfmt)


if __name__ == '__main__':
    main(parse_args())
