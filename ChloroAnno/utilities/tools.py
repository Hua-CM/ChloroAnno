# -*- coding: utf-8 -*-
# @File    :   tools.py
# @Time    :   2023/04/21 15:35:19
# @Author  :   Zhongyi Hua
# @Usage   :   The fucntion not related to chloroplast annotation directly. such as determing file format
# @Note    :   
# @E-mail  :   njbxhzy@hotmail.com

from pathlib import Path


def dft(file_path: Path):
    """Determine file type

    Args:
        file_path (Path): _description_

    Returns:
        file_type (str): The file type could be used by read_anno / write_anno 
    """
    if file_path.suffix in ['.gb', '.gbf', 'genbank']:
        return 'gb'
    elif file_path.suffix in ['.gff', '.gff3']:
        return 'gff'
    elif file_path.suffix in ['.tbl']:
        return 'tbl'
    elif file_path.suffix in ['.fa', '.fasta', '.fna']:
        return 'fasta'
    else:
        return None


def parse_meta(file_path: Path):
    """Parse meta info file

    Args:
        file_path (Path): The 

    Raises:
        ValueError: There are illegal meta file headers

    Returns:
        info_lst (list): [{'inpath1': xxx, 'refpath': xxx, ....}, {'inpath1': xxx, 'refpath': xxx, ....}]
    """
    init_method = {'inpath1': Path, 'inpath2': Path, 'refpath': Path, 'organism': str,
                   'prefix': str, 'informat1': str,  'infomart2':str }
    lines = file_path.read_text().strip().split('\n')
    name_lst = lines[0].split('\t')
    if set(name_lst) - set(init_method.keys()):
        raise ValueError('Illegal meta file header')
    info_lst = []
    for _line in lines[1:]:
        info_lst.append({name_lst[_idx]: init_method.get(name_lst[_idx])(_value) for _idx, _value in enumerate(_line.split('\t'))})
    return info_lst

def continue_on_error(func):
    """To make the loop could run when encounter errors

    Args:
        func (_type_): _description_
    """
    def wrapper(*args, **kwargs):
        while True:
            try:
                return func(*args, **kwargs)
            except Exception as e:
                print(e)
                continue
    return wrapper
