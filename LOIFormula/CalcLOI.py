# -*- coding:utf-8 -*-

import argparse, os
from chem_utils1 import Get_vocab, FRstr


script_dir = os.path.dirname(os.path.abspath(__file__))


def CalcLOI(smi, CR, coef_path, CR_coef, intercept, out_formula=False, digits=None):
    with open(coef_path, 'r') as fin:
        lines = fin.read().strip().split('\n')

    coef_dict = {}
    for line in lines:
        vocab, coef = line.strip().split('\t')
        coef_dict[vocab] = float(coef)

    vocab_num = Get_vocab(smi, list(coef_dict.keys()))
    LOI = intercept + CR_coef * float(CR) + sum([vocab_num[vocab] * coef_dict[vocab] for vocab in vocab_num])

    if digits:
        LOI = round(LOI, digits)

    if out_formula:
        formula = f'{LOI} = {intercept} + {CR_coef} * {CR} + ' + ' + '.join(
            [f'{vocab_num[vocab]} * {coef_dict[vocab]}' for vocab in vocab_num])
        print("LOI formula: ", formula)

    return LOI

def CalcNoTFRF(smi, CR, CR_coef=0.07, intercept=20.95, out_formula=False, digits=None):
    NoFR_coef_path = os.path.join(script_dir, 'NoFR_coef.txt')
    return CalcLOI(smi, CR, NoFR_coef_path, CR_coef, intercept, out_formula=out_formula, digits=digits)

def CalcTFRF(smi, CR, CR_coef=0.14, intercept=20.95, out_formula=False, digits=None):
    FR_coef_path = os.path.join(script_dir, 'FR_coef.txt')
    return CalcLOI(smi, CR, FR_coef_path, CR_coef, intercept, out_formula=out_formula, digits=digits)

def parse():
    parser = argparse.ArgumentParser(description='Calculate LOI for polymer repeat unit')
    parser.add_argument('--repeat_unit_SMILES', type=str, default='O=C(C1=CC=C(C(OCCO[*])=O)C=C1)[*]',
                        help='The SMILES string of monomer(example)')
    parser.add_argument('--CR', type=float, default=0.0,
                        help='The CR value (> 600°C)')
    parser.add_argument('--out_formula', action='store_true',
                        help='Print LOI formula')
    parser.add_argument('--digits', type=int, default=None,
                        help='LOI digits')
    return parser.parse_args()

def main():
    args = parse()
    smi = args.repeat_unit_SMILES
    CR = args.CR
    if FRstr(smi):
        CalcNoTFRF(smi, CR, out_formula=args.out_formula, digits=args.digits)
    else:
        CalcTFRF(smi, CR, out_formula=args.out_formula, digits=args.digits)

if __name__ == '__main__':
    main()

