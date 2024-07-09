# -*- coding:utf-8 -*-

import argparse, os
from chem_utils import smi2mol, Get_group
import csv


script_dir = os.path.dirname(os.path.abspath(__file__))


def CalcLOI(smi, coef_path, intercept, out_formula=False, digits=None):
    with open(coef_path, 'r') as fin:
        lines = fin.read().strip().split('\n')

    condensed_coef_dict = {}
    gas_coef_dict = {}
    for line in lines:
        group, condensed_coef, gas_coef = line.strip().split('\t')
        condensed_coef_dict[group] = float(condensed_coef)
        gas_coef_dict[group] = float(gas_coef)

    mol = smi2mol(smi)
    group_num = Get_group(mol, list(condensed_coef_dict.keys()))
    condensed_contri = sum([group_num[group] * condensed_coef_dict[group] for group in group_num])
    gas_contri = sum([group_num[group] * gas_coef_dict[group] for group in group_num])
    LOI = intercept + condensed_contri + gas_contri 

    if digits:
        LOI = round(LOI, digits)

    if out_formula:
        formula = f"{LOI} = {intercept} + " + ' + '.join([f'{group} * {condensed_contri[group]}' for group in group_num]) 
        + ' + ' + ' + '.join([f'{group} * {gas_contri[group]}' for group in group_num])
        print("LOI formula: ", formula)
        out_path = os.path.join(script_dir, 'LOI_formula.csv')
        with open(out_path, 'w') as fout:
            writer = csv.writer(fout)
            writer.writerow(['LOI', LOI])
            writer.writerow(['group', 'condensed_coef', 'gas_coef'])
            for group in group_num:
                writer.writerow([group, condensed_coef_dict[group], gas_coef_dict[group]])

    return LOI

def CalcRFF(smi, intercept=20.95, out_formula=None, digits=None):
    FR_coef_path = os.path.join(script_dir, 'FR_coef.txt')
    return CalcLOI(smi, FR_coef_path, intercept, out_formula=out_formula, digits=digits)

def parse():
    parser = argparse.ArgumentParser(description='Calculate LOI for polymer repeat unit')
    parser.add_argument('--repeat_unit_SMILES', type=str, default='O=C(C1=CC=C(C(OCCO[*])=O)C=C1)[*]',
                        help='The SMILES string of monomer(example)')
    parser.add_argument('--out_formula', action='store_true',
                        help='Print LOI formula')
    parser.add_argument('--digits', type=int, default=None,
                        help='LOI digits')
    return parser.parse_args()

def main():
    args = parse()
    smi = args.repeat_unit_SMILES
    CalcRFF(smi, out_formula=args.out_formula, digits=args.digits)

if __name__ == '__main__':
    main()

