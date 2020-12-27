#! /usr/bin/env python3
import argparse
from spring_package.Modeller import createModel


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create a 3D model from HH-search results.')
    parser.add_argument('-a', '--a_hhr', help='First HHR target file result', required=True)
    parser.add_argument('-b', '--b_hhr', help='Second HHR target file result', required=True)
    parser.add_argument('-i', '--index', help='PDB Database Index file (ffindex)', required=True)
    parser.add_argument('-d', '--database', help='PDB Database files (ffdata)', required=True)
    parser.add_argument('-c', '--cross', help='PDB Cross Reference', required=True)
    parser.add_argument('-o', '--output', help='Output model file', required=True)
    parser.add_argument('-g', '--log', help='Log file', required=True)
    parser.add_argument('-we', '--wenergy', help='Weight Energy term', type=float, default=-0.01, required=False)
    parser.add_argument('-ms', '--minscore', help='Minimum min-Z score threshold', type=float, default=10.0, required=False)
    parser.add_argument('-mt', '--maxtries', help='Maximum number of templates', type=int, default=20, required=False)
    parser.add_argument('-mc', '--maxclashes', help='Maximum fraction of clashes', type=float, default=0.1, required=False)
    parser.add_argument('-sr', '--showtemplate', help='Add reference template to model structure', required=False, default="true")
    args = parser.parse_args()
    createModel(args)
