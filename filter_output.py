#!/usr/bin/env python3

import sys
import re
from collections import defaultdict

# Catalytic residue functional atoms definitions
RESIDUE_DEFINITIONS = { 
    'Main': {'CA': [100, ''], 'C': [100, ''], 'O': [100, '']},
    'Ptm': {'CA': [100, ''], 'C': [100, ''],  'CB': [100, '']},
    'Ala': {'C' : [0, ''], 'CA' : [0, ''], 'CB' : [0, '']},
    'Cys': {'CA' : [0, ''], 'CB' : [0, ''], 'SG' : [0, '']},
    'Asp': {'CG': [3, 'E'], 'OD1': [3, 'E'], 'OD2': [3, 'E']},
    'Glu': {'CD': [3, 'D'], 'OE1': [3, 'D'], 'OE2': [3, 'D']},
    'Phe': {'CE1': [3, ''], 'CE2': [3, ''], 'CZ': [0, '']},
    'Gly': {'CA' : [0, ''], 'C' : [0, ''], 'O' : [0, '']},
    'His': {'CG': [0, ''], 'CD2':[8, ''], 'ND1': [8, '']},
    'Ile': {'CA': [0, 'VL'], 'CB': [0, 'VL'], 'CG1': [3, 'VL']},
    'Lys': {'CD' : [0, ''], 'CE' : [0, ''], 'NZ' : [0, '']},
    'Leu': {'CA': [0, 'VI'], 'CB': [0, 'VI'], 'CG': [3, 'VI']},
    'Met': {'CG' : [0, ''], 'SD' : [0, ''], 'CE' : [0, '']},
    'Asn': {'CG': [3, 'Q'], 'OD1': [1, 'Q'] , 'ND2': [1, 'Q']},
    'Pro': {'C' : [0, ''], 'CA' : [0, ''], 'O' : [0, '']},
    'Gln': {'CD': [3, 'N'], 'OE1': [1, 'N'], 'NE2': [1, 'N']},
    'Arg': {'CZ': [0, ''], 'NH1': [3, ''], 'NH2': [3, '']},
    'Ser': {'CA': [0, 'TY'], 'CB':[0, 'TY'], 'OG': [1, 'TY']},
    'Thr': {'CA': [0, 'SY'], 'CB': [0, 'SY'], 'OG1': [1, 'SY']},
    'Trp': {'NE1' : [0, ''], 'CZ2' : [0, ''], 'CH2' : [0, '']},
    'Tyr': {'CE1': [3, 'ST'], 'CZ': [3, 'ST'], 'OH': [1, 'ST']},
    'Val': {'CA': [0, 'LI'], 'CB': [0, 'LI'], 'CG': [3, 'LI']},
    }


class Match():
    """Jess match"""

    def __init__(self, target='', template='', rmsd=None, score=None, remark_string='REMARK MATCH'):
        self.target = target
        self.template = template
        self.rmsd = rmsd
        self.score = float(score)
        self.remark_string = remark_string
        self.atoms = []
        self.chains = set()

    def __repr__(self):
        return '{}\n{}\nENDMDL\n'.format(self.remark_string, '\n'.join(self.atoms))

    @property
    def id(self):
        return '{};{};{}'.format(target, template, ':'.join(list(chains)))

    @classmethod
    def from_remark(cls, remark_string):
        fields = remark_string.strip().split()
        return cls(fields[1], fields[2], fields[3], fields[-1], remark_string=remark_string)

    def add_atom(self, line):
        self.atoms.append(line)
        self.chains.add(line[20:22].strip())

    def is_from_same_structure(self, other):
        return (self.target == other.target and self.template == other.template)

    def is_from_same_chains(self, other):
        return self.chains == other.chains

def get_best_match(matches):
    min_score = 999
    best_matches = []
    for match in matches:
        if match.score < min_score:
            min_score = match.score
            best_match = match
    return best_match

def main(all_chains=False):

    prev_match = None
    matches = defaultdict(list)
    unique_matches = []
    for line in sys.stdin:
        line=line.strip()
        if line.startswith('REMARK'):
            match = Match.from_remark(line)
        if line.startswith('ATOM') or line.startswith('HETATM'):
            match.add_atom(line)
        if line.startswith('ENDMDL'):
            if prev_match and not match.is_from_same_structure(prev_match):
                for chain_matches in matches.values():
                    best_match = get_best_match(chain_matches)
                    unique_matches.append(best_match)
                if all_chains:
                    for unique_match in unique_matches:
                        print(unique_match)
                else:
                    print(get_best_match(unique_matches))
                matches = defaultdict(list)
                unique_matches = []
            matches[':'.join(list(match.chains))].append(match)
            prev_match = match


if __name__ == '__main__':
    main()
