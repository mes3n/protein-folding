import numpy as np
import molmod as mm  # load 

from enum import Enum


'''
Fold protein based of hydrogen bonds, sulfide bridges and charges
'''

class S(Enum):  # enumerator for atomic symbols
    H = 'hydrogen'
    O = 'oxygen'
    N = 'nitrogen'
    S = 'sulfide'

class AA(Enum):  # enumerator for amino acids
    A = 'alanine'
    B =	'asparagine'
    C =	'cysteine'
    D =	'aspartic acid'
    E =	'glutamic acid'
    F =	'phenylalanine'
    G =	'glycine'
    H =	'histidine'
    I =	'isoleucine'
    K =	'lysine'
    L =	'leucine'
    M =	'methionine'
    N =	'asparagine'
    P =	'proline'
    Q =	'glutamine'
    R =	'arginine'
    S =	'serine'
    T =	'threonine'
    U = 'selenocysteine'
    V =	'valine'
    W =	'tryptophan'
    Y =	'tyrosine'


class Atom:
    def __init__(self, position: np.array, symbol: str, *args, **kwargs):
        self.position: np.array = position
        self.symbol: str = symbol

        
# class AminoAcid:
#     def __init__(self, position: np.array, side_chain: list[Atom], *args, **kwargs):
#         self.chains: dict[str, list[Atom]] = {
#             'main': main_chain, 'side': side_chain
#         }
#         self.position: np.array = position  # position of the amino acids central carbon atom


class Protein:
    def __init__(self, sequence: str, *args, **kwargs):
        
        self.amino_acids = [getattr(AA, c) for c in sequence]

        # self.amino_acids: list[AminoAcid] = amino_acids


def main():
    TRP_CAGE_SEQ = "NLYIQWLKDGGPSSGRPPPS"

    trp_cage = Protein(TRP_CAGE_SEQ)
    print(trp_cage.amino_acids)

    for aa in AA:
        if aa not in trp_cage.amino_acids:
            print(aa)

    return


if __name__ == '__main__':
    main()
