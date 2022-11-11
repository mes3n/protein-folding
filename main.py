import numpy as np
import molmod as mm  # load 


'''
Fold protein based of hydrogen bonds, sulfide bridges and charges
'''

class S:  # enumerator for the potential atomic symbols
    H = "hydrogen"
    O = "oxygen"
    N = "nitrogen"
    S = "sulfide"

class AA:
    pass

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
        self.amino_acids: list[AminoAcid] = amino_acids


def main():
    TRP_CAGE_SEQ = "NLYIQWLKDGGPSSGRPPPS"

    print('[')
    for c in TRP_CAGE_SEQ:
        print(f'{c}, ', end='')
    print(']')

    return


if __name__ == '__main__':
    main()
