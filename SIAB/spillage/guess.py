def _elem_z(elem):
    """get the atomic number of the element
    
    Parameters
    ----------
    elem: str
        the element symbol
    
    Returns
    -------
    int: the atomic number
    """
    z = {"H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10,
         "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18, "K": 19,
         "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28,
         "Cu": 29, "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36, "Rb": 37,
         "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46,
         "Ag": 47, "Cd": 48, "In": 49, "Sn": 50, "Sb": 51, "Te": 52, "I": 53, "Xe": 54, "Cs": 55,
         "Ba": 56, "La": 57, "Ce": 58, "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64,
         "Tb": 65, "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71, "Hf": 72, "Ta": 73,
         "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80, "Tl": 81, "Pb": 82,
         "Bi": 83, "Po": 84, "At": 85, "Rn": 86, "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90, "Pa": 91,
         "U": 92, "Np": 93, "Pu": 94, "Am": 95, "Cm": 96, "Bk": 97, "Cf": 98, "Es": 99, "Fm": 100,
         "Md": 101, "No": 102, "Lr": 103, "Rf": 104, "Db": 105, "Sg": 106, "Bh": 107, "Hs": 108,
         "Mt": 109, "Ds": 110, "Rg": 111, "Cn": 112, "Nh": 113, "Fl": 114, "Mc": 115, "Lv": 116,
         "Ts": 117, "Og": 118, "Uue": 119}
    return z[elem]

def _nband_infer_from_hund(elem, nzeta_max, zval = None):
    """calculate the nband that should be calculated for including all radial functions
    according to Hund's rule that required by nzeta_max, which returned by _make_guess
    function.
    
    Parameters
    ----------
    elem: str
        the element symbol
    nzeta_max: list[int]
        the maximal number of zeta functions for each l
    zval: int, optional
        the valence charge, used for pseudopotential case, default is 0
    
    Returns
    -------
    int: the number of bands to be calculated

    Notes
    -----
    The function works in the following way (take Si as example):
    ```
    Si: 1s 2s 2p 3s 3p 4s ...
     |--- z = 14 --->|
               | <-- | zval = 4
               | --- nzeta_max --- >
    ```
    , therefore if zval = 0 (or not set), the number of bands to calculate will first
    include all 3p states (all partially occupied states), then count the number of bands
    starting from 4s. The zval should be read from pseudopotential.
    """
    assert len(nzeta_max) <= 4, "It is not possible to use atomic method initialize \
orbitals with l >= 3 (f orbitals), because g-orbitals are not observed in nature."

    # if zval not set, set it to 0
    zval = zval or 0

    # (n, l)-pairs
    ener_levels = [(1, 0), # 1s
                   (2, 0), (2, 1), # 2s, 2p 
                   (3, 0), (3, 1), # 3s, 3p
                   (4, 0), # 4s
                   (3, 2), (4, 1), (5, 0), # 3d, 4p, 5s
                   (4, 2), (5, 1), (6, 0), # 4d, 5p, 6s
                   (4, 3), (5, 2), (6, 1), (7, 0), # 4f, 5d, 6p, 7s
                   (5, 3), (6, 2), (7, 1), (8, 0)] # 5f, 6d, 7p, 8s -> Z = 120
    
    nband = 0

    zcore = _elem_z(elem) - zval
    z_accumu = 0

    i = 0
    for n, l in ener_levels:
        z_accumu += 2*(2*l + 1)
        i += 1
        if z_accumu >= zcore: # reach the core, if zval = 0, means all electrons are frozen
            break

    # if z_accumu != zcore, it means there are partially occupied states, first fill them
    nband += (z_accumu - zcore)//2

    # starting from the next level, fill the rest
    nzeta_rest = nzeta_max.copy()
    for n, l in ener_levels[i:]:
        nzeta_rest[l] -= 1
        nband += (2*l + 1)
        if all([nzeta <= 0 for nzeta in nzeta_rest]):
            break

    return max(nband, zval // 2)

import unittest
class TestGuess(unittest.TestCase):
    def test_elem_z(self):
        self.assertEqual(_elem_z("H"), 1)
        self.assertEqual(_elem_z("He"), 2)
        self.assertEqual(_elem_z("Li"), 3)
        self.assertEqual(_elem_z("C"), 6)
        self.assertEqual(_elem_z("Uue"), 119)
    
    def test_nband_infer_from_hund(self):

        # Si: 1s2 2s2 2p6 3s2 3p2

        # if zval is set to 0, will first fill 3p (2 bands), then start from 4s
        # require another [1, 0, 0] zeta functions, 2 + 1 = 3
        self.assertEqual(_nband_infer_from_hund("Si", [1, 0, 0]), 3)
        # require another [1, 1, 0] zeta functions, but 4p is after 3d, so 3d is also included
        # (2 +) 1 + 5 + 3 = 11
        self.assertEqual(_nband_infer_from_hund("Si", [1, 1, 0]), 11)
        # require another [2, 2, 1] zeta functions, so 4s, 5s and 4p, 5p and 3d are wanted
        # the energy levels arranges like 4s, 3d, 4p, 5s, 4d, 5p, ..., so 
        # (2 +) 1 + 5 + 3 + 1 + 5 + 3 = 20
        self.assertEqual(_nband_infer_from_hund("Si", [2, 2, 1]), 20)

        # if zval is set to 4, which means the 3s2 3p2 are valence electrons, 
        # require another [1, 0, 0], the 3s is enough
        self.assertEqual(_nband_infer_from_hund("Si", [1, 0, 0], 4), 2)
        # require another [0, 1, 0], but the 3p is partially occupied, should fill them up
        self.assertEqual(_nband_infer_from_hund("Si", [0, 1, 0], 4), 4)
        # require another [1, 1, 1], 3s, 3p and 3d, but 4s will appear before 3d
        # so 3s, 3p, 4s, 3d, are wanted: 1 + 3 + 1 + 5 = 10
        self.assertEqual(_nband_infer_from_hund("Si", [1, 1, 1], 4), 10)

if __name__ == "__main__":
    unittest.main()