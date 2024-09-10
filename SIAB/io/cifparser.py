"""a easy version of cifparser"""
def read(fcif: str):
    """read the Crystallographic Information File (CIF) and stores information
    as much as possible
    
    Parameter
    ---------
    fcif: str
        the filename of cif
    
    Technical details left here.
    The most common cif file contains several data blocks, each block is seperated
    by one 'loop_', so first split all the content of cif file with it, then for 
    each block, either parse line-by-line and make key-value pairs, or build a
    table if the structure is defined in that way."""
    import os
    out = {}
    if not os.path.exists(fcif):
        raise FileNotFoundError(f"{fcif} not found.")
    with open(fcif, "r") as f:
        content = f.readlines()
    content = [l.replace("\n", " ").strip() for l in content]
    content = " ".join([l for l in content if l and not l.startswith("#") and not l.startswith("data_")])
    blocks = content.split("loop_")
    for block in blocks:
        if len(block) == 0:
            continue
        out.update(_build_cif_block(block))
    return out

def _build_cif_tab(titles: list, value: str):
    """build the table in cif file, with a littlt bit hard-code here: contents in
    the range enclosed by ' will not be split
    
    Parameter
    ---------
    titles: list
        the titles of the table
    value: str
        the value of the table
    
    Return
    ------
    dict
        the table in dictionary format"""
    import numpy as np
    # first make sure there are no adjacent whitespaces in value
    value = " ".join(value.split())
    words = []
    vcache = ""

    if value.startswith(";") and value.endswith(";"):
        value = value.replace(";", "")
        assert len(titles) == 1, f"A paragraph should have only one title, but got {len(titles)}"
        return {titles[0]: [value.strip()]}
    
    for word in value.split():
        if word.startswith("'"):
            vcache += f" {word}"
        elif word.endswith("'"):
            vcache += f" {word}"
            words.append(vcache.strip())
            vcache = ""
        else:
            if vcache:
                vcache += f" {word}"
            else:
                words.append(word.strip())
    # then make sure the length of words is a multiple of len(titles)
    ncols = len(titles)
    nrows = len(words) // ncols
    if nrows * ncols != len(words):
        raise ValueError(f"the length of words is not a multiple of {ncols}")
    return dict(zip(titles, np.array(words).reshape(nrows, ncols).T.tolist()))

def _build_cif_block(raw: str):
    """build the block in cif file, by reading words of raw one-by-one. Once there
    are adjcent definition of keys (those words starts with a underscore), then it
    will be a table, its value will be read and organized as a list, otherwise it
    will be a plain key-value pair"""
    import re
    out, stack = {}, []
    raw = [l.replace("\n", " ").strip() for l in re.split(r"(\s\_\S+)", " " + raw)]
    raw = [l for l in raw if l] # after removal of empty string
    for word in raw:
        if word.startswith("_"):
            stack.append(word)
        else:
            tab = _build_cif_tab(stack, word)
            out.update(tab)
            stack = []
    return out

import unittest
class CifParserTest(unittest.TestCase):

    def test_build_cif_tab(self):
        titles = ["_atom_site_type_symbol", "_atom_site_label", "_atom_site_symmetry_multiplicity", 
                  "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z", "_atom_site_occupancy"]
        raw_table = """
  Se  Se0  1  0.83418970  0.83308821  0.91084093  1
  Se  Se1  1  0.03736952  0.85904742  0.18343119  1
  Se  Se2  1  0.63279384  0.23507614  0.31163980  1
  Se  Se3  1  0.14896496  0.70693955  0.47425805  1
  Se  Se4  1  0.64797437  0.98859060  0.09009684  1
  Se  Se5  1  0.88650771  0.72089911  0.01594700  1
  Se  Se6  1  0.40550744  0.35038179  0.47155027  1
  Se  Se7  1  0.53736952  0.64095258  0.18343119  1
  Se  Se8  1  0.38650771  0.77910089  0.01594700  1
  Se  Se9  1  0.33418970  0.66691179  0.91084093  1
"""
        tab = _build_cif_tab(titles, raw_table)
        self.assertEqual(len(tab), 7)
        self.assertEqual(tab.keys(), set(titles))
        self.assertEqual(tab["_atom_site_type_symbol"], ["Se"]*10)
        self.assertEqual(tab["_atom_site_label"], [f"Se{i}" for i in range(10)])
        self.assertEqual(tab["_atom_site_symmetry_multiplicity"], ["1"]*10)
        self.assertEqual(tab["_atom_site_fract_x"], ["0.83418970", "0.03736952", "0.63279384", "0.14896496", 
            "0.64797437", "0.88650771", "0.40550744", "0.53736952", "0.38650771", "0.33418970"])
        self.assertEqual(tab["_atom_site_fract_y"], ["0.83308821", "0.85904742", "0.23507614", "0.70693955", 
            "0.98859060", "0.72089911", "0.35038179", "0.64095258", "0.77910089", "0.66691179"])
        self.assertEqual(tab["_atom_site_fract_z"], ["0.91084093", "0.18343119", "0.31163980", "0.47425805", 
            "0.09009684", "0.01594700", "0.47155027", "0.18343119", "0.01594700", "0.91084093"])
        self.assertEqual(tab["_atom_site_occupancy"], ["1"]*10)

    def test_build_cif_block(self):
        raw_pair = """
_cell_length_a       3.87933
_cell_length_b       3.87933
_cell_length_c       3.87933
_cell_angle_alpha    60
_cell_angle_beta     60.0001
_cell_angle_gamma    60.0001
"""
        block = _build_cif_block(raw_pair)
        self.assertEqual(len(block), 6)
        self.assertEqual(block["_cell_length_a"][0], "3.87933")
        self.assertEqual(block["_cell_length_b"][0], "3.87933")
        self.assertEqual(block["_cell_length_c"][0], "3.87933")
        self.assertEqual(block["_cell_angle_alpha"][0], "60")
        self.assertEqual(block["_cell_angle_gamma"][0], "60.0001")
        self.assertEqual(block["_cell_angle_beta"][0], "60.0001")

        raw_table = """
 _atom_site_type_symbol
 _atom_site_label
 _atom_site_symmetry_multiplicity
 _atom_site_fract_x
 _atom_site_fract_y
 _atom_site_fract_z
 _atom_site_occupancy
  Se  Se0  1  0.83418970  0.83308821  0.91084093  1
  Se  Se1  1  0.03736952  0.85904742  0.18343119  1
  Se  Se2  1  0.63279384  0.23507614  0.31163980  1
  Se  Se3  1  0.14896496  0.70693955  0.47425805  1
  Se  Se4  1  0.64797437  0.98859060  0.09009684  1
  Se  Se5  1  0.88650771  0.72089911  0.01594700  1
  Se  Se6  1  0.40550744  0.35038179  0.47155027  1
  Se  Se7  1  0.53736952  0.64095258  0.18343119  1
  Se  Se8  1  0.38650771  0.77910089  0.01594700  1
  Se  Se9  1  0.33418970  0.66691179  0.91084093  1
"""
        block = _build_cif_block(raw_table)
        # the table has 7 columns, check each data column dimension
        self.assertEqual(len(block), 7)
        self.assertEqual(len(block["_atom_site_type_symbol"]), 10)
        self.assertEqual(len(block["_atom_site_label"]), 10)
        self.assertEqual(len(block["_atom_site_symmetry_multiplicity"]), 10)
        self.assertEqual(len(block["_atom_site_fract_x"]), 10)
        self.assertEqual(len(block["_atom_site_fract_y"]), 10)
        self.assertEqual(len(block["_atom_site_fract_z"]), 10)
        self.assertEqual(len(block["_atom_site_occupancy"]), 10)
        # the first row
        self.assertEqual(block["_atom_site_type_symbol"][0], "Se")
        self.assertEqual(block["_atom_site_label"][0], "Se0")
        self.assertEqual(block["_atom_site_symmetry_multiplicity"][0], "1")
        self.assertEqual(block["_atom_site_fract_x"][0], "0.83418970")
        self.assertEqual(block["_atom_site_fract_y"][0], "0.83308821")
        self.assertEqual(block["_atom_site_fract_z"][0], "0.91084093")
        self.assertEqual(block["_atom_site_occupancy"][0], "1")
        # the last row
        self.assertEqual(block["_atom_site_type_symbol"][-1], "Se")
        self.assertEqual(block["_atom_site_label"][-1], "Se9")
        self.assertEqual(block["_atom_site_symmetry_multiplicity"][-1], "1")
        self.assertEqual(block["_atom_site_fract_x"][-1], "0.33418970")
        self.assertEqual(block["_atom_site_fract_y"][-1], "0.66691179")
        self.assertEqual(block["_atom_site_fract_z"][-1], "0.91084093")
    
    def test_read(self):
        import uuid, os
        mp1018134 = """
# generated using pymatgen
data_Li
_symmetry_space_group_name_H-M   'P 1'
_cell_length_a   7.66540754
_cell_length_b   7.66540754
_cell_length_c   7.66540754
_cell_angle_alpha   23.16936004
_cell_angle_beta   23.16936004
_cell_angle_gamma   23.16936004
_symmetry_Int_Tables_number   1
_chemical_formula_structural   Li
_chemical_formula_sum   Li3
_cell_volume   61.20564072
_cell_formula_units_Z   3
loop_
 _symmetry_equiv_pos_site_id
 _symmetry_equiv_pos_as_xyz
  1  'x, y, z'
loop_
 _atom_site_type_symbol
 _atom_site_label
 _atom_site_symmetry_multiplicity
 _atom_site_fract_x
 _atom_site_fract_y
 _atom_site_fract_z
 _atom_site_occupancy
  Li  Li0  1  0.77702792  0.77702792  0.77702792  1
  Li  Li1  1  0.22297208  0.22297208  0.22297208  1
  Li  Li2  1  0.00000000  0.00000000  0.00000000  1
"""
        fcif = f"{str(uuid.uuid4())}.cif"
        with open(fcif, "w") as f:
            f.write(mp1018134)
        cif = read(fcif)
        os.remove(fcif)
        self.assertEqual(len(cif), 21) # in total there are 21 key-value pairs
        self.assertEqual(cif["_symmetry_space_group_name_H-M"][0], "'P 1'")
        self.assertEqual(cif["_cell_length_a"][0], "7.66540754")
        self.assertEqual(cif["_cell_length_b"][0], "7.66540754")
        self.assertEqual(cif["_cell_length_c"][0], "7.66540754")
        self.assertEqual(cif["_cell_angle_alpha"][0], "23.16936004")
        self.assertEqual(cif["_cell_angle_beta"][0], "23.16936004")
        self.assertEqual(cif["_cell_angle_gamma"][0], "23.16936004")
        self.assertEqual(cif["_symmetry_Int_Tables_number"][0], "1")
        self.assertEqual(cif["_chemical_formula_structural"][0], "Li")
        self.assertEqual(cif["_chemical_formula_sum"][0], "Li3")
        self.assertEqual(cif["_cell_volume"][0], "61.20564072")
        self.assertEqual(cif["_cell_formula_units_Z"][0], "3")
        self.assertEqual(len(cif["_symmetry_equiv_pos_site_id"]), 1)
        self.assertEqual(len(cif["_symmetry_equiv_pos_as_xyz"]), 1)
        self.assertEqual(len(cif["_atom_site_type_symbol"]), 3)
        self.assertEqual(len(cif["_atom_site_label"]), 3)
        self.assertEqual(len(cif["_atom_site_symmetry_multiplicity"]), 3)
        self.assertEqual(len(cif["_atom_site_fract_x"]), 3)
        self.assertEqual(len(cif["_atom_site_fract_y"]), 3)
        self.assertEqual(len(cif["_atom_site_fract_z"]), 3)
        # check the content of the table
        self.assertEqual(cif["_atom_site_type_symbol"], ["Li"]*3)
        self.assertEqual(cif["_atom_site_label"], ["Li0", "Li1", "Li2"])
        self.assertEqual(cif["_atom_site_symmetry_multiplicity"], ["1"]*3)
        self.assertEqual(cif["_atom_site_fract_x"], ["0.77702792", "0.22297208", "0.00000000"])
        self.assertEqual(cif["_atom_site_fract_y"], ["0.77702792", "0.22297208", "0.00000000"])
        self.assertEqual(cif["_atom_site_fract_z"], ["0.77702792", "0.22297208", "0.00000000"])
        self.assertEqual(cif["_atom_site_occupancy"], ["1"]*3)

        cod1527966 = """
#------------------------------------------------------------------------------
#$Date: 2018-09-27 07:13:35 +0300 (Thu, 27 Sep 2018) $
#$Revision: 211196 $
#$URL: file:///home/coder/svn-repositories/cod/cif/1/52/79/1527966.cif $
#------------------------------------------------------------------------------
#
# This file is available in the Crystallography Open Database (COD),
# http://www.crystallography.net/
#
# All data on this site have been placed in the public domain by the
# contributors.
#
data_1527966
loop_
_publ_author_name
'Radchenko, V.M.'
'Shushakov, V.D.'
'Seleznev, A.G.'
'Ryabinin, M.A.'
'Nikolaev, V.M.'
'Lebedeva, L.S.'
'Vasil'ev, V.Ya.'
_publ_section_title
;
 Pt5 An (AN= Am, Cm, Bk, Cf) intermetallic compounds
;
_journal_name_full               'Journal of the Less-Common Metals'
_journal_page_first              147
_journal_page_last               153
_journal_volume                  157
_journal_year                    1990
_chemical_formula_sum            'Cf Pt5'
_space_group_IT_number           191
_symmetry_space_group_name_Hall  '-P 6 2'
_symmetry_space_group_name_H-M   'P 6/m m m'
_cell_angle_alpha                90
_cell_angle_beta                 90
_cell_angle_gamma                120
_cell_formula_units_Z            1
_cell_length_a                   5.264
_cell_length_b                   5.264
_cell_length_c                   4.417
_cell_volume                     105.996
_citation_journal_id_ASTM        JCOMAH
_cod_data_source_file            Radchenko_JCOMAH_1990_329.cif
_cod_data_source_block           Cf1Pt5
_cod_original_cell_volume        105.9961
_cod_original_formula_sum        'Cf1 Pt5'
_cod_database_code               1527966
loop_
_symmetry_equiv_pos_as_xyz
x,y,z
x-y,x,z
-y,x-y,z
-x,-y,z
-x+y,-x,z
y,-x+y,z
-y,-x,-z
x-y,-y,-z
x,x-y,-z
y,x,-z
-x+y,y,-z
-x,-x+y,-z
-x,-y,-z
-x+y,-x,-z
y,-x+y,-z
x,y,-z
x-y,x,-z
-y,x-y,-z
y,x,z
-x+y,y,z
-x,-x+y,z
-y,-x,z
x-y,-y,z
x,x-y,z
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
_atom_site_U_iso_or_equiv
Cf1 Cf 0 0 0 1 0.0
Pt1 Pt 0.3333 0.6667 0 1 0.0
Pt2 Pt 0.5 0 0.5 1 0.0
"""
        fcif = f"{str(uuid.uuid4())}.cif"
        with open(fcif, "w") as f:
            f.write(cod1527966)
        cif = read(fcif)
        os.remove(fcif)
        self.assertEqual(cif["_publ_author_name"], ["'Radchenko, V.M.'", "'Shushakov, V.D.'", 
        "'Seleznev, A.G.'", "'Ryabinin, M.A.'", "'Nikolaev, V.M.'", "'Lebedeva, L.S.'", "'Vasil'ev, V.Ya.'"])
        self.assertEqual(cif["_symmetry_equiv_pos_as_xyz"], 
        ['x,y,z', 'x-y,x,z', '-y,x-y,z', '-x,-y,z', '-x+y,-x,z', 'y,-x+y,z', '-y,-x,-z', 
         'x-y,-y,-z', 'x,x-y,-z', 'y,x,-z', '-x+y,y,-z', '-x,-x+y,-z', '-x,-y,-z', 
         '-x+y,-x,-z', 'y,-x+y,-z', 'x,y,-z', 'x-y,x,-z', '-y,x-y,-z', 'y,x,z', 
         '-x+y,y,z', '-x,-x+y,z', '-y,-x,z', 'x-y,-y,z', 'x,x-y,z'])
        self.assertEqual(cif["_journal_name_full"], ["'Journal of the Less-Common Metals'"])
        self.assertEqual(cif["_publ_section_title"], ["Pt5 An (AN= Am, Cm, Bk, Cf) intermetallic compounds"])

        cod4031557 = """
#------------------------------------------------------------------------------
#$Date: 2018-09-27 07:13:35 +0300 (Thu, 27 Sep 2018) $
#$Revision: 211196 $
#$URL: file:///home/coder/svn-repositories/cod/cif/4/03/15/4031557.cif $
#------------------------------------------------------------------------------
#
# This file is available in the Crystallography Open Database (COD),
# http://www.crystallography.net/
#
# All data on this site have been placed in the public domain by the
# contributors.
#
data_4031557
loop_
_publ_author_name
'Baybarz, R.D.'
'Fahey, J.A.'
'Haire, R.G.'
_publ_section_title
;
 The preparation, crystal structures and some properties of californium
 oxysulfate and oxysulfide
;
_journal_name_full               'Journal of Inorganic and Nuclear Chemistry'
_journal_page_first              2023
_journal_page_last               2027
_journal_volume                  36
_journal_year                    1974
_chemical_formula_sum            'Cf2 O2 S'
_space_group_IT_number           164
_symmetry_space_group_name_Hall  '-P 3 2"'
_symmetry_space_group_name_H-M   'P -3 m 1'
_cell_angle_alpha                90
_cell_angle_beta                 90
_cell_angle_gamma                120
_cell_formula_units_Z            1
_cell_length_a                   3.844
_cell_length_b                   3.844
_cell_length_c                   6.656
_cell_volume                     85.175
_citation_journal_id_ASTM        JINCAO
_cod_data_source_file            Baybarz_JINCAO_1974_1804.cif
_cod_data_source_block           Cf2O2S1
_cod_original_cell_volume        85.17472
_cod_original_formula_sum        'Cf2 O2 S1'
_cod_database_code               4031557
loop_
_symmetry_equiv_pos_as_xyz
x,y,z
-y,x-y,z
-x+y,-x,z
y,x,-z
-x,-x+y,-z
x-y,-y,-z
-x,-y,-z
y,-x+y,-z
x-y,x,-z
-y,-x,z
x,x-y,z
-x+y,y,z
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
_atom_site_U_iso_or_equiv
Cf1 Cf+3 0.3333 0.6667 0.28 1 0.0
O1 O-2 0.3333 0.6667 0.64 1 0.0
S1 S-2 0 0 0 1 0.0
"""
        fcif = f"{str(uuid.uuid4())}.cif"
        with open(fcif, "w") as f:
            f.write(cod4031557)
        cif = read(fcif)
        os.remove(fcif)
        self.assertEqual(cif["_publ_author_name"], ["'Baybarz, R.D.'", "'Fahey, J.A.'", "'Haire, R.G.'"])
        self.assertEqual(cif['_publ_section_title'], ["The preparation, crystal structures and some properties of californium oxysulfate and oxysulfide"])
        self.assertEqual(cif["_cod_data_source_file"], ["Baybarz_JINCAO_1974_1804.cif"])

if __name__ == "__main__":
    unittest.main()