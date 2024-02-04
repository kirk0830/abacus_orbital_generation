"""a simple module for parsing main information from pseudopotential"""
"""With Python-xml parser, provide general parser for UPF format pseudopotential"""

def preprocess(fname: str):
    """ADC pseudopotential has & symbol at the beginning of line, which is not allowed in xml, replace & with &amp;"""
    with open(fname, "r") as f:
        lines = f.readlines()
    """GBRV pseudopotential does not startswith <UPF version="2.0.1">, but <PP_INFO>, 
    add <UPF version="2.0.1"> to the beginning of the file and </UPF> to the end of the file"""
    if not lines[0].startswith("<UPF version="):
        lines.insert(0, "<UPF version=\"2.0.1\">\n")
        lines.append("</UPF>")

    with open(fname, "w") as f:
        for line in lines:
            """if line starts with &, replace & with &amp;, 
            but if already &amp;, do not replace"""
            if line.strip().startswith("&") and not line.strip().startswith("&amp;"):
                line = line.replace("&", "&amp;")
            
            f.write(line)

            if line.strip() == "</UPF>":
                break

import re
def is_numeric_data(data):
    """judge if the data line is full of numbers (including scientific notation) separated by spaces, tabs or newlines"""
    if re.match(r"^\s*[+-]?\d+(\.\d*)?([eE][+-]?\d+)?(\s+[+-]?\d+(\.\d*)?([eE][+-]?\d+)?)*\s*$", data):
        return True
    else:
        return False

def decompose_data(data):
    """to decompose all numbers in one line, but need to judge whether int or float"""
    if re.match(r"^\s*[+-]?\d+(\.\d*)([eE][+-]?\d+)?(\s+[+-]?\d+(\.\d*)([eE][+-]?\d+)?)*\s*$", data):
        return [float(x) for x in data.split()] if " " in data else float(data.strip())
    elif re.match(r"^\s*[+-]?\d+(\s+[+-]?\d+)*\s*$", data):
        return [int(x) for x in data.split()] if " " in data else int(data.strip())
    else:
        raise ValueError("data is not numeric")

def postprocess(parsed: dict):

    for section in parsed:
        """first the data"""
        if parsed[section]["data"] is not None:
            parsed[section]["data"] = parsed[section]["data"].strip()
            if is_numeric_data(parsed[section]["data"]):
                parsed[section]["data"] = decompose_data(parsed[section]["data"])
        """then the attributes"""
        if parsed[section]["attrib"] is not None:
            for attrib in parsed[section]["attrib"]:
                parsed[section]["attrib"][attrib] = parsed[section]["attrib"][attrib].strip()
                if is_numeric_data(parsed[section]["attrib"][attrib]):
                    parsed[section]["attrib"][attrib] = decompose_data(parsed[section]["attrib"][attrib])
                elif parsed[section]["attrib"][attrib] == "T":
                    parsed[section]["attrib"][attrib] = True
                elif parsed[section]["attrib"][attrib] == "F":
                    parsed[section]["attrib"][attrib] = False
    return parsed

import json
import xml.etree.ElementTree as ET

def iter_tree(root: ET.Element):
    """iterate through the tree, return a dictionary flattened from the tree"""
    return {child.tag: {"attrib": child.attrib, "data": child.text} for child in list(root.iter())}

def parse(fname: str):
    preprocess(fname)
    tree = ET.ElementTree(file=fname)
    root = tree.getroot()
    parsed = iter_tree(root)
    parsed = postprocess(parsed)
    return parsed

def parse_to_json(fname: str):
    parsed = parse(fname)
    with open(fname[:-4]+'.json', 'w') as f:
        json.dump(parsed, f, indent=4)

def determine_type(parsed: dict):
    """pseudopotentials can be generated by not only one codes, to extract information 
    from PP_INFO, PP_INPUTFILE, need to know the exact way how information is organized
    """

    """ONCVPSP
    ONCVPSP is the format of pseudopotential most seen in norm-conserving pseudopotential,
    such as SG15, PD (developed by pwmat team?) and DOJO"""
    if "ONCVPSP" in parsed["PP_INFO"]["data"]:
        return "ONCVPSP"
    if "ONCVPSP" in parsed["PP_HEADER"]["attrib"]["generated"]:
        return "ONCVPSP"

    """ADC
    ADC is the format of pseudopotential collected in pslibrary, including
    pslnc, rrkjus and kjpaw, most collected in QE website the pptable"""
    if "ADC" in parsed["PP_INFO"]["data"]:
        return "ADC"
    if "ADC" in parsed["PP_HEADER"]["attrib"]["generated"]:
        return "ADC"
    if "ADC" in parsed["PP_HEADER"]["attrib"]["author"]:
        return "ADC"
    if "Generated using \"atomic\" code by A. Dal Corso" in parsed["PP_INFO"]["data"]:
        return "ADC"
    if "Generated using \"atomic\" code by A. Dal Corso" in parsed["PP_HEADER"]["attrib"]["generated"]:
        return "ADC"
    if "Generated using \"atomic\" code by A. Dal Corso" in parsed["PP_HEADER"]["attrib"]["author"]:
        return "ADC"
    
    """GTH
    this is the kind developed by CP2K developers, Goedecker, Hartwigsen, Hutter and Teter
    et al. However, this kind of pseudopotential has non-diagonal element in DIJ matrices,
    which is not supported by ABACUS yet."""
    if "Goedecker/Hartwigsen/Hutter/Teter" in parsed["PP_HEADER"]["attrib"]["author"]:
        return "GTH"
        raise NotImplementedError("GTH pseudopotential is not supported by ABACUS yet because of non-diagonal DIJ matrices")
    
    """GBRV
    It is one of the most efficient pseudopotential presently, ABACUS pw supports this kind
    of pseudopotential, ABACUS lcao not yet.
    """
    if "Generated using Vanderbilt code" in parsed["PP_INFO"]["data"]:
        return "GBRV"

    """ATOMPAW
    atompaw looks like ADC but not quite the same in occupation information
    Comparatively the uni_marburg is actually more similar to ADC"""
    if "ATOMPAW" in parsed["PP_INFO"]["data"]:
        return "ATOMPAW"
    if "ATOMPAW" in parsed["PP_HEADER"]["attrib"]["generated"]:
        return "ATOMPAW"
    
    raise ValueError("Pseudopotential type not recognized")

def valelec_config(fname: str):
    """extract valence electron configuration from pseudopotential file
    return a list of lists, 
    [
        ["1s", "2s"] # for s
        ["2p"] # for p
        ["3d"] # for d
        ["4f"] # for f
    ]"""
    parsed = parse(fname)
    pseudo_type = determine_type(parsed)
    if pseudo_type == "ONCVPSP":
        return ONCV_parser(parsed)
    elif pseudo_type == "ADC":
        return ADC_parser(parsed)
    elif pseudo_type == "GBRV":
        return GBRV_parser(parsed)
    elif pseudo_type == "ATOMPAW":
        return ATOMPAW_parser(parsed)
    elif pseudo_type == "GTH":
        raise ValueError("GTH pseudopotential does not have valence electron configuration information")
    else:
        raise ValueError("Pseudopotential type not recognized")

def ONCV_parser(parsed: dict) -> list:

    result = {}

    zval = parsed["PP_HEADER"]["attrib"]["z_valence"]
    
    content = parsed["PP_INPUTFILE"]["data"]
    lines = [line.strip() for line in content.split("\n")]

    reference_config = []
    read_valence_config = False
    for line in lines:
        if line.startswith("#   n    l    f        energy (Ha)"):
            read_valence_config = True
            continue

        if read_valence_config:
            if line.startswith("#"):
                break
            else:
                if len(line.split()) >= 3:
                    reference_config.append(line)
    # reversely read the reference_config
    sequence = ["S", "P", "D", "F", "G", "H", "I", "K", "L", "M", "N"]
    reference_config.reverse()
    for line in reference_config:
        if zval <= 0:
            break
        else:
            words = line.split()
            index = int(words[1])
            symbol = sequence[index]
            if symbol not in result:
                result[symbol] = []
            result[symbol].append(words[0]+symbol)
            zval -= float(words[2])
    # then convert to list
    result_list = []
    for isym, symbol in enumerate(sequence):
        if symbol in result:
            result_list.append(result[symbol])
        else:
            if isym >= len(result):
                break
            else:
                result_list.append([])
    return result_list

def GBRV_parser(parsed: dict) -> list:

    contents = parsed["PP_HEADER"]["data"]
    lines = [line.strip() for line in contents.split("\n")]

    sequence = ["S", "P", "D", "F", "G", "H", "I", "K", "L", "M", "N"]

    result = {}
    read_valence_config = False
    for line in lines:
        if line.startswith("Wavefunctions"):
            read_valence_config = True
            continue
        if read_valence_config:
            if not line[0].isdigit():
                break
            else:
                words = line.split()
                symbol = words[0][-1]
                if symbol not in result:
                    result[symbol] = []
                result[symbol].append(words[0])
    
    # then convert to list
    result_list = []
    for isym, symbol in enumerate(sequence):
        if symbol in result:
            result_list.append(result[symbol])
        else:
            if isym >= len(result):
                break
            else:
                result_list.append([])
    return result_list

def ATOMPAW_parser(parsed: dict) -> list:

    result = {}

    contents = parsed["PP_INFO"]["data"]
    lines = [line.strip() for line in contents.split("\n")]
    
    sequence = ["S", "P", "D", "F", "G", "H", "I", "K", "L", "M", "N"]

    nmax = []
    norb = []
    iorb = 0
    read_valence_config = False
    for il, line in enumerate(lines):
        if il == 3: # this is the line specifying how many orbitals considered for s, p, d, f and g
            nmax = [int(x) for x in line.split()]
            norb = [(nmax[i] - i) if nmax[i] > i else 0 for i in range(len(nmax))]
            continue
        if line == "0 0 0":
            read_valence_config = True
            continue
        if read_valence_config:
            if not line[0].isalpha():
                break
            else:
                if line == "c": # core electron
                    pass
                elif line == "v": # valence electron
                    present_l = 0
                    present_n = 1
                    while iorb >= sum(norb[:present_l]):
                        present_l += 1
                    present_n += iorb - sum(norb[:present_l - 1]) + present_l - 1
                    symbol = sequence[present_l - 1]
                    if symbol not in result:
                        result[symbol] = []
                    result[symbol].append(str(present_n)+symbol)
                else:
                    raise ValueError("Unknown line in valelec_config: {}".format(line))
                iorb += 1
    # then convert to list
    result_list = []
    for isym, symbol in enumerate(sequence):
        if symbol in result:
            result_list.append(result[symbol])
        else:
            if isym >= len(result):
                break
            else:
                result_list.append([])
    return result_list

def ADC_parser(parsed: dict) -> list:

    result = {}

    content = parsed["PP_INFO"]["data"]
    lines = content.split("\n")

    read_valence_config = False
    for line in lines:
        line = line.strip()
        if line.startswith("Valence configuration:"):
            read_valence_config = True
            continue

        if line.startswith("Generation configuration:"):
            break

        if read_valence_config and line[0].isdigit():
            words = line.split()
            if len(words) == 7:
                symbol = words[0][-1]
                if symbol not in result:
                    result[symbol] = []
                result[symbol].append(words[0])
    # then convert to list
    sequence = ["S", "P", "D", "F", "G", "H", "I", "K", "L", "M", "N"]
    result_list = []

    for isym, symbol in enumerate(sequence):
        if symbol in result:
            result_list.append(result[symbol])
        else:
            if isym > len(result):
                break
            else:
                result_list.append([])
    return result_list

def zeta_notation_toorbitalconfig(zeta_notation: str, minimal_basis: list = None):

    pattern = r"([SDTQPH]Z)([SDTQ5-9]?P)?"
    symbols = ["s", "p", "d", "f", "g", "h", "i", "k", "l", "m", "n", "o"]
    multiplier = {"S": 1, "D": 2, "T": 3, "Q": 4, "5": 5, "6": 6, "7": 7, "8": 8, "9": 9}
    _match = re.match(pattern, zeta_notation)
    if _match is None:
        raise ValueError("zeta_notation is not valid")
    nzeta = multiplier[_match.group(1)[0]]
    basis = [nzeta*i for i in minimal_basis]
    result = ""
    for i in range(len(minimal_basis)):
        if basis[i] != 0:
            result += str(basis[i]) + symbols[i]
    if _match.group(2) is not None:
        if len(_match.group(2)) > 1:
            result += str(multiplier[_match.group(2)[0]]) + symbols[len(minimal_basis)]
        else:
            result += "1" + symbols[len(minimal_basis)]
    return result
