from SIAB.structure.atom_species_and_cell import \
    AtomSpeciesGeneartor, CellGenerator, AtomSpecies, Cell
##############################################
#    prepare based on AtomSpecies and Cell   #
##############################################
def export(paramset, atomset, cell):
    keys = ["INPUT", "STRU", "KPT"]
    # here it is possible to support the converged ecutwfc value auto-set for INPUT.
    # first get the fpp from atomset, then get the max, set to ecutwfc in paramset
    ecutwfc_set = paramset.get("ecutwfc") # if ecut_set is None or "auto", then set it to the max of fpp
    ecutrho_set = paramset.get("ecutrho")
    # get a copy of paramset to let the original unchanged
    paramset = paramset.copy()
    if ecutwfc_set is None or ecutwfc_set == "auto":
        candidates = [as_.ecutwfc for as_ in atomset if as_.ecutwfc is not None]
        ecutwfc = 100 if len(candidates) == 0 else max(candidates)
        print(f"AUTOSET: ecutwfc is autoset to {ecutwfc} Ry")
        paramset["ecutwfc"] = ecutwfc
    if isinstance(ecutrho_set, (int, float)) and ecutrho_set < 0: # negative value indicates the dual
        ecutrho = ecutwfc * abs(ecutrho_set)
        paramset["ecutrho"] = ecutrho
    vals = [write_abacus_input(paramset), write_abacus_stru(atomset, cell), write_abacus_kpt(cell)]
    return dict(zip(keys, vals))

def write_abacus_input(paramset: dict):
    # first precondition such that all parameters are set as string
    paramset = {k: str(v) if not isinstance(v, list) else " ".join([str(_v) for _v in v]) for k, v in paramset.items()}
    result = "INPUT_PARAMETERS\n"
    result += "\n".join([f"{k:<20s} {v:<s}" for k, v in paramset.items()])
    return result

def write_abacus_stru(atomset: list[AtomSpecies], cell: Cell):
    from os.path import basename
    result = "ATOMIC_SPECIES\n"
    # need a map from cell.labels to index of AtomSpecies in atomset list!
    temp_ = [a.symbol for a in atomset]
    uniquelabels_atomspecies_map = [temp_.index(cell.kinds[cell.labels_kinds_map[cell.labels.index(ulbl)]]) for ulbl in dict.fromkeys(cell.labels)]
    result += "\n".join([f"{label:<4s} \
{atomset[uniquelabels_atomspecies_map[il]].mass:<8.4f} {basename(atomset[uniquelabels_atomspecies_map[il]].pp)}" \
            for il, label in enumerate(dict.fromkeys(cell.labels))])
    
    if not all([as_.nao is None for as_ in atomset]):
        result += "\n\nNUMERICAL_ORBITAL\n"
        result += "\n".join([f"{basename(atomset[uniquelabels_atomspecies_map[il]].nao)}" for il in range(len(set(cell.labels)))])

    result += f"\n\nLATTICE_CONSTANT\n{cell.lat0:<20.10f}\n"
    result += "\nLATTICE_VECTORS\n"
    latvec = CellGenerator.abc_angles_to_vec([cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma], True)
    result += "\n".join([f"{latvec[i][0]:<20.10f} {latvec[i][1]:<20.10f} {latvec[i][2]:<20.10f}" for i in range(3)])
    
    coord = "Direct" if cell.periodic else "Cartesian"
    result += f"\n\nATOMIC_POSITIONS\n{coord}\n"
    for label in dict.fromkeys(cell.labels):
        ind = [i for i, l in enumerate(cell.labels) if l == label]
        result += f"{label}\n{cell.magmoms[ind[0]]:<4.2f}\n{len(ind)}\n"
        for i in ind:
            result += f"{cell.coords[i][0]:<20.10f}{cell.coords[i][1]:<20.10f}{cell.coords[i][2]:<20.10f} \
m {cell.mobs[i][0]:<2d}{cell.mobs[i][1]:<2d}{cell.mobs[i][2]:<2d}\n"
    return result

def write_abacus_kpt(cell: Cell):
    """it is not clarified how to perform a two-step calculation like scf-nscf, for example the
    band structure calculation, therefore, the band structure calculation is not supported yet."""
    result = f"K_POINTS\n0\nGamma\n\
{cell.mpmesh_nks[0]} {cell.mpmesh_nks[1]} {cell.mpmesh_nks[2]} 0 0 0"
    return result