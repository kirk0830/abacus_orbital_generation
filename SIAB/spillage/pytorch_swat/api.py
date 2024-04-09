import SIAB.interface.old_version as siov
import SIAB.spillage.pytorch_swat.main as sspsm
import SIAB.spillage.orbscreen as sso
def run(params: dict = None, ilevel: int = 0, nlevel: int = 3):
    """Run the spillage calculation
    
    Args:
        params (dict): parameters for the spillage calculation
        ilevel (int): the current level of the calculation
        nlevel (int): the total number of levels of the calculation
    
    Returns:
        tuple: a tuple containing the path of the orbital file and the screen values
    """
    
    """convert-back the information organized in the way that is acceptable
    for the original version of SIAB to the following format:
    ```python
    return {
        "element": element,
        "ecutwfc": ecutwfc,
        "rcut": rcut,
        "zeta_notation": zeta_notation,
    }
    ```
    """
    chkpt = siov.unpack(orb_gen=params)
    folder = siov.folder(unpacked_orb=chkpt)
    if is_duplicate(folder):
        return
    
    if params is None:
        sspsm.main()
    else:
        sspsm.main(params)
    
    refresh = True if ilevel == nlevel-1 else False
    # refresh is for controlling whether leave the ORBITAL_RESULTS.txt in the original folder
    # if so, orbitals will be optimized based on the previous results, say hierarchical optimization
    # but for the last level, we should remove the ORBITAL_RESULTS.txt.
    files = checkpoint(src="./", dst=folder, this_point=chkpt, refresh=refresh)

    # analysis the newly generated orbitals' kinetic energies of all components
    screen_vals = sso.screen(fnao=files[4], item="T")

    return files[4], screen_vals

import os
import re
def is_duplicate(folder: str):
    """check if the siab calculation is skipped"""
    if not os.path.isdir(folder):
        return False
    orbital_u = r"^(ORBITAL_)([0-9]+)(U\.dat)$"
    files = os.listdir(folder)
    print("Checking files in %s..."%folder)
    if "Spillage.dat" in files:
        print("    Spillage.dat exists")
        if "ORBITAL_RESULTS.txt" in files:
            print("    ORBITAL_RESULTS.txt exists")
            if "ORBITAL_PLOTU.dat" in files:
                print("    ORBITAL_PLOTU.dat exists")
                if "SIAB_INPUT" in files:
                    print("    SIAB_INPUT exists")
                    for file in files:
                        if re.match(orbital_u, file):
                            print("    ORBITAL_*U.dat exists\n=> Restart check pass.")
                            return True
    return False

import SIAB.interface.env as sienv
import SIAB.data.interface as sdi
def checkpoint(src: str,
               dst: str,
               this_point: dict,
               refresh: bool = False,
               env: str = "local",):
    """After optimization of numerical orbitals' coefficients,
        move generated orbitals to the folder named as:
        [element]_gga_[Ecut]Ry_[Rcut]au_[orbital_config]

        ONCE ONE OPTIMIZATION TASK COMPLETES, CALL THIS FUNCTION.
    Design:
        all information should be included in user_settings
        rather than externally defined additionally.
    
    Args:
        user_settings (dict): user settings
        rcut (float): cutoff radius
        orbital_config (str): orbital configuration, e.g. 1s1p
    
    Returns:
        all files path after mv/cp
    """
    # first check if the folder exists, if not, create it
    element = this_point["element"]
    ecutwfc = this_point["ecutwfc"]
    rcut = this_point["rcut"]
    orbital_config = this_point["zeta_notation"]

    files = []
    if not os.path.isdir(dst):
        sienv.op("mkdir", dst, additional_args=["-p"], env=env)
    """backup input file, unlike the original version, we fix it must be named as SIAB_INPUT"""
    sienv.op("cp", "%s/SIAB_INPUT"%src, "%s/SIAB_INPUT"%dst, env=env)
    files.append("%s/SIAB_INPUT"%dst)
    """move spillage.dat"""
    sienv.op("mv", "%s/Spillage.dat"%src, "%s/Spillage.dat"%dst, env=env)
    files.append("%s/Spillage.dat"%dst)
    """move ORBITAL_PLOTU.dat and ORBITAL_RESULTS.txt"""
    sienv.op("mv", "%s/ORBITAL_PLOTU.dat"%src, "%s/ORBITAL_PLOTU.dat"%dst, env=env)
    files.append("%s/ORBITAL_PLOTU.dat"%dst)
    if not refresh:
        sienv.op("cp", "%s/ORBITAL_RESULTS.txt"%src, "%s/ORBITAL_RESULTS.txt"%dst, env=env)
        results_backup(src=src, env=env)
    else:
        sienv.op("mv", "%s/ORBITAL_RESULTS.txt"%src, "%s/ORBITAL_RESULTS.txt"%dst, env=env)
        sienv.op("rm", "%s/Level*.ORBITAL_RESULTS.txt"%src, env=env)
    files.append("%s/ORBITAL_RESULTS.txt"%dst)
    """move ORBITAL_[element]U.dat to [element]_gga_[Ecut]Ry_[Rcut]au.orb"""
    forb = "%s_gga_%sRy_%sau_%s.orb"%(element, str(ecutwfc), str(rcut), orbital_config)
    index = sdi.PERIODIC_TABLE_TOINDEX[element]
    sienv.op("cp", "%s/ORBITAL_%sU.dat"%(src, index), "%s/%s"%(dst, forb), env=env)
    files.append("%s/%s"%(dst, forb))
    print("Orbital file %s generated."%forb)
    """and directly move it to the folder"""
    sienv.op("mv", "%s/ORBITAL_%sU.dat"%(src, index), "%s/ORBITAL_%sU.dat"%(dst, index), env=env)
    files.append("%s/ORBITAL_%sU.dat"%(dst, index))

    return files

def results_backup(src: str = "./", env: str = "local"):

    ilevel = 0
    while True:
        if os.path.isfile("%s/Level%s.ORBITAL_RESULTS.txt"%(src, ilevel)):
            ilevel += 1
        else:
            break
    sienv.op("mv", "%s/ORBITAL_RESULTS.txt"%src, "%s/Level%s.ORBITAL_RESULTS.txt"%(src, ilevel), env=env)
