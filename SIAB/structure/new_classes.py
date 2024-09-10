from typing import Optional


class AtomSpeciesGenerator:
    """A physical element corresponds to one AtomSpeciesGenerator instance. An iteration
    on it will generate a series of AtomSpecies instances.
    
    During one orbital generation run, the pseudopotential is always fixed, so the
    pseudopot file can be included in the param to instantiate this class"""


    elem: str
    fpseudo: str

    def __init__(self,
                 elem: str,
                 fpseudo: str) -> None:
        self.elem = elem
        self.fpseudo = fpseudo

    def forb(elem: str, ecut: float, rcut: float, nzeta: list):
        """Generate the filename of the numerical atomic orbital file.
        
        Parameter
        ---------
        elem: str
            The element symbol.
        ecut: float
            The energy cutoff.
        nzeta: list
            The list of zeta values.
        """
        syms = "spdfghiklmnortuvwxyz"
        suffix = "".join([f"{nzeta[j]}{syms[j]}" for j in range(len(nzeta))])
        ecut = f"{ecut}Ry"
        rcut = f"{rcut}au"
        return "_".join([elem, rcut, ecut, suffix])

    def __call__(self,
                 ecut: float,
                 lmax: int,
                 orbital_dir = None,
                 rcut = None):
        """Generate a series of AtomSpecies instances."""
        import SIAB.data.interface as database
        name = self.elem
        fullname = database.PERIODIC_TABLE_TOFULLNAME[self.elem]
        index = database.PERIODIC_TABLE_TOINDEX[self.elem]
        rcovalent = database.RCOVALENT[self.elem]
        mass = 1.0000 # hard code for now, cannot be used on MD simulation
        magmom = 0.0
        init_set = dict(zip(["name", "fullname", "symbol", "index", "rcovalent", 
                             "mass", "magmom", "pp", "nao"],
                            [name, fullname, self.elem, index, rcovalent, 
                             mass, magmom, self.fpseudo, None]))
        if rcut is None:
            yield AtomSpecies(**init_set)
        from SIAB.spillage.api import _coef_gen, _save_orb
        import os
        for rc in rcut:
            coef_ = _coef_gen(rc, ecut, lmax, "eye") # jY
            nzeta = [len(coef_[l]) for l in range(lmax+1)]
            forb = AtomSpeciesGenerator.forb(name, ecut, rc, nzeta)
            # first check if this file exists, if so, will not generate it again
            if not os.path.exists(os.path.join(orbital_dir, forb)):
                # generate the numerical atomic orbital file
                _save_orb([[coef_]], self.elem, ecut, rc, nzeta)

class AtomSpecies:
    """A element with specific pseudopotential and numerical atomic orbital file corresponds to
    one AtomSpecies instance"""
    label, fullname, symbol, index, rcovalent, mass, magmom = \
        None, None, None, None, None, None, None
    pp, nao = None, None
    ecutwfc = None
    
    def __init__(self, **kwargs) -> None:
        self.label = kwargs.get('label', None)
        self.fullname = kwargs.get('fullname', None)
        self.symbol = kwargs.get('symbol', None)
        self.index = kwargs.get('index', None)
        self.rcovalent = kwargs.get('rcovalent', None)
        self.mass = kwargs.get('mass', None)
        self.magmom = kwargs.get('magmom', None)
        self.pp = kwargs.get('pp', None)
        self.nao = kwargs.get('nao', None)
        self.ecutwfc = kwargs.get('ecutwfc', None)

    def __str__(self) -> str:
        """overload the print function to print the AtomSpecies information"""
        return (f"AtomSpecies: name: {self.label}, fullname: {self.fullname}, "
                f"symbol: {self.symbol}, index: {self.index}, covalent radius: {self.rcovalent}, "
                f"mass: {self.mass}, magmom: {self.magmom}, pseudopotential: {self.pp}, "
                f"Recommended ecutwfc: {self.ecutwfc} Ry, numerical atomic orbital: {self.nao}")

    def as_dict(self) -> dict:
        return {"name": self.label, "fullname": self.fullname, "symbol": self.symbol,
                "index": self.index, "rcovalent": self.rcovalent, "mass": self.mass,
                "magmom": self.magmom, "pp": self.pp, "ecutwfc": self.ecutwfc, "nao": self.nao}
