from SIAB.structure.atom_species_and_cell \
    import AtomSpecies, Cell, AtomSpeciesGeneartor, CellGenerator

class OrbitalMatrixCalculator:
    """orbital matrix calculator is designed as an ABACUS job manager to care for how to execute the series of
    ABACUS jobs to calculate the orbital matrix for orbital generation. 
    For orbital generated from PW calculation, from one PW calculation any number of rcut orbmat can be 
    calculated, however, for the jY basis, the rcut of orbital generated from it, should have the same
    rcut, so one jY calculation can only support generation of orbitals with the same rcut.
    The difference between PW and jY defines two different way to run ABACUS, one starts an SCF calculation,
    then calculate orbmats of rcuts one-by-one. The other starts an SCF calculation, then calculate orbmat
    for only one rcut, then start another SCF calculation for the next rcut:
        
        1. PW:
            |-------> during the runtime of ABACUS PW SCF calculation ------------------>|
            |  SCF  -> cal_orbmat(rcut1) -> cal_orbmat(rcut2) -> cal_orbmat(rcut3) -> ...|
        2. jY:
            |-------> during the runtime of ABACUS jY SCF calculation ------------------>|
            |  SCF  -> cal_orbmat(rcut1) -> SCF -> cal_orbmat(rcut2) -> SCF ->    ...    |
    
    However, actually the PW can also proceed in the same way as jY, but it is not recommended, because
    the low efficiency. We will call the first way as 'one-shot' and the second way as 'for-each'. Because
    the 'one-shot' is not supported by both jY and PW, therefore the default will be 'for-each'.

    The concept of this class? 
    given one structure and a series of rcuts, use the same DFTParamSet to calculate the orbital matrix 
    for each rcut.

    What should be used to initialize/instantiate this class?
    All needed to calculate the orbmat are: DFTParamSet, structure (Cell instance), rcuts (list of float).
    But 
    """

    mode: str = "local-serial" # can be 'local-serial', 'local-parallel', 'remote-serial', 'remote-parallel'

    def __init__(self,
                 asgens: list[AtomSpeciesGeneartor],
                 proto: str,
                 pertkind: str,
                 pertmag: list|str,
                 mode: str = "local-serial"):
        """instantiate the class with qsub to submit the ABACUS job and qstat to check the job status.
        
        Parameters
        ----------
        qsub : str
            the command to submit the ABACUS job.
        qstat : str
            the command to check the job status.
        """
        if mode not in ["local-serial", "local-parallel", "remote-serial", "remote-parallel"]:
            raise ValueError("mode must be one of 'local-serial', 'local-parallel', 'remote-serial', 'remote-parallel'")
        if mode != "local-serial":
            raise NotImplementedError("only 'local-serial' is supported now.")
        self.mode = mode

    def run(self,
            cell,
            dftparam: dict,
            rcuts: list):
        for rcut in rcuts:
            # new a copy of dftparam
            _dftparam = dftparam.copy().update({"bessel_nao_rcut": rcut})
            if self.mode.endswith("serial"):
                pass
            elif self.mode.endswith("parallel"):
                pass

    def _run_serial(cell,
                    dftparam: dict):
        pass