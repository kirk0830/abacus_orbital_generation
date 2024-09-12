from SIAB.structure.atom_species_and_cell \
    import AtomSpecies, Cell, AtomSpeciesGeneartor, CellGenerator
import subprocess as sp
"""
The ABACUS numerical atomic orbital generation task requires all calculations
of orbmat ready, while the calculation requires all required jy basis ready.
Here a general structure of this process is given.

A practical example:

Given orb1, orb2, orb3, orb4
orb1 requires dft1, dft2
orb2 requires dft2, dft3
orb3 requires orb2
orb4 requires dft4, orb1

dft1 requires atomspecies1, atomspecies2, atomspecies3
dft2 requires atomspecies2
dft3 requires atomspecies3

the generation of orbs can be represented as:
orb1: [[atomspecies1, atomspecies2, atomspecies3], [dft1, dft2], []]
orb2: [[atomspecies2], [dft2, dft3], []]
orb3: [[atomspecies2], [dft2, dft3], [orb2]]
orb4: [[atomspecies1, atomspecies2, atomspecies3], [dft1, dft2], [orb1]]

we define each orb generation task as a TreeProcess, there are several stages
in it, in which all tasks can be run in parallel (or say they are independent
and equivalent), which is represented as a TaskGroup, and in one TaskGroup,
there are many Tasks defined.

during the running of a TreeProcess, if one TaskGroup is finished, then call
the finalize() function of it, then initialize() the next TaskGroup and call
its run(). If all TaskGroups are finished, then the TreeProcess is finished,
call finalize() and return.

for each task, should have a function to check its status, initialized/
running/finished.
"""

class Task:
    """a base class, should define the run() function"""
    
    status: str = "initialized" # can be new/running/finished

    def __init__(self) -> None:
        pass

    def run() -> None:
        pass

    def check_status(self) -> str:
        return self.status

class AtomSpeciesGeneration(Task):
    """generate an AtomSpecies instance"""
    elem: str
    fpseudo: str
    ecutwfc: float
    orbital_dir: str
    rcut: float
    lmax: int

    def __init__(self,
                 elem: str,
                 fpseudo: str, 
                 ecutwfc, 
                 orbital_dir = None,
                 rcut = None, 
                 lmax = None) -> None:
        """configure the AtomSpeciesGeneration task. It is okay to
        leave orbital_dir, rcut, lmax as None, because they are conditionally
        used in the run() function (for jy basis generation).
        
        Parameters
        ----------
        elem : str
            the element symbol
        fpseudo : str
            the pseudopotential file path
        ecutwfc : float
            the plane wave cutoff energy
        orbital_dir : str
            the directory to store the generated jy basis
        rcut : float
            the cutoff radius of the orbital
        lmax : int
            the maximum angular momentum of jy basis
        """
        self.elem = elem
        self.fpseudo = fpseudo
        self.ecutwfc = ecutwfc
        self.orbital_dir = orbital_dir
        self.rcut = rcut
        self.lmax = lmax

    def run(self):
        """Generate the corresponding AtomSpecies instance"""
        import os
        asgen = AtomSpeciesGeneartor(self.elem, os.path.dirname(self.fpseudo), os.path.basename(self.fpseudo), self.ecutwfc)
        as_ = [a for a in asgen([self.rcut], self.lmax, self.orbital_dir)]
        assert len(as_) == 1 # actually there is only one, because AtomSpecies with different rcut will be different
        self.status = "finished"
        return as_[0]

    def __eq__(self, value: object) -> bool:
        if not isinstance(value, AtomSpeciesGeneration):
            return False
        return self.fpseudo == value.fpseudo and \
            self.ecutwfc == value.ecutwfc and \
            self.rcut == value.rcut and \
            self.lmax == value.lmax

class OrbitalMatrixCalculation(Task):
    """task of orbital matrix calculation, spillage optimzier will read all
    of them needed. There is a little bit tricky task that, for PW, calculation
    of orbital matrix of different rcut can be calculated in-one-shot, that
    is because the PW itself does not have rcut, while for jy basis employed
    LCAO calculation which generates orbital matrix, for one jy SCF run, only
    one rcut can be calculated."""
    # physically...
    dftparam: dict
    atomspecies: list[AtomSpecies]
    cell: Cell
    command: list[str]

    # internally...
    handle: sp.Popen
    def __init__(self, 
                 dftparam: dict, 
                 atomspecies: list[AtomSpecies], 
                 cell: Cell,
                 command: list[str]) -> None:
        self.dftparam = dftparam
        self.atomspecies = atomspecies
        self.cell = cell
        self.command = command

    def run(self):
        """run the ABACUS dft job to calculate orbital matrix. Because this task
        is time-consuming, new strategy will be implemented to let it run in parallel.
        However, it is not cared in this level of class because this task itself can
        only know how to run itself, not how to run in parallel."""
        self.status = "running"
        # integrate the dftparam as INPUT, atomspecies and cell as STRU...
        self.handle = sp.Popen(args=self.command, stdout=sp.PIPE, stderr=sp.PIPE)

    def check_status(self) -> str:
        if not hasattr(self, "handle"):
            self.status = "initialized"
        elif self.handle.poll() is None:
            self.status = "running"
        else:
            self.status = "finished"
        return self.status

    def __eq__(self, value: object) -> bool:
        if not isinstance(value, OrbitalMatrixCalculation):
            return False
        return self.dftparam == value.dftparam and \
            self.atomspecies == value.atomspecies and \
            self.cell == value.cell and \
            self.command == value.command

class SpillageOptimization(Task):
    def __init__(self,
                 optimizer: str,
                 max_steps: int,
                 nzeta: list[int]) -> None:
        pass

class TaskGroup:
    """a group of tasks"""
    tasks: list[Task]
    status: str

    def __init__(self, recipe: list) -> None:
        self.status = "initialized"
        self.tasks = self._build_tasks(recipe)

    def _build_tasks(self, recipe: list) -> list[Task]:
        """build tasks based on the recipe"""
        return []

    def append(self, task: Task) -> None:
        self.tasks.append(task)

    def initialize(self, upstream) -> None:
        """initialize the taskgroup based on the previous taskgroup"""
        pass

    def check_status(self):
        """check the status of the taskgroup"""
        if all([task.check_status() == "finished" for task in self.tasks]):
            self.status = "finished"
        elif any([task.check_status() == "running" for task in self.tasks]):
            self.status = "running"
        else:
            self.status = "initialized"
        return self.status

    def finalize(self) -> None:
        """after complete all tasks defined, run this function"""
        pass

class ASGenStage(TaskGroup):
    """run all AtomSpeciesGeneration tasks required by some OMCal tasks in parallel.
    For each AtomSpecies, the following keys are required:
    Compulsory: elem (str), fpseudo (str), ecutwfc (float), 
    Optional: orbital_dir (str), rcut (float), lmax (int).
    
    Organize them in a dict for each AtomSpecies, and store in a list as `recipe`
    to instantiate the ASGenStage."""

    results: list[AtomSpecies]

    def _build_tasks(self, recipe: list[dict]):
        keys = ["elem", "fpseudo", "ecutwfc", "orbital_dir", "rcut", "lmax"]
        extractor = lambda x: dict(zip(keys, [x.get(k, None) for k in keys]))
        recipe_ = [extractor(r) for r in recipe]
        return [AtomSpeciesGeneration(**r) for r in recipe_]

    def run(self):
        """the run of AtomSpeciesGeneration will be trivial in time-consuming,
        so it is okay to directly run them in sequence, directly implement"""
        self.results = [task.run() for task in self.tasks]
        assert all([task.check_status() == "finished" for task in self.tasks])

    def finalize(self):
        return self.results

class OMCalStage(TaskGroup):
    """run all independent Orbital Matrix Calculation tasks, register to
    the TaskServer"""

    limit: int
    
    def _build_tasks(self, recipe: list):
        """build tasks for the TaskGroup to run and manage, for OMCalStage,
        the recipe is a list of fixed dict like
        {}
        """

    def initialize(self, upstream: ASGenStage):
        pass

    def run(self) -> None:
        """submit a series of orbital matrix calculation tasks. For this
        specific task, usually it is quite time-consuming, so the task"""
        pass

class SOptStage(TaskGroup):
    """prerequisite for spillage"""

    def __init__(self, tasks: list[SpillageOptimization]) -> None:
        super().__init__(tasks)

    def run(self) -> None:
        pass

    def finalize(self) -> None:
        """optimize itself in finalize function"""
        pass

class TreeProcess:
    """a tree-like process"""

    plan: list[TaskGroup]
    status: str
    def __init__(self, groups: list[TaskGroup]) -> None:
        
        pass

    def run(self) -> None:
        pass

import unittest
class TestAtomSpeciesGeneration(unittest.TestCase):
    
    def test_constructor(self):
        asgen_pw = AtomSpeciesGeneration("Si", "/root/pseudo/Si.pbe-n-rrkjus_psl.1.0.0.UPF", 100.0)
        self.assertEqual(asgen_pw.elem, "Si")
        self.assertEqual(asgen_pw.fpseudo, "/root/pseudo/Si.pbe-n-rrkjus_psl.1.0.0.UPF")
        self.assertEqual(asgen_pw.ecutwfc, 100.0)
        self.assertIsNone(asgen_pw.orbital_dir)
        self.assertIsNone(asgen_pw.rcut)
        self.assertIsNone(asgen_pw.lmax)
        self.assertEqual(asgen_pw.check_status(), "initialized")
        asgen_jy = AtomSpeciesGeneration("Si", "/root/pseudo/Si.pbe-n-rrkjus_psl.1.0.0.UPF", 100.0, "/root/orbital", 10.0, 3)
        self.assertEqual(asgen_jy.elem, "Si")
        self.assertEqual(asgen_jy.fpseudo, "/root/pseudo/Si.pbe-n-rrkjus_psl.1.0.0.UPF")
        self.assertEqual(asgen_jy.ecutwfc, 100.0)
        self.assertEqual(asgen_jy.orbital_dir, "/root/orbital")
        self.assertEqual(asgen_jy.rcut, 10.0)
        self.assertEqual(asgen_jy.lmax, 3)
        self.assertEqual(asgen_jy.check_status(), "initialized")

    def test_run(self):
        import os
        asgen_pw = AtomSpeciesGeneration("Si", "/root/pseudo/Si.pbe-n-rrkjus_psl.1.0.0.UPF", 100.0)
        a_pw = asgen_pw.run()
        self.assertEqual(a_pw.as_dict(), 
            {'name': None, 'fullname': 'Silicon', 'symbol': 'Si', 'index': 14, 
             'rcovalent': 1.11, 'mass': 1.0, 'magmom': 0.0, 'pp': '/root/pseudo/Si.pbe-n-rrkjus_psl.1.0.0.UPF', 'ecutwfc': 100.0, 'nao': None})
        asgen_jy = AtomSpeciesGeneration("Si", "/root/pseudo/Si.pbe-n-rrkjus_psl.1.0.0.UPF", 100.0, os.getcwd(), 10.0, 3)
        a_jy = asgen_jy.run()
        forb = os.path.join(os.getcwd(), a_jy.nao)
        fpng = os.path.join(os.getcwd(), a_jy.nao[:-4] + ".png")
        fparam = os.path.join(os.getcwd(), "ORBITAL_RESULTS.txt")
        os.remove(forb)
        os.remove(fpng)
        os.remove(fparam)
        self.assertEqual(a_jy.as_dict(), 
            {'name': None, 'fullname': 'Silicon', 'symbol': 'Si', 'index': 14, 
             'rcovalent': 1.11, 'mass': 1.0, 'magmom': 0.0, 'pp': '/root/pseudo/Si.pbe-n-rrkjus_psl.1.0.0.UPF', 'ecutwfc': 100.0, 
             'nao': f'{forb}'})

class TestASGenStage(unittest.TestCase):

    def test_constructor(self):
        recipe = [
            {"elem": "Si", "fpseudo": "/root/pseudo/Si.pbe-n-rrkjus_psl.1.0.0.UPF", "ecutwfc": 100.0},
            {"elem": "Si", "fpseudo": "/root/pseudo/Si.pbe-n-rrkjus_psl.1.0.0.UPF", "ecutwfc": 100.0, "orbital_dir": "/root/orbital", "rcut": 10.0, "lmax": 3}
        ] # one for pw, one for jy
        stage = ASGenStage(recipe)
        self.assertEqual(len(stage.tasks), 2)
        self.assertIsInstance(stage.tasks[0], AtomSpeciesGeneration)
        self.assertIsInstance(stage.tasks[1], AtomSpeciesGeneration)
    
    def test_run(self):
        import os
        recipe = [
            {"elem": "Si", "fpseudo": "/root/pseudo/Si.pbe-n-rrkjus_psl.1.0.0.UPF", "ecutwfc": 100.0},
            {"elem": "Si", "fpseudo": "/root/pseudo/Si.pbe-n-rrkjus_psl.1.0.0.UPF", "ecutwfc": 100.0, "orbital_dir": os.getcwd(), "rcut": 10.0, "lmax": 3}
        ]
        stage = ASGenStage(recipe)
        stage.run()
        self.assertEqual(len(stage.results), 2)
        self.assertIsInstance(stage.results[0], AtomSpecies)
        self.assertIsInstance(stage.results[1], AtomSpecies)
        # remove the generated files
        forb = stage.results[1].nao
        fpng = forb[:-4] + ".png"
        fparam = os.path.dirname(forb) + "/ORBITAL_RESULTS.txt"
        os.remove(forb)
        os.remove(fpng)
        os.remove(fparam)

        asref = [
            {'name': None, 'fullname': 'Silicon', 'symbol': 'Si', 'index': 14, 
             'rcovalent': 1.11, 'mass': 1.0, 'magmom': 0.0, 
             'pp': '/root/pseudo/Si.pbe-n-rrkjus_psl.1.0.0.UPF', 'ecutwfc': 100.0, 'nao': None},
            {'name': None, 'fullname': 'Silicon', 'symbol': 'Si', 'index': 14, 
             'rcovalent': 1.11, 'mass': 1.0, 'magmom': 0.0, 
             'pp': '/root/pseudo/Si.pbe-n-rrkjus_psl.1.0.0.UPF', 'ecutwfc': 100.0, 
             'nao': '/root/abacus-develop/ABACUS-ORBGEN/Si_gga_10.0au_100.0Ry_31s31p30d30f.orb'}
        ]
        for ia, a in enumerate(stage.results):
            self.assertEqual(a.as_dict(), asref[ia])
        
class TestOrbitalMatrixCalculation(unittest.TestCase):
    """the test on orbitam matrix calculation is actually test the use of subprocess,
    let a task can instantly run and check its status"""
    def test_constructor(self):
        asgen_pw = AtomSpeciesGeneration("Si", "/root/pseudo/Si.pbe-n-rrkjus_psl.1.0.0.UPF", 100.0)
        a_pw = asgen_pw.run()
        cell = Cell()
        dftparam = {"ecutwfc": 100.0}
        abacus = OrbitalMatrixCalculation(dftparam, [a_pw], cell, "which abacus")
        self.assertEqual(abacus.dftparam, dftparam)
        self.assertEqual(abacus.atomspecies, [a_pw])
        self.assertEqual(abacus.cell, cell)
        self.assertEqual(abacus.command, "which abacus")
        self.assertEqual(abacus.check_status(), "initialized")

    def test_run(self):
        import os, uuid, time
        asgen_pw = AtomSpeciesGeneration("Si", "/root/pseudo/Si.pbe-n-rrkjus_psl.1.0.0.UPF", 100.0)
        a_pw = asgen_pw.run()
        cell = Cell()
        dftparam = {"ecutwfc": 100.0}
        script = '''
import time
if __name__ == "__main__":
    print("hello world")
    time.sleep(5)
    print("goodbye world")
'''
        fscr = os.path.join(os.getcwd(), f"{uuid.uuid4()}.py")
        with open(fscr, "w") as f:
            f.write(script)
        abacus = OrbitalMatrixCalculation(dftparam, [a_pw], cell, ["python3", fscr])
        abacus.run()
        
        status = []
        for _ in range(10):
            status.append(abacus.check_status())
            time.sleep(1)
        self.assertEqual(status,
            ['running', 'running', 'running', 'running', 'running', 'running', 
             'finished', 'finished', 'finished', 'finished'])
        
        os.remove(fscr)



if __name__ == "__main__":
    unittest.main()

    asgen = AtomSpeciesGeneration("pseudo", 10.0, 10.0, 10)
    abacus = OrbitalMatrixCalculation({"ecutwfc": 10.0}, [asgen], Cell(10.0), "abacus")
    spillage = SpillageOptimization("optimizer", 10, [10, 10])
    dftinit = ASGenStage([asgen])
    spillinit = OMCalStage([abacus])
    spillopt = SOptStage([spillage])
    orbgen = TreeProcess([dftinit, spillinit, spillopt])
    orbgen.run()