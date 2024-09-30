'''Read ABACUS-ORBGEN input script
The logic of ABACUS-ORBGEN input script is a little complicated, because 
the denpendecy between orbitals and ABACUS run is not direct enough.

There are mainly two parts needed to be defined in input, the set of
orbitals and set of structures (reference systems).

orbitals
--------
The orbitals can be uniquely defined by `ecut`, `rcut`, `nzeta`. Its value
can be uniquely defined if with a set of reference systems.

structures
----------
The generation of orbitals needs extract information from the calculation
of structures. For pw cases, the maximal angular momentum over all atom
types lmaxmax should be explicitly defined, this introduce the dependency
of calculation on orbitals.

'''

from SIAB.io.read_input import _cal_nzeta, abacus_params, cal_nbands_fill_lmax, nbands_from_str
import unittest
def _set_orbparam(**kwargs):
    '''set the parameters for orbital to generate
    '''
    out = {}

    of_interest = ['shape', 'zeta_notation', 'nbands_ref', 'orb_ref']
    of_interest = {key: kwargs[key] for key in of_interest}
    
    # nzeta: can be a list of integers, [S]s[P]p[D]d[F]f..., or 'auto'
    nz = of_interest['zeta_notation']
    minbas = kwargs.get('minbas')
    if minbas is None: # the case that nzeta is not needed to be calculated
        # really?
        general = isinstance(nz, list) and all([isinstance(i, int) for i in nz])
        infer = nz == 'auto'
        assert general or infer, f'invalid zeta_notation: {nz}'
    else:
        nz = _cal_nzeta(nz, minbas)
    out['nzeta'] = nz

    # nzeta_from: can be an integer or 'none'
    orbr = of_interest['orb_ref']
    assert isinstance(orbr, int) or orbr == "none", f'invalid orb_ref: {orbr}'
    out['nzeta_from'] = orbr

    # nbands_ref: can be an integer, list of integers, or 'occ', 'all', ...
    nbndr = of_interest['nbands_ref']
    assert isinstance(nbndr, int) or\
           (isinstance(nbndr, list) and all([isinstance(i, int) for i in nbndr])) or\
           nbndr in ['occ', 'all'], f'invalid nbands_ref: {nbndr}'
    out['nbands_ref'] = nbndr

    # shape: can be a string, integer or list of either of them
    shape = of_interest['shape']
    assert isinstance(shape, str) or isinstance(shape, int) or\
           (isinstance(shape, list) and all([isinstance(i, (str, int)) for i in shape])),\
           f'invalid shape: {shape}'
    out['shape'] = shape

    return out

def _parse_orb_block(block):
    '''parse the orbitals block in ABACUS-ORBGEN input file
    
    Parameters
    ----------
    block : list
        list of dict, each dict defines an orbital

    Returns
    -------
    list
        list of dict, each dict defines an orbital
    '''
    return [_set_orbparam(**i) for i in block]

def _collect_abacus_keys(data):
    '''collect all keys that are related to ABACUS INPUT
    
    Parameters
    ----------
    data : dict
        input data

    Returns
    -------
    dict
        dict of keys with the corresponding values
    '''
    keys = abacus_params() + ['pseudo_name']
    return {k: v for k, v in data.items() if k in keys}

def _lmaxmax_infer(orbparams, sysparams):
    '''set the lmaxmax based on the nzeta and dependencies the orbital
    on systems
    
    Parameters
    ----------
    orbparams : list
        list of dict, each dict defines an orbital, should be the output of
        _parse_orb_block
    sysparams : list
        list of dict, each dict defines a system.
        
    Returns
    -------
    list
        list of dict, each dict defines a system, with lmaxmax set
    '''
    shapes = [s['shape'] for s in sysparams]
    for orb in orbparams:
        s = orb['shape']
        nz = orb['nzeta']
        lmaxmax = len(nz) - 1 if isinstance(nz, list) else None
        lmaxmax = "__infer__" if nz == "auto" else lmaxmax
        if isinstance(s, str):
            l_ = sysparams[shapes.index(s)].get('lmaxmax', -1)
            sysparams[shapes.index(s)]['lmaxmax'] = max(l_, lmaxmax)
        elif isinstance(s, int):
            l_ = sysparams[s].get('lmaxmax', -1)
            sysparams[s]['lmaxmax'] = max(l_, lmaxmax)
        elif isinstance(s, list) and all([isinstance(i, int) for i in s]):
            for i in s:
                l_ = sysparams[i].get('lmaxmax', -1)
                sysparams[i]['lmaxmax'] = max(l_, lmaxmax)
        else:
            assert False, f'invalid shape in orbparam: {s}'
    for i, s in enumerate(sysparams):
        if s['lmaxmax'] is None:
            assert False, f'lmaxmax is not set for the {i}-th system'
    return sysparams

def _nbands_infer(sysparam, zval = None):
    '''infer the nbands from zval
    
    Parameters
    ----------
    sysparam : dict
        system parameters
    zval : float
        valence charge of the system
    '''
    if zval is None:
        return sysparam
    nbnd = sysparam.get('nbands', "auto")
    sysparam['nbands'] = nbands_from_str(nbnd, sysparam['shape'], zval)
    return sysparam

def _render_abacus_input(general, sysparam):
    '''generate abacus input parameters based on general and specific
    setting on system
    
    Parameters
    ----------
    general : dict
        general settings for abacus
    sysparam : dict
        specific settings for abacus, if a key is not defined in sysparam,
        the value in general will be used, otherwise the value in sysparam
        will be used.
    
    Returns
    -------
    dict
        abacus input parameters
    '''
    out = general.copy()
    overwrite = _collect_abacus_keys(sysparam)
    return out|overwrite

def _collect_env_keys(data):
    '''collect settings for calculation'''
    ext = ['environment', 'mpi_command', 'abacus_command']
    return {k: data.get(k, "") for k in ext}

def _make_monomer(orbparams, zval = None, zcore = None):
    '''make a guess of spill based on the orbitals'''
    out = {
        'shape': 'monomer',
        'nbands': None,
        'nspin': 1,
        'lmaxmax': None
    }

    lmaxmax = None
    for orb in orbparams:
        if orb['nzeta'] == 'auto':
            out['nbands'] == '__infer__'
            out['lmaxmax'] == '__infer__'
            break
        else:
            assert isinstance(orb['nzeta'], list), 'nzeta should be a list'
            lmaxmax = max(lmaxmax, len(orb['nzeta']) - 1)

    assert lmaxmax is not None or out['lmaxmax'] == '__infer__'
    if lmaxmax is not None:
        out['lmaxmax'] = lmaxmax
        assert zval is not None and zcore is not None, 'zval and zcore should be provided'
        out['nbands'] = cal_nbands_fill_lmax(zval, zcore, lmaxmax)
    return out

def parse(data, zval = None, zcore = None):
    '''main interface to parse the JSON input of ABACUS-ORBGEN
    
    Parameters
    ----------
    data : dict
        input data
    '''
    orb = _parse_orb_block(data['orbitals'])
    sys = _lmaxmax_infer(orb, data['reference_systems'])
    sys = [_nbands_infer(i, zval) for i in sys]
    abacus_general = _collect_abacus_keys(data)
    abacus = [_render_abacus_input(abacus_general, i) for i in sys]
    env = _collect_env_keys(data)

class TestUnpack(unittest.TestCase):

    def test_set_orbparam(self):

        orb = {
            'shape': 'dimer',
            'zeta_notation': 'auto',
            'nbands_ref': 'occ',
            'orb_ref': 'none'
        }
        self.assertEqual(_set_orbparam(**orb), {'nzeta': 'auto', 'nzeta_from': 'none', 'nbands_ref': 'occ', 'shape': 'dimer'})

if __name__ == '__main__':
    unittest.main()