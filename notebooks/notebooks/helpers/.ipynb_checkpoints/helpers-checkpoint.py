from pyscf import gto, scf, lo, tools, fci


def xyz_to_array(xyz):
    
    split = xyz.split()
    n = int(split[0])
    
    arr = []
    
    for i in range(n):
        
        
        atom = split[2+i*4]
        geom = [float(j) for j in split[3+i*4:3+i*4+3]]
        
        arr.append([atom, geom])
        
        
    return arr

def get_h2_geom(r):
    
    h2_xyz = '2\n*\nH 0.0 0.0 0.0\nH 0.0 0.0 ' + str(r) + '\n'
    
    return h2_xyz

def run_pyscf_hartfree_fock_calculation(xyz, basis="sto-3g", do_fci=False, do_uhf=False, do_break=False):
    """Calculate the energy (+ additional things like MO coefficients) with pyscf."""
    mol = gto.M(
        atom=xyz,
        basis=basis,
        unit="ANG",
        symmetry=True,
    )
    mol.build()
    if do_uhf:
        if do_break:
            mf = scf.uhf.UHF(mol)
            dm_alpha, dm_beta = mf.get_init_guess()
            dm_beta[:2,:2] = 0
            dm = (dm_alpha,dm_beta)
            mf.run(dm, verbose=0)
        else:
            mf = scf.uhf.UHF(mol).run(verbose=0)
    else:
        mf = scf.RHF(mol).run(verbose=0)
    
    if do_fci:
        cisolver = fci.FCI(mf)
        
        fci_energy, ci_coeffs = cisolver.kernel()

        return mf, mol, fci_energy, ci_coeffs

    else:
        return mf, mol


