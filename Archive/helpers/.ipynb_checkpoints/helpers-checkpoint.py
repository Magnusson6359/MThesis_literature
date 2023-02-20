from pyscf import gto, scf, lo, tools


def xyz_to_array(xyz):
    
    split = xyz.split()
    n = int(split[0])
    
    arr = []
    
    for i in range(n):
        
        
        atom = split[2+i*4]
        geom = [float(j) for j in split[3+i*4:3+i*4+3]]
        
        arr.append([atom, geom])
        
        
    return arr


def run_pyscf_hartfree_fock_calculation(xyz, basis="sto-3g"):
    """Calculate the energy (+ additional things like MO coefficients) with pyscf."""
    mol = gto.M(
        atom=xyz,
        basis=basis,
        unit="ANG",
        symmetry=True,
    )
    mol.build()
    mf = scf.RHF(mol).run()
    return mf, mol


