import bonds, angles, impropers, omega_dihedrals, phi_psi
import chi1, chi2, chi3, chi4

def make_restraints(atmsel, restraints, num_selected):
    for a in (bonds, angles, impropers, omega_dihedrals, phi_psi,
              chi1, chi2, chi3, chi4):
        a.make_restraints(atmsel, restraints, num_selected)
