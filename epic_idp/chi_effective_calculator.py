import numpy as np
from .interaction_scale_functions import interaction_functions
import sys
import matplotlib.pyplot as plt

class chi_effective_calculator:

    def __init__(self, rho0=1. , lB=9437.9 , kappa=11.861 , a=0.24651, interaction_matrix='KH-D', electrostatics_MFT=True, hydrophobicity_MFT=False, electrostatics_RPA=True):
        '''
        Options for the interaction_matrix are:
            'KH-D'
            'Mpipi'
            'Mpipi_RNA' (24-by-24 matrix, includes RNA bases)
            'CALVADOS1'
            'CALVADOS2'
            'HPS'
            'URRY'
            'FB'
        '''

        # All inputs are given in units of the bond length b

        self.rho0  = rho0    # Total concentration. Only provides an over-all scaling of effective chi-parameters.
        self.lB    = lB      # Bjerrum length
        self.kappa = kappa   # Inverse Debye screening length. kappa = sqrt( 8*pi*lB*[NaCl] )
        self.a     = a       # Smearing length for RPA integral. 

        # Interactions to use
        self.electrostatics_MFT = electrostatics_MFT
        self.hydrophobicity_MFT = hydrophobicity_MFT
        self.electrostatics_RPA = electrostatics_RPA

        self.isf = interaction_functions(interaction_matrix)

        self.species = {}  # Dictionary of all molecular species to consider for chi_eff calculations

        # RPA related. y = b*k
        #self.nk = int(4e3)
        self.nk = int(3e2)
        if self.kappa == 0:
            eps = 1e-15
        else:
            eps = 0

        kmax = 2./self.a
        self.k = np.linspace(eps,kmax,self.nk)
        self.k2 = self.k**2
        self.dk = self.k[1] - self.k[0]

        self.Gamma4 = np.exp(-2.*self.a**2 * self.k2)

    # Set params
    def set_rho0(self, rho0):
        self.rho0 = rho0
    
    def set_lB(self, lB):
        self.lB = lB
    
    def set_kappa(self, kappa):
        self.kappa = kappa
    
    def set_a(self, a):
        self.a = a
        self.Gamma4 = np.exp(-2.*a**2 * self.k2)

    # Add a sequence
    def add_IDP(self, name, sequence, g0 = None, verbose=True):
        if verbose:
            print("Adding",name,"... (N="+str(len(sequence))+")")
        self.species[name] = _GaussianChain(self, name, sequence, g0=g0)
        if verbose:
            print("  Done!")

    def remove_IDP(self, name):
        self._check_if_molecules_added(name)
        del self.species[name]

    def _check_if_molecules_added(self, *names):
        for name in names:
            if name not in self.species:
                raise ValueError("[ERROR] Molecule "+name+" has not been added. Added molecules are: "+', '.join( list( self.species.keys()) ) )

    # MFT electrostatics
    def calc_chi_e_MFT(self, name1, name2):
        self._check_if_molecules_added(name1,name2)
        mol1 = self.species[name1]
        mol2 = self.species[name2]

        chi = - 2. * np.pi * self.lB / self.kappa**2 * mol1.sig_tot * mol2.sig_tot / mol1.N / mol2.N
        return chi * self.rho0

    # RPA electrostatics
    def calc_chi_e_RPA(self, name1, name2, ax=None):
        self._check_if_molecules_added(name1,name2)
        mol1 = self.species[name1]
        mol2 = self.species[name2]

        integrand = self.k2 / (self.k2 + self.kappa**2)**2 * self.Gamma4 * mol1.g0 * mol2.g0
        integrand_approx = self.k2 /self.kappa**4 * self.Gamma4 * mol1.g0 * mol2.g0

        chi = 2. * np.pi * self.lB**2 * np.sum(integrand)*self.dk

        if ax is not None:
            ax.plot(self.k, integrand,'-o')
            ax.plot(self.k, integrand_approx,'-')

        return chi * self.rho0
    
    # MFT hydrophobicity (not yet implemented)
    def calc_chi_h(self, name1, name2):
        self._check_if_molecules_added(name1,name2)

        seq1 = self.species[name1].sequence
        seq2 = self.species[name2].sequence

        chi_h = self.isf.chi_eff_pair(seq1,seq2)

        return chi_h * self.rho0
    
    # Calculate total chi_eff
    def calc_chi_eff(self, name1, name2, ax=None):
        chi_e_MFT = self.calc_chi_e_MFT(name1, name2) if self.electrostatics_MFT else 0
        chi_e_RPA = self.calc_chi_e_RPA(name1, name2, ax) if self.electrostatics_RPA else 0
        chi_h     = self.calc_chi_h(name1, name2)     if self.hydrophobicity_MFT else 0

        #print("Electrostatics: RPA/MFT = ",chi_e_RPA/chi_e_MFT)
        #print("MFT:",chi_e_MFT,"RPA:",chi_e_RPA,"Hydrophobicity:",chi_h)

        chi_eff = chi_e_MFT + chi_e_RPA + chi_h
        return chi_eff

    # Calculates to chi_eff matrix of all added species
    def calc_all_chi_eff(self):
        names = np.array( list(self.species.keys()) )

        chi_eff = np.zeros((len(names),len(names)))
        for i in range(len(names)):
            for j in range(i,len(names)):
                chi_eff[i,j] = self.calc_chi_eff(names[i],names[j])
        
        for i in range(len(names)-1):
            for j in range(i+1,len(names)):
                chi_eff[j,i] = chi_eff[i,j]
        
        return chi_eff


class _GaussianChain:
    def __init__(self, cec, name, sequence, g0 = None):
        self.cec      = cec              # chi_effective_calculator object
        self.name     = name             # Name of the molecule
        self.sequence = sequence         # Fasta sequence 
        self.N        = len(sequence)    # Number of residues

        self.sig     = self.cec.isf.charge_from_fasta(self.sequence)
        self.sig_tot = np.sum(self.sig)

        if g0 is None:
            b2k2_over_6 = self.cec.k2 / 6.
            res_diff = np.array( [[ np.abs(alpha-beta) for alpha in range(self.N) ] for beta in range(self.N) ])
            connection_tensor = np.exp( - np.einsum('ab,i->abi', res_diff, b2k2_over_6) )
            self.g0 = np.einsum( 'a,b,abi->i', self.sig, self.sig, connection_tensor ) / self.N
        else:
            self.g0 = g0

        # plt.plot(self.cec.y, self.g)
        # plt.plot(self.cec.y, self.g0)
        # plt.show()
