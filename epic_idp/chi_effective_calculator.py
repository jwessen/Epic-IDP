import numpy as np
from .interaction_scale_functions import interaction_functions

class chi_effective_calculator:

    """
    A class for calculating effective Flory-Huggins interaction parameters. 

    Parameters
    ----------
    rho0 : float, optional
        Reference concentration. Only provides an over-all scaling of effective chi-parameters. The default is 1.
    lB : float, optional
        Bjerrum length in units of bond length. The default is 9437.9.
    kappa : float, optional
        Inverse Debye screening length in units of inverse bond length. The default is 11.861.
    a : float, optional
        Smearing length for RPA integral in units of bond length. The default is 0.24651.
    Vh0 : float, optional
        Volume integral of short-range interaction potential. Setting Vh0=0 means that non-electrostatic interactions are neglected. The default is 0.
    interaction_matrix : str, optional
        Options for the interaction_matrix are:
            'KH-D'
            'Mpipi'
            'Mpipi_RNA' (24-by-24 matrix, includes RNA bases)
            'CALVADOS1'
            'CALVADOS2'
            'HPS'
            'URRY'
            'FB'
        The default is 'KH-D'.
    electrostatics_RPA : bool, optional
        Whether to include RPA electrostatics. The default is True.

    Methods
    ----------
    add_IDP(name, sequence, g0=None, verbose=True)
        Add a sequence to the chi_effective_calculator object.
    remove_IDP(name)
        Remove a sequence from the chi_effective_calculator object.
    calc_chi_eff(name1, name2)
        Calculate the effective Flory-Huggins interaction parameter between two sequences. The sequences must be added to the chi_effective_calculator object.
    calc_all_chi_eff()
        Calculate the matrix of effective Flory-Huggins interaction parameters between all added sequences.
    set_rho0(rho0)
        Set the reference concentration.
    set_lB(lB)
        Set the Bjerrum length.
    set_kappa(kappa)
        Set the inverse Debye screening length.
    set_a(a)
        Set the smearing length for RPA integral.
    set_Vh0(Vh0)
        Set the volume integral of short-range interaction potential.
    """

    def __init__(self, rho0=1. , lB=9437.9 , kappa=11.861 , a=0.24651, Vh0 = 0, interaction_matrix='KH-D', electrostatics_RPA=True):
        """
        Initialize the chi_effective_calculator object. 

        Parameters
        ----------
        rho0 : float, optional
            Reference concentration. Only provides an over-all scaling of effective chi-parameters. The default is 1.
        lB : float, optional
            Bjerrum length in units of bond length. The default is 9437.9.
        kappa : float, optional
            Inverse Debye screening length in units of inverse bond length. The default is 11.861.
        a : float, optional
            Smearing length for RPA integral in units of bond length. The default is 0.24651.
        Vh0 : float, optional
            Volume integral of short-range interaction potential. Setting Vh0=0 means that non-electrostatic interactions are neglected. The default is 0.
        interaction_matrix : str, optional
            Only used when Vh0 is not zero. Options for the interaction_matrix are:
                'KH-D'
                'Mpipi'
                'Mpipi_RNA' (24-by-24 matrix, includes RNA bases)
                'CALVADOS1'
                'CALVADOS2'
                'HPS'
                'URRY'
                'FB'
            The default is 'KH-D'.
        electrostatics_RPA : bool, optional
            Whether to include the RPA contribution electrostatics. The default is True.
        """

        # All inputs are given in units of the bond length b

        self.rho0  = rho0    # Total concentration. Only provides an over-all scaling of effective chi-parameters.
        self.lB    = lB      # Bjerrum length
        self.kappa = kappa   # Inverse Debye screening length. kappa = sqrt( 8*pi*lB*[NaCl] )
        self.a     = a       # Smearing length for RPA integral.
        self.Vh0   = Vh0     # Volume integral of short-range interaction potential (Setting Vh0=0 means that non-electrostatic interactions are neglected)

        # Interactions to use
        if Vh0 == 0:
            self.hydrophobicity_MFT = False
        else:
            self.hydrophobicity_MFT = True
        self.electrostatics_MFT = True
        self.electrostatics_RPA = electrostatics_RPA

        self.isf = interaction_functions(interaction_matrix)

        self.species = {}  # Dictionary of all molecular species to consider for chi_eff calculations

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
        """
        Set the reference concentration rho0.
        """
        self.rho0 = rho0
    
    def set_lB(self, lB):
        """
        Set the Bjerrum length lB.
        """
        self.lB = lB
    
    def set_kappa(self, kappa):
        """
        Set the inverse Debye screening length kappa.
        """
        self.kappa = kappa
    
    def set_a(self, a):
        """
        Set the smearing length for RPA integral a.
        """
        self.a = a
        self.Gamma4 = np.exp(-2.*a**2 * self.k2)
    
    def set_Vh0(self, Vh0):
        """
        Set the volume integral of short-range interaction potential Vh0.
        """
        self.Vh0 = Vh0
        if Vh0 == 0:
            self.hydrophobicity_MFT = False
        else:
            self.hydrophobicity_MFT = True

    # Add a sequence
    def add_IDP(self, name, sequence, g0 = None, verbose=True):
        """
        Add a sequence to the chi_effective_calculator object.

        Parameters
        ----------
        name : str
            Name of the sequence.
        sequence : str
            Fasta sequence (one-letter code).
        g0 : numpy.ndarray, optional
            Electric charge form factor $g(k) = N^{-1} \sum_{i,j=1}^N \sigma_i \sigma_j \exp(-|i-j|b^2 k^2/ 6. )$. Computes it if not provided.
        verbose : bool, optional
            Print progress. The default is True.

        Returns
        -------
        None.

        """
        if verbose:
            print("Adding",name,"... (N="+str(len(sequence))+")")
        self.species[name] = _GaussianChain(self, name, sequence, g0=g0)
        if verbose:
            print("  Done!")

    def remove_IDP(self, name):
        """
        Remove a sequence from the chi_effective_calculator object.

        Parameters
        ----------
        name : str
            Name of the sequence to remove.

        Returns
        -------
        None.

        """
        self._check_if_molecules_added(name)
        del self.species[name]

    def _check_if_molecules_added(self, *names):
        for name in names:
            if name not in self.species:
                raise ValueError("[ERROR] Molecule "+name+" has not been added. Added molecules are: "+', '.join( list( self.species.keys()) ) )

    # MFT electrostatics
    def _calc_chi_e_MFT(self, name1, name2):
        self._check_if_molecules_added(name1,name2)
        mol1 = self.species[name1]
        mol2 = self.species[name2]

        if mol1.sig_tot * mol2.sig_tot != 0:
            chi = - 2. * np.pi * self.lB / self.kappa**2 * mol1.sig_tot * mol2.sig_tot / mol1.N / mol2.N
        else:
            chi = 0
        return chi * self.rho0

    # RPA electrostatics
    def _calc_chi_e_RPA(self, name1, name2):
        self._check_if_molecules_added(name1,name2)
        mol1 = self.species[name1]
        mol2 = self.species[name2]

        integrand = self.k2 / (self.k2 + self.kappa**2)**2 * self.Gamma4 * mol1.g0 * mol2.g0
        chi = 2. * np.pi * self.lB**2 * np.sum(integrand)*self.dk

        return chi * self.rho0
    
    # MFT hydrophobicity (or other short-range non-electrostatic interactions)
    def _calc_chi_h(self, name1, name2):
        self._check_if_molecules_added(name1,name2)

        seq1 = self.species[name1].sequence
        seq2 = self.species[name2].sequence

        chi_h = -0.5*self.Vh0 * self.isf.calculate_average_interaction_energy(seq1,seq2)

        return chi_h * self.rho0
    
    def calc_chi_eff(self, name1, name2):
        """
        Calculate the effective Flory-Huggins interaction parameter between two sequences. 

        Parameters
        ----------
        name1 : str
            Name of the first sequence. (Must be added to the chi_effective_calculator object)
        name2 : str
            Name of the second sequence. (Must be added to the chi_effective_calculator object)

        Returns
        -------
        chi_eff : float
            The effective Flory-Huggins interaction parameter between the two sequences.
        """

        chi_e_MFT = self._calc_chi_e_MFT(name1, name2) if self.electrostatics_MFT else 0
        chi_e_RPA = self._calc_chi_e_RPA(name1, name2) if self.electrostatics_RPA else 0
        chi_h     = self._calc_chi_h(name1, name2)     if self.hydrophobicity_MFT else 0

        chi_eff = chi_e_MFT + chi_e_RPA + chi_h
        return chi_eff

    # Calculates to chi_eff matrix of all added species
    def calc_all_chi_eff(self):
        """
        Calculate the effective Flory-Huggins interaction parameter matrix between all added sequences.

        Returns
        -------
        chi_eff : numpy.ndarray
            The effective Flory-Huggins interaction parameter matrix between all added sequences. 
        """
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
