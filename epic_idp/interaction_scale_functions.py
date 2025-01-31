import numpy as np
import os

class interaction_functions:
    """
    This class contains the functions to calculate the average short-range non-electrostatic interaction strength between two sequences.

    Parameters
    ----------
    scheme_str : str
        The name of the interaction scheme to use. The available schemes are:
        - 'KH-D' : Taken from S3 Data in "Sequence determinants of protein phase behavior from a coarse-grained model", Dignon et al (2018), PLoS Computational Biology, https://doi.org/10.1371/journal.pcbi.1005941
        - 'Mpipi' : Interaction parameters from Joseph et al. Physics-driven coarse-grained model for biomolecular phase separation with near-quantitative accuracy. Nat Comput Sci 1, 732–743 (2021). https://doi.org/10.1038/s43588-021-00155-3. Use this scheme if you are only working with amino-acid sequences and no RNA bases.
        - 'Mpipi_RNA' : Same as 'Mpipi', but with additional interaction parameters for RNA bases (denoted using lower-case symbols 'a', 'c', 'g' and 'u'.).
        - 'CALVADOS1' : CALVADOS1 column of "residues.csv" at https://github.com/KULL-Centre/CALVADOS (27 September, 2022)
        - 'CALVADOS2' : CALVADOS2 column of "residues.csv" at https://github.com/KULL-Centre/CALVADOS (27 September, 2022)
        - 'HPS' : Lambda_i from Table S1 in in "Sequence determinants of protein phase behavior from a coarse-grained model", Dignon et al (2018), PLoS Computational Biology, https://doi.org/10.1371/journal.pcbi.1005941. This corresponds to their “HPS” model.
        - 'URRY' : Urry normalised hydropathy scale given in Table S2 of Regy, Thompson, Kim and Mittal, Protein Science, 2021 (Improved coarse-grained model for studying sequence dependent phase separation of disordered proteins).
        - 'HPS' : Lambda_i from Table S1 in in "Sequence determinants of protein phase behavior from a coarse-grained model", Dignon et al (2018), PLoS Computational Biology, https://doi.org/10.1371/journal.pcbi.1005941. This corresponds to their “HPS” model.

    Attributes
    ----------
    scheme_str : str
        The name of the interaction scheme to use.
    scheme_int : int
        The integer index that identifies the interaction scheme.
    all_res : list
        A list of all residues in the interaction scheme.
    Nres : int
        The number of residues in the interaction scheme.
    res_to_ind : dict
        A dictionary that maps residue names to their integer index.
    residue_charges : np.array
        The electric charges of the residues in the interaction scheme.
    eps_matr : np.array
        The short-range interaction matrix for the interaction scheme.

    Methods
    -------
    charge_from_fasta(seq)
        Calculate the electric charge sequence given a one-letter amino-acid sequence.
    calculate_average_interaction_energy(seq1, seq2)
        Calculate the average short-range non-electrostatic interaction strength between two sequences. Returns the value of $\sum_{i=1}^{N_1} \sum_{j=1}^{N_2} \epsilon_{ij} / (N_1 N_2) $, where $N_1$ and $N_2$ are the lengths of the two sequences, and $\epsilon_{ij}$ is the interaction strength between residues $i$ and $j$.
    """

    def __init__( self, scheme_str ):
        """
        Initialize the interaction_functions object. 

        Parameters
        ----------
        scheme_str : str
            This has to be one of the following strings: 'KH-D', 'Mpipi', 'Mpipi_RNA', 'CALVADOS1', 'CALVADOS2', 'HPS', 'URRY', 'FB'.

        Returns
        -------
        None
        """
        
        all_schemes = {'KH-D'      : 0 ,\
                       'Mpipi'     : 1 ,\
                       'Mpipi_RNA' : 2 ,\
                       'CALVADOS1' : 3 ,\
                       'CALVADOS2' : 4 ,\
                       'HPS'       : 5 ,\
                       'URRY'      : 6 ,\
                       'FB'        : 7
                       }

        base_directory = os.path.dirname(__file__) + '/interaction_matrices/'
        files = {'KH-D'      : 'eps_KH-D.txt'        ,\
                 'Mpipi'     : 'eps_Mpipi.txt'       ,\
                 'Mpipi_RNA' : 'eps_Mpipi_RNA.txt'   ,\
                 'CALVADOS1' : 'lambda_CALVADOS.txt' ,\
                 'CALVADOS2' : 'lambda_CALVADOS.txt' ,\
                 'HPS'       : 'lambda_HPS.txt'      ,\
                 'URRY'      : 'lambda_Urry-HPS.txt' ,\
                 'FB'        : 'lambda_FB.txt'
                 }

        self.scheme_str = scheme_str
        self.scheme_int = int(all_schemes[scheme_str])

        if self.scheme_int == int(all_schemes['Mpipi_RNA']):
            self.all_res = ['R','H','K','D','E','S','T','N','Q','C','G','P','A','V','I','L','M','F','Y','W','a','c','g','u']
        else:
            self.all_res = ['R','H','K','D','E','S','T','N','Q','C','G','P','A','V','I','L','M','F','Y','W']

        res3_to_res1 = { 'ALA' : 'A' ,\
                         'ARG' : 'R' ,\
                         'ASN' : 'N' ,\
                         'ASP' : 'D' ,\
                         'CYS' : 'C' ,\
                         'GLN' : 'Q' ,\
                         'GLU' : 'E' ,\
                         'GLY' : 'G' ,\
                         'HIS' : 'H' ,\
                         'ILE' : 'I' ,\
                         'LEU' : 'L' ,\
                         'LYS' : 'K' ,\
                         'MET' : 'M' ,\
                         'PHE' : 'F' ,\
                         'PRO' : 'P' ,\
                         'SER' : 'S' ,\
                         'THR' : 'T' ,\
                         'TRP' : 'W' ,\
                         'TYR' : 'Y' ,\
                         'VAL' : 'V' }

        self.Nres = len(self.all_res)
        self.res_to_ind = {self.all_res[i] : int(i) for i in range(self.Nres)}

        # Electric charges:
        self.residue_charges = np.zeros(self.Nres)
        if self.scheme_int == int(all_schemes['Mpipi']):
            self.residue_charges[ self.res_to_ind['R'] ] = +0.75
            self.residue_charges[ self.res_to_ind['K'] ] = +0.75
            self.residue_charges[ self.res_to_ind['D'] ] = -0.75
            self.residue_charges[ self.res_to_ind['E'] ] = -0.75
            self.residue_charges[ self.res_to_ind['H'] ] = +0.375
        elif self.scheme_int == int(all_schemes['Mpipi_RNA']):
            self.residue_charges[ self.res_to_ind['R'] ] = +0.75
            self.residue_charges[ self.res_to_ind['K'] ] = +0.75
            self.residue_charges[ self.res_to_ind['D'] ] = -0.75
            self.residue_charges[ self.res_to_ind['E'] ] = -0.75
            self.residue_charges[ self.res_to_ind['H'] ] = +0.375
            self.residue_charges[ self.res_to_ind['a'] ] = -0.75
            self.residue_charges[ self.res_to_ind['c'] ] = -0.75
            self.residue_charges[ self.res_to_ind['g'] ] = -0.75
            self.residue_charges[ self.res_to_ind['a'] ] = -0.75
        else:
            self.residue_charges[ self.res_to_ind['R'] ] = +1.
            self.residue_charges[ self.res_to_ind['K'] ] = +1.
            self.residue_charges[ self.res_to_ind['D'] ] = -1.
            self.residue_charges[ self.res_to_ind['E'] ] = -1
            self.residue_charges[ self.res_to_ind['H'] ] = 0.

        # Short-range interactions
        self.eps_matr = np.zeros((self.Nres,self.Nres))
        file = base_directory + files[self.scheme_str]
        
        if self.scheme_int == int(all_schemes['KH-D']):
            q1, q2, val = np.loadtxt(file,dtype=str).T
            for i in range(len(q1)):
                r1 = res3_to_res1[ q1[i] ]
                r2 = res3_to_res1[ q2[i] ]
                self.eps_matr[ self.res_to_ind[ r1 ] , self.res_to_ind[ r2 ] ] = float(val[i])
                self.eps_matr[ self.res_to_ind[ r2 ] , self.res_to_ind[ r1 ] ] = float(val[i])
        elif self.scheme_int == int(all_schemes['Mpipi']) or self.scheme_int == int(all_schemes['Mpipi_RNA']):
            q1, q2, val, _ = np.loadtxt(file,dtype=str).T
            #matr = np.zeros((self.Nres,self.Nres))
            for i in range(len(q1)):
                self.eps_matr[ self.res_to_ind[ q1[i] ] , self.res_to_ind[ q2[i] ] ] = -float(val[i])
                self.eps_matr[ self.res_to_ind[ q2[i] ] , self.res_to_ind[ q1[i] ] ] = -float(val[i])
        elif self.scheme_int == int(all_schemes['CALVADOS1']) or self.scheme_int == int(all_schemes['CALVADOS2']):
            q, _, lambd1, lambd2 = np.loadtxt(file,dtype=str).T
            if self.scheme_int == int(all_schemes['CALVADOS1']): 
                lambd = lambd1
            else:
                lambd = lambd2
            mu, Delta = 1,0

            val = np.array( [float(l) for l in lambd] )
            for i in range(self.Nres):
                for j in range(self.Nres):
                    r1 = self.res_to_ind[ q[i] ]
                    r2 = self.res_to_ind[ q[j] ]
                    self.eps_matr[r1,r2] = - ( 0.5 * (val[i] + val[j]) * mu - Delta)
        elif self.scheme_int == int(all_schemes['HPS']):
            q, lambd = np.loadtxt(file,dtype=str).T
            mu, Delta = 1,0
            val = np.array( [float(l) for l in lambd] )
            for i in range(self.Nres):
                for j in range(self.Nres):
                    r1 = self.res_to_ind[ res3_to_res1[ q[i] ] ]
                    r2 = self.res_to_ind[ res3_to_res1[ q[j] ] ]
                    self.eps_matr[ r1 , r2 ] = - ( 0.5 * (val[i] + val[j]) * mu - Delta)
        elif self.scheme_int == int(all_schemes['URRY']):
            q, lambd = np.loadtxt(file,dtype=str).T
            mu, Delta = 1, 0.08
            val = np.array( [float(l) for l in lambd] )
            for i in range(self.Nres):
                for j in range(self.Nres):
                    r1 = self.res_to_ind[ res3_to_res1[ q[i] ] ]
                    r2 = self.res_to_ind[ res3_to_res1[ q[j] ] ]
                    self.eps_matr[ r1 , r2 ] = - ( 0.5 * (val[i] + val[j]) * mu - Delta)
        elif self.scheme_int == int(all_schemes['FB']):
            q, lambd = np.loadtxt(file,dtype=str).T
            mu, Delta = 1, 0
            val = np.array( [float(l) for l in lambd] )
            for i in range(self.Nres):
                for j in range(self.Nres):
                    r1 = self.res_to_ind[ res3_to_res1[ q[i] ] ]
                    r2 = self.res_to_ind[ res3_to_res1[ q[j] ] ]
                    self.eps_matr[ r1 , r2 ] = - ( 0.5 * (val[i] + val[j]) * mu - Delta)

    def charge_from_fasta(self,seq):
        """
        Calculate the electric charge sequence given a one-letter amino-acid sequence.
        """
        sig = np.array([ self.residue_charges[self.res_to_ind[r]] for r in seq])
        return sig

    def calculate_average_interaction_energy(self, seq1, seq2):
        """
        Calculate the average short-range non-electrostatic interaction strength between two sequences. Returns the value of $\sum_{i=1}^{N_1} \sum_{j=1}^{N_2} \epsilon_{ij} / (N_1 N_2) $, where $N_1$ and $N_2$ are the lengths of the two sequences, and $\epsilon_{ij}$ is the interaction strength between residues $i$ and $j$.
        """
        chi = np.mean( [[ self.eps_matr[ self.res_to_ind[r1], self.res_to_ind[r2] ] for r1 in seq1] for r2 in seq2] ) 
        return chi









