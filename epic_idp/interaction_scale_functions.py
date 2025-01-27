import numpy as np
import os

class interaction_functions:
    def __init__( self, scheme_str ):
        
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
                 'Mpipi_RNA' : 'eps_Mpipi.txt'       ,\
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
            #matr = np.zeros((self.Nres,self.Nres))
            for i in range(self.Nres):
                for j in range(self.Nres):
                    r1 = self.res_to_ind[ res3_to_res1[ q[i] ] ]
                    r2 = self.res_to_ind[ res3_to_res1[ q[j] ] ]
                    self.eps_matr[ r1 , r2 ] = - ( 0.5 * (val[i] + val[j]) * mu - Delta)
        elif self.scheme_int == int(all_schemes['URRY']):
            q, lambd = np.loadtxt(file,dtype=str).T
            mu, Delta = 1, 0.08
            val = np.array( [float(l) for l in lambd] )
            #matr = np.zeros((self.Nres,self.Nres))
            for i in range(self.Nres):
                for j in range(self.Nres):
                    r1 = self.res_to_ind[ res3_to_res1[ q[i] ] ]
                    r2 = self.res_to_ind[ res3_to_res1[ q[j] ] ]
                    self.eps_matr[ r1 , r2 ] = - ( 0.5 * (val[i] + val[j]) * mu - Delta)
        elif self.scheme_int == int(all_schemes['FB']):
            q, lambd = np.loadtxt(file,dtype=str).T
            mu, Delta = 1, 0
            val = np.array( [float(l) for l in lambd] )
            #matr = np.zeros((self.Nres,self.Nres))
            for i in range(self.Nres):
                for j in range(self.Nres):
                    r1 = self.res_to_ind[ res3_to_res1[ q[i] ] ]
                    r2 = self.res_to_ind[ res3_to_res1[ q[j] ] ]
                    self.eps_matr[ r1 , r2 ] = - ( 0.5 * (val[i] + val[j]) * mu - Delta)

        self.lambd, self.Q = np.linalg.eig( self.eps_matr ) # lambd - eigenvalues,  Q - eigenvectors

        # sort according to |lambda| (largest first)
        I = np.flip( np.argsort(np.abs(self.lambd)) )
        self.lambd = self.lambd[I]
        self.Q     = self.Q[:,I]


    def charge_from_fasta(self,seq):
        sig = np.array([ self.residue_charges[self.res_to_ind[r]] for r in seq])
        return sig

    def chi_eff_pair(self, seq1, seq2, norm=1.):
        chi = np.sum( [[ self.eps_matr[ self.res_to_ind[r1], self.res_to_ind[r2] ] for r1 in seq1] for r2 in seq2] ) / np.abs(norm)
        return chi

    def chi_eff(self, seq, norm=1.):
        return self.chi_eff_pair(seq,seq,norm)

    def chi_matr(self, seqs, norm=1.):
        chi = np.array([[ self.chi_eff_pair(s1,s2,norm) for s1 in seqs] for s2 in seqs])
        return chi
    
    def eps_pair(self,res1,res2):
        return self.eps_matr[ self.res_to_ind[res1], self.res_to_ind[res2] ]









