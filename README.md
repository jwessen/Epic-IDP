## EPIC-IDP: A Python tool for calculating effective interactions between intrinsically disordered proteins

### About

EPIC-IDP (Effective Protein Interaction Calculator for Intrinsically Disordered Proteins) is a Python package for calculating interaction strengths between intrinsically disordered proteins (IDPs), as quantified by a matrix of effective Flory-Huggins $\chi_{ij}$ parameters, using the IDP amino-acid sequences as input. The $\chi_{ij}$ parameters give information about the propensity to form phase-separated biomolecular condensates and about the partitioning of client proteins inside these condensates. The program accounts for short-range interactions using a mean-field treatment and for long-range electrostatic interactions using the random phase approximation (RPA) theory.

### Background

The main function of the program is to compute $\chi_{ij}$ using the following equation:

$$
\chi_{ij} = \left( \chi_{\rm h}^{(0)} \right)_{ij} + \left( \chi_{\rm e}^{(0)} \right)_{ij} + \left( \chi_{\rm e}^{(1)} \right)_{ij} 
$$
where the three terms correspond to:
- $\left( \chi_{\rm h}^{(0)} \right)_{ij}$: Effective $\chi_{ij}$ parameter following from a mean-field treatment of short-range non-electrostatic interactions (e.g., hydrophobic interactions or cation-$\pi$ interactions). This only depends on the amino-acid content (composition), but not the residue order (sequence), of the involved proteins.
- $\left( \chi_{\rm e}^{(0)} \right)_{ij}$: Effective $\chi_{ij}$ parameter following from a mean-field treatment of long-range electrostatic interactions. This only depends on the net charge per chain of the two proteins.
- $\left( \chi_{\rm e}^{(1)} \right)_{ij}$: The first order correction from electrostatic interactions that follows from RPA theory. This term accounts for charge sequence patterns in the amino-acid sequences, and can thus distinguish between IDPs with same composition but different sequences.

### Usage

To use the package, you first define a `chi_effective_calculator` instance which takes interaction parameter values as input, e.g.

```python
from epic_idp import chi_effective_calculator

cec = chi_effective_calculator( rho0  = 1.0,
                                lB    = 2.0,
                                kappa = 0.1,
                                a     = 0.15,
                                Vh0   = 0,
                                interaction_matrix = 'KH-D' )
```

The arguments correspond to:

- `rho0` $\leftrightarrow \rho_0 b^3$
- `lB` $\leftrightarrow l_{\rm B} / b$
- `kappa` $\leftrightarrow \kappa b$
- `a` $\leftrightarrow a$
- `Vh0` $\leftrightarrow \int {\rm d} {\bf r} V_{\rm h}(|{\bf r}|)$

The `interaction_matrix` argument refers to the data-set used for contact energies $\epsilon_{r,r'}$ and has to be one of the following:

- `'KH-D'` (S3 Data in Dignon et al., 2018)
- `'Mpipi'` (using the 20-by-20 matrix for amino-acid pairs, Joseph et al., 2021)
- `'Mpipi_RNA'` (Using the 24-by-24 matrix including RNA bases, Joseph et al., 2021)
- `'CALVADOS1'` (Using the original CALVADOS dataset, Tesei et al., 2021)
- `'CALVADOS2'` (Using the updated CALVADOS dataset, Tesei et al., 2022)
- `'HPS'` (Table S1 in Dignon et al., 2018)
- `'URRY'` (Table S2 in Regy et al., 2021)
- `'FB'` (Table S7 in Dannenhoffer-Lafage et al., 2021)

After creating the `chi_effective_calculator` instance, we define which IDR sequences we are interested in. The IDR sequences are added sequentially with no upper limit on how many sequences can be added. For example:

```python
# Sequences in FASTA format
sequences = [ 'EHHSGSQGPLLTTGDLGKEKTQ', 'RKQDEEERELRAKQEQEKELLRQKLLKEQEEK', 'AGREAKRR' ]

# Names of the IDPs
names = [ 'seq1', 'seq2', 'seq3' ]

for seq, name in zip(sequences, names):
   cec.add_IDP(name, seq)
```

When adding an IDP, the program directly computes the double sum over residues in $g_i(k)$, defined in Eq.~\eqref{eq:g_def} which constitutes the most computationally heavy step of the effective $\chi$ parameter calculations. After all IDPs have been added, the $\chi_{ij}$ elements (with an overall factor $\rho_0$) are computed using the functions of the `chi_effective_calculator` instance, e.g.

```python
# Returns the effective chi-parameter for the seq1-seq2 pair
cec.calc_chi_eff('seq1', 'seq2')

# Returns the full M-by-M matrix of chi parameter (M is the number of added IDPs)
cec.calc_all_chi_eff()
```