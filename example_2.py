from epic_idp import chi_effective_calculator
import numpy as np

# Simulation data taken from Joseph et al., Nat Comput Sci 1, 732–743 (2021). https://doi.org/10.1038/s43588-021-00155-3
# Critical temperatures (Mpipi)
Tc_sim = {  '+7R+12D' : 379.1862321467743  ,\
            'WT'      : 340.1549941494668  ,\
            '-12F+12Y': 339.7728588121274  ,\
            '+7K+12D' : 335.3883654686791  ,\
            '+7F-7Y'  : 322.4232219343142  ,\
            '-3R+3K'  : 311.708309889124   ,\
            '-4F-2Y'  : 290.01317430057094 ,\
            '-6R+6K'  : 283.9045195238529  ,\
            '-9F+3Y'  : 295.4001291149048
            }
# Standard deviations of the critical temperatures (Mpipi)
Tc_sim_std = {  '+7R+12D' : 1.09 ,\
                'WT'      : 1.19 ,\
                '-12F+12Y': 2.30  ,\
                '+7K+12D' : 0.86 ,\
                '+7F-7Y'  : 1.47 ,\
                '-3R+3K'  : 2.58 ,\
                '-4F-2Y'  : 0.74 ,\
                '-6R+6K'  : 5.20 ,\
                '-9F+3Y'  : 0.36
            }

seqs = {}
seqs['WT']      =   'MASASSSQRG'+'RSGSGNFGGG'+'RGGGFGGNDN'+'FGRGGNFSGR'+'GGFGGSRGGG'+\
                    'GYGGSGDGYN'+'GFGNDGSNFG'+'GGGSYNDFGN'+'YNNQSSNFGP'+'MKGGNFGGRS'+\
                    'SGPYGGGGQY'+'FAKPRNQGGY'+'GGSSSSSSYG'+'SGRRF'
seqs['-3R+3K']  =   'MASASSSQRG'+'KSGSGNFGGG'+'RGGGFGGNDN'+'FGRGGNFSGR'+'GGFGGSKGGG'+\
                    'GYGGSGDGYN'+'GFGNDGSNFG'+'GGGSYNDFGN'+'YNNQSSNFGP'+'MKGGNFGGRS'+\
                    'SGGSGGGGQY'+'FAKPRNQGGY'+'GGSSSSSSYG'+'SGRKF'
seqs['-4F-2Y']  =   'MASASSSQRG'+'RSGSGNSGGG'+'RGGGFGGNDN'+'FGRGGNSSGR'+'GGFGGSRGGG'+\
                    'GYGGSGDGYN'+'GFGNDGSNSG'+'GGGSSNDFGN'+'YNNQSSNFGP'+'MKGGNFGGRS'+\
                    'SGGSGGGGQY'+'SAKPRNQGGY'+'GGSSSSSSSG'+'SGRRF'
seqs['-6R+6K']  =   'MASASSSQKG'+'KSGSGNFGGG'+'RGGGFGGNDN'+'FGKGGNFSGR'+'GGFGGSKGGG'+\
                    'GYGGSGDGYN'+'GFGNDGSNFG'+'GGGSYNDFGN'+'YNNQSSNFGP'+'MKGGNFGGKS'+\
                    'SGGSGGGGQY'+'FAKPRNQGGY'+'GGSSSSSSYG'+'SGRKF'
seqs['+7F-7Y']  =   'MASASSSQRG'+'RSGSGNFGGG'+'RGGGFGGNDN'+'FGRGGNFSGR'+'GGFGGSRGGG'+\
                    'GFGGSGDGFN'+'GFGNDGSNFG'+'GGGSFNDFGN'+'FNNQSSNFGP'+'MKGGNFGGRS'+\
                    'SGGSGGGGQF'+'FAKPRNQGGF'+'GGSSSSSSFG'+'SGRRF'
seqs['+7K+12D'] =   'MASADSSQRD'+'RDDKGNFGDG'+'RGGGFGGNDN'+'FGRGGNFSDR'+'GGFGGSRGDG'+\
                    'KYGGDGDKYN'+'GFGNDGKNFG'+'GGGSYNDFGN'+'YNNQSSNFDP'+'MKGGNFKDRS'+\
                    'SGPYDKGGQY'+'FAKPRNQGGY'+'GGSSSSKSYG'+'SDRRF'
seqs['+7R+12D'] =   'MASADSSQRD'+'RDDRGNFGDG'+'RGGGFGGNDN'+'FGRGGNFSDR'+'GGFGGSRGDG'+\
                    'RYGGDGDRYN'+'GFGNDGRNFG'+'GGGSYNDFGN'+'YNNQSSNFDP'+'MKGGNFRDRS'+\
                    'SGPYDRGGQY'+'FAKPRNQGGY'+'GGSSSSRSYG'+'SDRRF'
seqs['-9F+3Y']  =   'MASASSSQRG'+'RSGSGNFGGG'+'RGGGYGGNDN'+'GGRGGNYSGR'+'GGFGGSRGGG'+\
                    'GYGGSGDGYN'+'GGGNDGSNYG'+'GGGSYNDSGN'+'GNNQSSNFGP'+'MKGGNYGGRS'+\
                    'SGGSGGGGQY'+'GAKPRNQGGY'+'GGSSSSSSYG'+'SGRRS'
seqs['-12F+12Y']=   'MASASSSQRG'+'RSGSGNYGGG'+'RGGGYGGNDN'+'YGRGGNYSGR'+'GGYGGSRGGG'+\
                    'GYGGSGDGYN'+'GYGNDGSNYG'+'GGGSYNDYGN'+'YNNQSSNYGP'+'MKGGNYGGRS'+\
                    'SGGSGGGGQY'+'YAKPRNQGGY'+'GGSSSSSSYG'+'SGRRY'
#seq_names = list(seqs.keys()) # All sequence names

# Order seq_names according to T_c
idx = np.argsort( np.array(list(Tc_sim.values())) )[::-1]
seq_names = np.array(list(Tc_sim.keys()))[idx]

######################################################################
######### Use EPIC-IDP to calculate effective chi parameters #########
######################################################################

# Create the chi_effective_calculator object
cec = chi_effective_calculator(rho0=1., lB=1.7, kappa=0.75, a=0.1, Vh0=3.0, interaction_matrix='Mpipi')

# Add the A1 sequences
for seq_name in seq_names:
    cec.add_IDP(seq_name,seqs[seq_name])

chi_eff_matrix = cec.calc_all_chi_eff()


######################################################################
###################### Visualize the results #########################
######################################################################

import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['text.usetex'] = True

fig, axes = plt.subplots(figsize=(8,3.),ncols=4, gridspec_kw={'width_ratios':[1,0.27,1,0.06], 'wspace':0.01})

ax = axes[0]
x_TC = np.array([Tc_sim[seq_name] for seq_name in seq_names])
std_Tc = np.array([Tc_sim_std[seq_name] for seq_name in seq_names])
chi_ii = np.diag(chi_eff_matrix)

# Annotated points
for i, sn in enumerate(seq_names):
    ax.errorbar(chi_ii[i], Tc_sim[sn], yerr=Tc_sim_std[sn], fmt='o', color='C'+str(i),mec='k',mew=0.5, capsize=3, capthick=2.)

    # Annotate the points
    X = chi_ii[i]+0.007
    Y = Tc_sim[sn]-4
    if sn == '+7R+12D':
        X -= 0.07
    if sn == 'WT':
        X -= 0.033
        Y += 3.
    if sn == '-12F+12Y':
        Y += 2
    if sn == '-4F-2Y':
        X += 0.01
    if sn == '-9F+3Y':
        Y += 1.
    ax.text(X, Y, sn, fontsize=10, color='C'+str(i), va='bottom')

ax.set_xlabel(r'$\chi_{ii}$')
ax.set_ylabel(r'$T_c / {\rm K}$ (Mpipi)')

axes[1].axis('off')

ax = axes[2]
cax = ax.imshow(chi_eff_matrix,cmap='rainbow',aspect='equal',origin='lower')
cbar = fig.colorbar(cax,cax=axes[3])
cbar.set_label(r'$\chi_{ij}$')

# Ticks and tick labels are the sequence names
ax.set_xticks(np.arange(len(seq_names)))
ax.set_yticks(np.arange(len(seq_names)))
ax.set_xticklabels(seq_names,rotation=90, fontsize=8)
ax.set_yticklabels(seq_names, fontsize=8)

ax.set_xlabel(r'Sequence $i$')
ax.set_ylabel(r'Sequence $j$')

plt.savefig('example_2.png',dpi=300, bbox_inches='tight')
plt.show()
