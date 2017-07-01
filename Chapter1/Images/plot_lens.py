import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import rc, rcParams
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from matplotlib.ticker import NullFormatter
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FormatStrFormatter

from IPython import embed

from classy import Class
import camb
from camb import model, initialpower
from astropy.cosmology import FlatLambdaCDM, LambdaCDM 

sns.set(rc={"figure.figsize": (8, 6)})

def SetPlotStyle():
	rc('text',usetex=True)
	rc('font',**{'family':'serif','serif':['Computer Modern']})
	plt.rcParams['axes.linewidth']  = 3.
	plt.rcParams['axes.labelsize']  = 28
	plt.rcParams['axes.titlesize']  = 22
	plt.rcParams['xtick.labelsize'] = 20
	plt.rcParams['ytick.labelsize'] = 18
	plt.rcParams['xtick.major.size'] = 7
	plt.rcParams['ytick.major.size'] = 7
	plt.rcParams['xtick.minor.size'] = 3
	plt.rcParams['ytick.minor.size'] = 3
	plt.rcParams['legend.fontsize']  = 22
	plt.rcParams['legend.frameon']  = False

	plt.rcParams['xtick.major.width'] = 1
	plt.rcParams['ytick.major.width'] = 1
	plt.rcParams['xtick.minor.width'] = 1
	plt.rcParams['ytick.minor.width'] = 1
	# plt.clf()
	sns.set(rc('font',**{'family':'serif','serif':['Computer Modern']}))
	sns.set_style("ticks", {'figure.facecolor': 'grey'})

SetPlotStyle()

lmax = 3000
l    = np.arange(lmax+1)

# Evaluating CMB T&P spectra
params_lin = {
    'output': 'tCl pCl lCl',
    'l_max_scalars': lmax+500,
    'lensing': 'yes',
    'A_s': 2.3e-9,
    'n_s': 0.9624, 
    'h': 0.6711,
    'omega_b': 0.022068,
    'omega_cdm': 0.12029,
    'Omega_k':0.,
    'z_reio':0}

params_nl = {
    'output': 'tCl pCl lCl',
    'l_max_scalars': lmax+500,
    'lensing': 'yes',
    'A_s': 2.3e-9,
    'n_s': 0.9624, 
    'h': 0.6711,
    'omega_b': 0.022068,
    'omega_cdm': 0.12029,
    'Omega_k':0.,
    'non linear':'halofit'}

cosmo = Class()

# Linear 
cosmo.set(params_lin)
cosmo.compute()
cls_lin   = cosmo.raw_cl(lmax)
lcls_lin  = cosmo.lensed_cl(lmax)

cosmo.struct_cleanup()
cosmo.empty()

# NonLinear Corrections
cosmo.set(params_nl)
cosmo.compute()
cosmo.compute()
cls_nl   = cosmo.raw_cl(lmax)
lcls_nl  = cosmo.lensed_cl(lmax)


fig = plt.figure(figsize=(18,6))
ax1 = fig.add_subplot(121)

ax1.plot(l, (l*(l+1))**2/4.*lcls_lin['pp'], 'grey', ls='--', label='Linear')
ax1.plot(l, (l*(l+1))**2/4.*lcls_nl['pp'], 'black', ls='-', label='NonLinear')

ax1.legend(loc='best')
ax1.set_xlabel(r'$\ell$')
ax1.set_ylabel(r'$C_{\ell}^{\kappa\kappa}$')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim([2,lmax])

# plt.savefig('cmblens.pdf', bbox_inches='tight')

# fig = plt.figure(figsize=(8,6))
ax2 = fig.add_subplot(122)

ax2.plot(lcls_lin['tt']/cls_lin['tt']-1, label='TT')
ax2.plot(lcls_lin['ee']/cls_lin['ee']-1, label='EE')
ax2.axhline(ls='--', color='grey')
ax2.set_xlabel(r'$\ell$')
ax2.set_ylabel(r'$\Delta C_{\ell}/C_{\ell}$')
ax2.legend(loc='best')

plt.savefig('cmblens2.pdf', bbox_inches='tight')







