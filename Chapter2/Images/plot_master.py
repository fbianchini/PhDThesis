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

import healpy as hp
from master import Master
import utils

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

# Params
lmax = 1000
l    = np.arange(lmax+1)
delta_ell = 20

# Masks 
mask_her = hp.read_map('/Users/federicobianchini/Dropbox/Needlet/mask_hatlas_convergence_planck_2015_512.fits', verbose=False)
# mask     = utils.GetGalMask(512, lat=20.)
# fsky     = np.mean(mask)
fsky_her = np.mean(mask_her)
# wl       = hp.anafast(mask, lmax=lmax)
wl_her   = hp.anafast(mask_her, lmax=lmax)
# print 'Galactic mask fsky is:', fsky
print 'H-ATLAS  mask fsky is:', fsky_her

# Master classes
# Gal = Master(2, lmax, delta_ell, mask, fsky_approx=False)
Her = Master(2, lmax, delta_ell, mask_her, fsky_approx=False)

# fig = plt.figure(figsize=(8,6))
# ax1 = fig.add_subplot(111)

# ax1.plot(l, wl, label='Gal Mask')
# ax1.plot(l, wl_her, label='H-ATLAS Mask')
# ax1.set_xlabel(r'$\ell$')
# ax1.set_ylabel(r'$\mathcal{W}_{\ell}$')
# ax1.set_xlim([2,lmax-200])
# ax1.set_yscale('log')
# ax1.legend(loc='best')
# plt.savefig('wl_masks.pdf', bboxes_inches='tight')

# embed()

sns.set_palette("coolwarm",5)

fig  = plt.figure(figsize=(18,8))
ax1  = fig.add_subplot(121)
mat  = ax1.imshow(np.log10(np.abs(Her.M_ll)), interpolation='nearest', cmap='Greys')
cbar = fig.colorbar(mat)
ax1.set_xlabel(r'$\ell_1$')
ax1.set_ylabel(r'$\ell_2$')
ax1.set_title(r'$\log_{10}(|M_{\ell_1\ell_2}|)$')

# sns.palplot(sns.cubehelix_palette(8, start=.5, rot=-.75))
ax2 = fig.add_subplot(122)
for l_ in [5,50,100,300,500]:
	ax2.plot(Her.M_ll[l_,:], label=r'$\ell_2 = %d$' %l_)
ax2.set_xlim([2,lmax-200])
ax2.set_xlabel(r'$\ell_1$')
ax2.set_ylabel(r'$\log_{10}(|M_{\ell_1\ell_2}|)$')
ax2.set_yscale('log')
ax2.set_ylim([1e-9,1e-3])
ax2.legend(loc='best')
plt.savefig('master.pdf', bboxes_inches='tight')

