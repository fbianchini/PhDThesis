import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
import seaborn as sns
from matplotlib import rc, rcParams
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from matplotlib.ticker import NullFormatter
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FormatStrFormatter
from IPython import embed

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

kappa_file = '/Volumes/SAMSUNG/Work/SimMap_Lensing/maps_planck/convergence_planck_2015_512.fits'
mask_file  = '/Volumes/SAMSUNG/Work/SimMap_Lensing/masks/mask_convergence_planck_2015_512.fits' 

kappa = hp.read_map(kappa_file)
mask  = hp.read_map(mask_file)

ell   = np.loadtxt('nlkk.dat', usecols=[0], unpack=True)
cl_kk = np.loadtxt('nlkk.dat', usecols=[2], unpack=True)
nl_kk = np.loadtxt('nlkk.dat', usecols=[1], unpack=True)
cl_kk -= nl_kk

kappa_lm      = hp.map2alm(kappa)
phi_lm        = hp.almxfl(kappa_lm, 2./(ell*(ell+1)))
kappa_lm_WF   = hp.almxfl(phi_lm, cl_kk/(cl_kk+nl_kk))
kappa_WF      = hp.ma(hp.alm2map(np.nan_to_num(kappa_lm_WF), 512))
kappa_WF.mask = np.logical_not(mask)

hp.mollview(kappa_WF, title=r'$2015 \quad \rm{Planck}\quad \hat{\phi}^{\rm{WF}}$', min=-7e-5, max=6e-5)
hp.graticule()
plt.savefig('Planck_kappa_CMB.pdf', bbox_inches='tight')

