import numpy as np
import matplotlib.pyplot as plt
import camb
from camb import model, initialpower
from astropy.cosmology import FlatLambdaCDM, LambdaCDM 
import seaborn as sns
from matplotlib import rc, rcParams
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from matplotlib.ticker import NullFormatter
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FormatStrFormatter
from IPython import embed
import healpy as hp
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

# Cosmo pars
Omega_r = (2.59e-5/.67**2+1.69e-5/.67**2)
Omega_b = 0.05
Omega_c = 0.3
Omega_m = Omega_c + Omega_b
Omega_v = 1. - Omega_r - Omega_m

# Flat LCDM cosmology 
pars = camb.CAMBparams()
pars.set_cosmology(H0=70, ombh2=Omega_b*.7**2, omch2=Omega_c*.7**2, mnu=0., omk=0)
pars.InitPower.set_params(ns=0.965, r=0)

results = camb.get_results(pars)

# embed()

# Distances
z   = np.logspace(-4,6,10000)
a   = 1./(1.+z)
eta = results.conformal_time(z)

E2_z = lambda x: Omega_m*(1+z)**3 + Omega_r*(1+z)**4 + Omega_v

# plt.subplot(121)
plt.loglog(eta, a, lw=2, color='royalblue')

plt.loglog(eta, Omega_r*(1+z)**4/E2_z(z), 'k--')
plt.loglog(eta, Omega_m*(1+z)**3/E2_z(z), 'k--')
plt.loglog(eta, Omega_v/E2_z(z), 'k--')

plt.axvline(1e2, ls=':', color='indianred')
plt.axvline(results.conformal_time(1./.44 - 1), ls=':', color='indianred')

plt.text(20, 1.3e-6, r'Radiation', fontsize=20)
plt.text(150, 1.3e-6, r'Matter', fontsize=20)
plt.text(10500, 1.3e-6, r'$\Lambda$', fontsize=20)

plt.text(10, 1.3, r'$\Omega_{\rm r}$', fontsize=20)
plt.text(10, 0.02, r'$\Omega_{\rm m}$', fontsize=20)
plt.text(3000, 2e-5, r'$\Omega_{\Lambda}$', fontsize=20)

plt.text(20, 3e-4, r'$a(\eta)\propto \eta$', fontsize=20, rotation=25)
plt.text(300, 1.4e-2, r'$a(\eta)\propto \eta^2$', fontsize=20, rotation=35)

# plt.setp(plt.gca(), 'xticklabels',[0.01, 0.1, 1, 10])

plt.xlabel(r'Conformal time $\eta$ [Mpc]')
plt.ylabel(r'Scale factor $a$')
plt.legend(loc='best')
plt.ylim([1e-6,4])
plt.xlim([5,14000])
# plt.show()
plt.savefig('dens.pdf', bbox_inches='tight')


