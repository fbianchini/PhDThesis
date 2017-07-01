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

# Astropy Flat LCDM cosmology 
cosmoLCDM = FlatLambdaCDM(H0=70, Om0=0.35, Ob0=0.05)
print 'cosmo', cosmoLCDM, '\n'

# Astropy Flat SCDM cosmology 
cosmoSCDM = FlatLambdaCDM(H0=70, Om0=0.95, Ob0=0.05)
print 'cosmo', cosmoSCDM, '\n'

# embed()

z = np.logspace(-2,1,1000)
# plt.subplot(121)
plt.loglog(z, cosmoLCDM.comoving_distance(z).value, 'b', label=r'$f_K(\chi)$')
plt.loglog(z, cosmoLCDM.luminosity_distance(z).value, 'r', label=r'$D_L$')
plt.loglog(z, cosmoLCDM.angular_diameter_distance(z).value, 'g', label=r'$D_A$')
plt.loglog(z, cosmoSCDM.comoving_distance(z).value, 'b--')#, label=r'$f_K(\chi)$')
plt.loglog(z, cosmoSCDM.luminosity_distance(z).value, 'r--')#, label=r'$D_L$')
plt.loglog(z, cosmoSCDM.angular_diameter_distance(z).value, 'g--')#, label=r'$D_A$')
plt.text(0.8, 1e4, r'$\Omega_{\Lambda}$', fontsize=22)
plt.text(2.6, 1e4, r'No $\Omega_{\Lambda}$', fontsize=22)
# plt.setp(plt.gca(), 'xticklabels',[0.01, 0.1, 1, 10])
plt.xlabel(r'$z$')
plt.ylabel(r'Distances [Mpc]')
plt.legend(loc='best')
plt.ylim([30,1e5])
# plt.show()
plt.savefig('dists.pdf', bbox_inches='tight')


