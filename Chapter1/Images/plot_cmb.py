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

# Evaluating CMB T&P spectra
params = {
    'output': 'tCl',
    'l_max_scalars': 3000,
    # 'lensing': 'yes',
    'A_s': 2.3e-9,
    'n_s': 0.9624, 
    'h': 0.6711,
    'omega_b': 0.022068,
    'omega_cdm': 0.12029,
    'Omega_k':0.}

cosmo = Class()

# embed()

fig = plt.figure(figsize=(18,5))
ax1 = fig.add_subplot(121)

for source in ['tsw', 'eisw', 'lisw', 'dop', '']:
	if source == 'tsw':  
		lab = 'SW'
		col = sns.xkcd_rgb['marigold']
	if source == 'eisw': 
		lab = 'e-iSW'
		col = sns.xkcd_rgb['turquoise']
	if source == 'lisw': 
		lab = 'l-iSW'
		col = sns.xkcd_rgb['burnt orange']
	if source == 'dop':
		lab = 'Doppler'
		col = sns.xkcd_rgb['greenish grey']
	
	params['temperature contributions'] = source 
	
	if source == '':
		lab = 'Total'
		col = sns.xkcd_rgb['black']
		del(params['temperature contributions'])

	print source, lab, col

	cosmo.set(params)
	cosmo.compute()

	cls = cosmo.raw_cl(2500)
	l   = np.arange(cls['tt'].shape[0])
	
	ax1.plot(l*(l+1)*cls['tt']/2/np.pi*(cosmo.T_cmb()*1e6)**2, color=col, label=lab)
	# plt.plot(cls['tt'], label=lab)

	cosmo.struct_cleanup()
	cosmo.empty()

ax1.legend(loc='best')
ax1.set_xlabel(r'$\ell$')
ax1.set_ylabel(r'$\mathcal{D}_{\ell}^{TT} [\mu K]^2$')
ax1.set_xscale('log')
ax1.set_xlim([2,2500])

ax2 = fig.add_subplot(122)

params = {
    'output': 'tCl pCl',
    'l_max_scalars': 3000,
    'lensing': 'no',
    'A_s': 2.3e-9,
    'n_s': 0.9624, 
    'h': 0.6711,
    'omega_b': 0.022068,
    'omega_cdm': 0.12029,
    'Omega_k':0.,
    'modes':'s t',
    'z_reio':0.,
    'r':0.2}

cosmo = Class()
cosmo.set(params)
cosmo.compute()
cls  = cosmo.raw_cl(2500)
l    = np.arange(cls['tt'].shape[0])

cosmo.struct_cleanup()
cosmo.empty()
params['z_reio']=8
cosmo.set(params)
cosmo.compute()
clsr  = cosmo.raw_cl(2500)

ax2.loglog(l*(l+1)*clsr['tt']/2/np.pi*(cosmo.T_cmb()*1e6)**2, 'b', label='TT')
ax2.loglog(l*(l+1)*cls['tt']/2/np.pi*(cosmo.T_cmb()*1e6)**2, 'b--')
# ax2.loglog(l*(l+1)*cls['te']/2/np.pi*(cosmo.T_cmb()*1e6)**2, 'r', label='TE')
# ax2.loglog(-l[cls['te']<0]*(l[cls['te']<0]+1)*cls['te'][cls['te']<0]/2/np.pi*(cosmo.T_cmb()*1e6)**2, 'r--')
ax2.loglog(l*(l+1)*clsr['ee']/2/np.pi*(cosmo.T_cmb()*1e6)**2,'g', label='EE')
ax2.loglog(l*(l+1)*cls['ee']/2/np.pi*(cosmo.T_cmb()*1e6)**2, 'g--')
ax2.loglog(l*(l+1)*clsr['bb']/2/np.pi*(cosmo.T_cmb()*1e6)**2, 'r',label='BB')
ax2.loglog(l*(l+1)*cls['bb']/2/np.pi*(cosmo.T_cmb()*1e6)**2, 'r--')
ax2.set_xlabel(r'$\ell$')
ax2.set_ylabel(r'$\mathcal{D}_{\ell} [\mu K]^2$')
ax2.legend(loc='best')

plt.savefig('cmb.pdf', bbox_inches='tight')







