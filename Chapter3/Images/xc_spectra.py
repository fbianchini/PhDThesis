import numpy as np
import mycosmo
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import time
import seaborn as sns
from matplotlib import rc, rcParams
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from matplotlib.ticker import NullFormatter
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FormatStrFormatter

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

cosmo = {'H0':67.94,
         'omegab':0.0480,
         'omegac':0.24225,#0.255,
         'omegak':0.,
         'omegav':0.683,         
         'omegan':0.,
         'TCMB':2.7255,
         'yhe':0.24,
         'Num_Nu_massless':3.046,
         'Num_Nu_massive':0.,
         'scalar_index':0.9624,
         'reion__redshift':11.42,
         'reion__optical_depth':0.0943,
         'scalar_index':0.9582,
         'scalar_amp':2.2107e-9,
         'scalar_running':0,
         'reion__use_optical_depth':1,
         'DoLensing':0
         }
print "...Evaluating cosmo module..."
pf = mycosmo.lcdm(cosmo)
print "...done..."

# k = np.logspace(-3,1,100)
# plt.loglog(k, pf.pkz.P(0,k))
# plt.loglog(k, pf.pkz.P(0.4,k))
# plt.loglog(k, pf.pkz.P(1,k))
# plt.show()

# z_ph = np.loadtxt('hatlas_SGP_psf_zsissa_new.dat', usecols=[13], skiprows=3, unpack=True)
# print 'start'
# dndz_005, z = mycosmo.catalog2dndz(z_ph, zbins=(1.5, 2.1), sigma=0.05)
# print 'first done'
# dndz_010, z = mycosmo.catalog2dndz(z_ph, zbins=(1.5, 2.1), sigma=0.10)
# print 'second done'
# dndz_020, z = mycosmo.catalog2dndz(z_ph, zbins=(1.5, 2.1), sigma=0.20)

# plt.plot(z, dndz_020, label='0.20')
# plt.plot(z, dndz_010, label='0.10')
# plt.plot(z, dndz_005, label='0.05')
# plt.legend()
# plt.show()

# print "...Reading and interpolating dN/dz..."
# z_herschel_low, dndz_herschel_low   = np.loadtxt('hatlas_ALL_convolved_dndz_1.5_2.1.dat', unpack=True)
# z_herschel_high, dndz_herschel_high = np.loadtxt('dndz_GAMA_conv_only35_FIRbased_all_0.8_1.5_5.dat', unpack=True)
# dndz_low  = interp1d(z_herschel_low, dndz_herschel_low, bounds_error=False, fill_value=0.)
# dndz_high = interp1d(z_herschel_high, dndz_herschel_high, bounds_error=False, fill_value=0.)
# print "...done..."



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# cl_gg   = pf.power_spectra_gg(dndz_high)
# cl_gmu  = pf.power_spectra_gmu(dndz_high)
# cl_mumu = pf.power_spectra_mumu(dndz_high)
# z_andrea, b_andrea_ = np.loadtxt('/var/folders/SS/SSq5iu+DExywCYiAN6Q3FU+++TI/-Tmp-/com.apple.mail.drag-T0x1005200a0.tmp.NPhdP6/bias_z_fede.dat', unpack=True)
# b_andrea = interp1d(z_andrea, b_andrea_, bounds_error=False, fill_value=0.)

z_test = np.linspace(0,5,1000)

# plt.plot(z_test, b_andrea(z_test))
# plt.show()

b_0 = lambda z: 1
b_1 = lambda z: (1+z)**1
b_2 = lambda z: (1+z)**2
# b_a = lambda z: b_andrea(z)

# plt.plot(z_test, pf.w_kappa_cmb(z_test)*pf.w_cl(z_test, dndz_high, b=b_1), label=r'$W^{\kappa}\frac{dN}{dz}\, b=1$')
# plt.plot(z_test, pf.w_cl(z_test, dndz_high, b=b_1)**2, label=r'$(\frac{dN}{dz})^2\, b=1$')
# plt.plot(z_test, pf.w_kappa_cmb(z_test)*pf.w_cl(z_test, dndz_high, b=b_andrea), label=r'$W^{\kappa}\frac{dN}{dz}$ Andrea bias')
# plt.plot(z_test, pf.w_cl(z_test, dndz_high, b=b_andrea)**2, label=r'$(\frac{dN}{dz})^2$  Andrea bias')
# plt.legend()
# plt.xlabel(r'$z$')
# plt.show()

def n_z(z,z0=1,sigma=0.3): 
   return 1./np.sqrt(2.*np.pi*sigma**2) * np.exp(-(z-z0)**2/2./sigma**2)
   # return lambda z: 1./np.sqrt(2.*np.pi*sigma**2) * np.exp(-(z-z0)**2/2./sigma**2)

# plt.plot(z_test, n_z(z_test,))
# plt.show()

cl_kg_tot_b_0 = pf.power_spectra_kg_tot(n_z, b=b_0, alpha=1)*1e7
cl_kg_tot_b_1 = pf.power_spectra_kg_tot(n_z, b=b_1, alpha=1)*1e7
cl_kg_tot_b_2 = pf.power_spectra_kg_tot(n_z, b=b_2, alpha=1)*1e7
# cl_kg_tot_b_a = pf.power_spectra_kg_tot(dndz_high, b=b_a, alpha=1)

cl_gg_tot_b_0 = pf.power_spectra_gg_tot(n_z, b=b_0, alpha=1)*1e5
cl_gg_tot_b_1 = pf.power_spectra_gg_tot(n_z, b=b_1, alpha=1)*1e5
cl_gg_tot_b_2 = pf.power_spectra_gg_tot(n_z, b=b_2, alpha=1)*1e5
# cl_gg_tot_b_1 = pf.power_spectra_gg_tot(dndz_high, b=b_a, alpha=1)


f, (ax1,ax2) = plt.subplots(2, sharex=True)#, figsize=(10,8))

# plt.subplot(1,2,1)
ax1.plot(cl_kg_tot_b_0, label=r'$b=1$')
ax1.plot(cl_kg_tot_b_1, label=r'$b=(1+z)$')
ax1.plot(cl_kg_tot_b_2, label=r'$b=(1+z)^2$')
# plt.plot(cl_kg_tot_b_z, label=r'Andrea bias')
ax1.legend(loc='best')
ax1.yaxis.set_major_locator(MaxNLocator(nbins=5))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax1.set_ylabel(r'$C_{\ell}^{\kappa g} (\times 10^{-7})$')

# ax2.axhline(ls='--', color='k')
ax2.plot(cl_kg_tot_b_1/cl_kg_tot_b_0-1, label=r'$b=(1+z)$', color='g')
ax2.plot(cl_kg_tot_b_2/cl_kg_tot_b_0-1, label=r'$b=(1+z)^2$', color='b')
ax2.set_xlabel(r'$\ell$')
ax2.set_ylabel(r'$\Delta C_{\ell}^{\kappa g}/C_{\ell}^{\kappa g}$')
ax2.legend(loc='best')
ax2.set_xlim([2,600])
ax2.yaxis.set_major_locator(MaxNLocator(nbins=4, prune='upper'))
f.subplots_adjust(hspace=0)
plt.savefig('cl_kg_bz.pdf', bbox_inches='tight')


f, (ax1,ax2) = plt.subplots(2, sharex=True)#, figsize=(10,8))

# plt.subplot(1,2,1)
ax1.plot(cl_gg_tot_b_0, label=r'$b=1$')
ax1.plot(cl_gg_tot_b_1, label=r'$b=(1+z)$')
ax1.plot(cl_gg_tot_b_2, label=r'$b=(1+z)^2$')
# plt.plot(cl_kg_tot_b_z, label=r'Andrea bias')
ax1.legend(loc='best')
ax1.yaxis.set_major_locator(MaxNLocator(nbins=5))
# plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax1.set_ylabel(r'$C_{\ell}^{gg} (\times 10^{-5})$')

# ax2.axhline(ls='--', color='k')
ax2.plot(cl_gg_tot_b_1/cl_gg_tot_b_0-1, label=r'$b=(1+z)$', color='g')
ax2.plot(cl_gg_tot_b_2/cl_gg_tot_b_0-1, label=r'$b=(1+z)^2$', color='b')
ax2.set_xlabel(r'$\ell$')
ax2.set_ylabel(r'$\Delta C_{\ell}^{gg}/C_{\ell}^{gg}$')
ax2.legend(loc='best')
ax2.set_xlim([2,600])
ax2.yaxis.set_major_locator(MaxNLocator(nbins=4, prune='upper'))
f.subplots_adjust(hspace=0)
plt.savefig('cl_gg_bz.pdf', bbox_inches='tight')

f

#~~~~~~~~~~~~~~~~~
cl_kg  = pf.power_spectra_kg(dndz_high)
cl_kmu = pf.power_spectra_kmu(dndz_high)

l, kg_, kmu_ = np.loadtxt('KG_hatlas_ALL_convolved_z2.1_10_mcmc_magbias_alpha_PLANCK1+highL+WP+Lensing_params.dat', unpack=1)
l, gg_, gmu_, mumu_= np.loadtxt('GG_hatlas_ALL_convolved_z2.1_10_mcmc_magbias_alpha_PLANCK1+highL+WP+Lensing_params.dat', unpack=1)

# cl_gg_tot_nomag = pf.power_spectra_gg_tot(dndz_high, b=2, alpha=1)
# cl_gg_tot_mag   = pf.power_spectra_gg_tot(dndz_high, b=2, alpha=3)

cl_kg_tot_nomag = pf.power_spectra_kg_tot(dndz_high, b=2, alpha=1)
cl_kg_tot_mag   = pf.power_spectra_kg_tot(dndz_high, b=2, alpha=3)

# plt.plot(2**2 * cl_gg + 2.*2*(3-1)*cl_gmu + (3-1)**2 * cl_mumu, '--', label='Old b=2 a=3')
# plt.plot(2**2 * cl_gg, '--', label='Old b=2 a=1')
# plt.plot(cl_gg_tot_mag, label='New w/ mag')
# plt.plot(cl_gg_tot_nomag, label='New w/o mag')
# plt.legend()
# plt.title(r'$C_{\ell}^{gg}$')
# plt.show()

plt.plot(2*cl_kg + (3-1)*cl_kmu, '--', label='Old b=2 a=3')
# plt.plot(2*kg_,'--', label='Old b=2 a=1')
plt.plot(cl_kg_tot_mag, label='New w/ mag')
# plt.plot(cl_kg_tot_nomag, label='New w/o mag')
plt.legend()
plt.title(r'$C_{\ell}^{\kappa g}$')
plt.show()

# plt.plot(cl_gg, label='gg')
# plt.plot(l, gg_, label='gg old')
# plt.plot(cl_gmu, label='gmu')
# plt.plot(l, gmu_, label='gmu old')
# plt.plot(cl_mumu, label='mumu')
# plt.plot(l, mumu_, label='mumu old')
# plt.legend()
# plt.show()

# plt.plot(cl_kg, label='kg')
# plt.plot(l, kg_, label='kg other')
# plt.plot(l, kmu_, label='kmu other')
# plt.plot(cl_kmu, label='kmu')
# plt.legend()
# plt.show()

# plt.show()