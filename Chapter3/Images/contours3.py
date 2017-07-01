import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
from matplotlib.ticker import NullFormatter
from matplotlib.ticker import MaxNLocator
from getdist import plots, MCSamples
from IPython import embed

# Matplotlib defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rc('text',usetex=True)
rc('font',**{'family':'serif','serif':['Computer Modern']})
plt.rcParams['axes.linewidth']  = 3.
plt.rcParams['axes.labelsize']  = 30
plt.rcParams['xtick.labelsize'] = 22
plt.rcParams['ytick.labelsize'] = 22
plt.rcParams['xtick.major.size'] = 10
plt.rcParams['ytick.major.size'] = 10
plt.rcParams['xtick.minor.size'] = 4
plt.rcParams['ytick.minor.size'] = 4
plt.rcParams['legend.fontsize']  = 22
plt.rcParams['legend.frameon']  = False

plt.rcParams['xtick.major.width'] = 1
plt.rcParams['ytick.major.width'] = 1
plt.rcParams['xtick.minor.width'] = 1
plt.rcParams['ytick.minor.width'] = 1

# names  = ['b', 'A']
# labels = names

# zbins     = ['1.5_10', '1.5_2.1', '2.1_10']
# zbins_leg = {'1.5_10':'$z > 1.5$', '1.5_2.1':'$1.5 < z < 2.1$', '2.1_10':'$z > 2.1$'}

# samps = {}
# for zbin in zbins:
# 	samps[zbin] = np.loadtxt('b_A_FINAL_SMM/ALL/samples/samples_kggg_'+zbin+'.dat', unpack=False)

# MCsamps = {}
# for zbin in zbins:
# 	MCsamps[zbin] = MCSamples(samples=samps[zbin], names=names, labels=labels)

# g = plots.getSinglePlotter()
# g.settings.figure_legend_frame = False
# g.settings.legend_frame        = False
# g.settings.legend_fontsize     = 20
# g.settings.lab_fontsize        = 30
# g.settings.axes_fontsize       = 22
# g.plot_2d([MCsamps['1.5_10'], MCsamps['1.5_2.1'], MCsamps['2.1_10']], 'b', 'A', filled=True, colors=['red', 'royalblue', 'green'], alphas=[0.6,0.6,0.6])
# plt.plot([3.54],[1.45], 'x', color='darkred', lw=2)
# plt.plot([2.89],[1.48], 'x', color='darkblue', lw=2)
# plt.plot([4.75],[1.37], 'x', color='darkgreen', lw=2)
# g.add_y_marker(1, ls='--', lw=2, color='k')
# g.add_legend(zbins_leg.values(), legend_loc='upper right')
# g.export(fname='b_A_contours_getdist_SMM_3sigma.pdf')
# plt.show()


# # ~~~~~~~~~~~~~~~~~ A - A_bias ~~~~~~~~~~~~~~~~~~~~~~~~~~~

# names  = ['A_bias', 'A']
# labels = ['$\mathcal{A}_{bias}$', 'A']

# zbins     = ['1.5_10', '1.5_2.1', '2.1_10']
# zbins_leg = {'1.5_10':'$z > 1.5$', '1.5_2.1':'$1.5 < z < 2.1$', '2.1_10':'$z > 2.1$'}

# samps = {}
# for zbin in zbins:
# 	samps[zbin] = np.loadtxt('A_A_bias_FINAL_SMM/ALL/samples/samples_kggg_'+zbin+'.dat', unpack=False)

# MCsamps = {}
# for zbin in zbins:
# 	MCsamps[zbin] = MCSamples(samples=samps[zbin], names=names, labels=labels)

# g = plots.getSinglePlotter()
# g.settings.figure_legend_frame = False
# g.settings.legend_frame        = False
# g.settings.legend_fontsize     = 20
# g.settings.lab_fontsize        = 30
# g.settings.axes_fontsize       = 22

# g.plot_2d([MCsamps['2.1_10'], MCsamps['1.5_2.1'], MCsamps['1.5_10']], 'A_bias', 'A', filled=True,  
# 				colors=['green', 'royalblue','red'], alphas=[0.6,0.6,0.6], lims=[.35, 1.2, .7, 3.])
# g.add_y_marker(1, ls='--', lw=2, color='k')
# g.add_x_marker(1, ls='--', lw=2, color='k')
# plt.plot([0.82],[1.49], 'x', color='darkred', lw=2)
# plt.plot([0.77],[1.51], 'x', color='darkblue', lw=2)
# plt.plot([1.02],[1.43], 'x', color='darkgreen',lw=2)
# g.add_legend(legend_labels=[zbins_leg['2.1_10'], zbins_leg['1.5_2.1'], zbins_leg['1.5_10']], legend_loc='upper left')
# g.finish_plot( label_order='-1')
# g.export(fname='A_Abias_contours_getdist_SMM_3sigma.pdf')
# plt.show()

# ~~~~~~~~~~~~~~~~~ A - b SMM 3 vs 5 sigma ~~~~~~~~~~~~~~~~~~~~~~~~~~~

names  = ['b', 'A']
labels = names

samps_gg   = np.loadtxt('b_A_FAKE_nomag/ALL/samples/samples_gg_1.5_10.dat', unpack=False)
samps_kg   = np.loadtxt('b_A_FAKE_nomag/ALL/samples/samples_kg_1.5_10.dat', unpack=False)
samps_kggg = np.loadtxt('b_A_FAKE_nomag/ALL/samples/samples_kggg_1.5_10.dat', unpack=False)

MCsamps_gg   = MCSamples(samples=samps_gg, names=['b'], labels=['b'])
MCsamps_kg   = MCSamples(samples=samps_kg, names=['b'], labels=['b'])
MCsamps_kggg = MCSamples(samples=samps_kggg, names=names, labels=labels)
MCsamps_kggg_dummy = MCSamples(samples=np.zeros((100,2))+np.random.rand(200).reshape((100,2)), names=names, labels=labels)

g = plots.getSinglePlotter()
g.settings.figure_legend_frame = False
g.settings.legend_frame        = False
g.settings.legend_fontsize     = 20
g.settings.lab_fontsize        = 30
g.settings.axes_fontsize       = 22
g.plot_2d([MCsamps_kggg,MCsamps_kggg_dummy,MCsamps_kggg_dummy], 'b', 'A', filled=True, colors=['royalblue','tomato','green'], alphas=[0.6], lims=[2.6,3.2,0.75,1.35])
plt.errorbar([3], [1.01], xerr=[0.05], fmt='o', color='tomato', elinewidth=2, label=r'$gg$')
plt.errorbar([2.92], [.99], xerr=[0.27], fmt='x', color='green', elinewidth=2, label=r'$\kappa g$')
g.add_y_marker(1, ls='--', lw=2, color='k')
g.add_x_marker(3, ls='--', lw=2, color='k')
g.add_legend(['$\kappa g + gg$', '$gg$', '$\kappa g$'], legend_loc='best')
# plt.legend()
g.export(fname='b_A_correct.pdf')
# plt.show()

# ~~~~~~~~~~~~~~~~~ A - b SMM vs PEARSON (3sigma) ~~~~~~~~~~~~~~~~~~~~~~~~~~~

# names  = ['b', 'A']
# labels = names

# zbins     = ['1.5_10', '1.5_2.1', '2.1_10']
# zbins_leg = {'1.5_10':'$z > 1.5$', '1.5_2.1':'$1.5 < z < 2.1$', '2.1_10':'$z > 2.1$'}

# sampsSMM  = {}
# sampsPaer = {}
# for zbin in zbins:
# 	sampsSMM[zbin]  = np.loadtxt('b_A_FINAL_SMM/ALL/samples/samples_kggg_'+zbin+'.dat', unpack=False)
# 	sampsPaer[zbin] = np.loadtxt('b_A_FINAL_Pearson/ALL/samples/samples_kggg_'+zbin+'.dat', unpack=False)

# MCsampsSMM = {}
# MCsampsPaer = {}
# for zbin in zbins:
# 	MCsampsSMM[zbin] = MCSamples(samples=sampsSMM[zbin], names=names, labels=labels)
# 	MCsampsPaer[zbin] = MCSamples(samples=sampsPaer[zbin], names=names, labels=labels)

# g = plots.getSinglePlotter()
# g.settings.figure_legend_frame = False
# g.settings.legend_frame        = False
# g.settings.legend_fontsize     = 20
# g.settings.lab_fontsize        = 30
# g.settings.axes_fontsize       = 22
# g.plot_2d([MCsampsPaer['1.5_10'], MCsampsPaer['1.5_2.1'], MCsampsPaer['2.1_10']], 'b', 'A', filled=True, colors=['red', 'royalblue', 'green'], alphas=[0.6,0.6,0.6], lims=[2.1,7,0.6,2.4])
# g.plot_2d([MCsampsSMM['1.5_10'], MCsampsSMM['1.5_2.1'], MCsampsSMM['2.1_10']], 'b', 'A', colors=['red', 'royalblue', 'green'], line_args=[{'ls':'--'},{'ls':'--'},{'ls':'--'}], lims=[2.1,7,0.6,2.4])
# # plt.plot([3.54],[1.45], 'x', color='darkred', lw=2)
# # plt.plot([2.89],[1.48], 'x', color='darkblue', lw=2)
# # plt.plot([4.75],[1.37], 'x', color='darkgreen', lw=2)
# g.add_y_marker(1, ls='--', lw=2, color='k')
# g.add_legend(zbins_leg.values(), legend_loc='upper right')
# g.export(fname='b_A_contours_getdist_SMM_Pearson.pdf')
# # plt.show()

