"""
Created on 24th of April 2023
Matéo Bédès

Last edited 24th of April 2023
Matéo Bédès

"""

import matplotlib.pyplot as plt
import numpy as np
import random as rand
from random import gauss
import lib_oscil as osc_mod
import params as osc_par
from scipy import interpolate
from scipy.interpolate import interp1d
import distribchi2 as osc_stat 

Nbin=25  #Has to be the same in lib_oscil.py
Nmultiply=30
Reco_nue = 20.2
bias_nue = 3.1
Reco_numu = 16.3
bias_numu = -3.4
Emax=4
Npot=1.1e21
#For 2 parameters fitted
Dchi2_1s = 2.3 
Dchi2_3s = 11.8
Dchi2_5s = 29.8

# Output data are extracted and smoothed
#nu_e

far_detector_spectrum = np.loadtxt('data/output4GeV.txt') #One extracts the output (WITH THE NOISE)	
ene_FD_article, spec_FD_article = np.array([i[0] for i in far_detector_spectrum]), np.array([i[1] for i in far_detector_spectrum])

fd_spectrum_noise = np.loadtxt('data/output_noise4GeV.txt') #Noise on the energy distribution of the Nve
ene_FD_art_noise, spec_FD_art_noise = np.array([i[0] for i in fd_spectrum_noise]), np.array([i[1] for i in fd_spectrum_noise])

list_FD_nue = interp1d(ene_FD_article, spec_FD_article, kind='cubic')
E_lisseND = np.linspace(0.5, Emax,Nbin)
N_FD_hist_nue= list_FD_nue(E_lisseND)
#Output-noise data smoothing
list_FD_noisenue = interp1d(ene_FD_art_noise, spec_FD_art_noise, kind='cubic')
N_FD_hist_noisenue = list_FD_noisenue(E_lisseND)
#(Output-with-noise data) - (Output-noise data) = Output-without-noise data
N_FD_hist_nue = np.array(N_FD_hist_nue) - np.array(N_FD_hist_noisenue)

#nu_mu

far_detector_spectrum_mu = np.loadtxt('data/out_numu4GeV.txt') #One extracts also the output
Emu_FD_article, specmu_FD_article = np.array([i[0] for i in far_detector_spectrum_mu]), np.array([i[1] for i in far_detector_spectrum_mu])*1e03

fd_spectrum_mu_noise = np.loadtxt('data/out_numu_noise4GeV.txt') #Noise on the energy distribution of the Nvmu
ene_FD_art_mu_noise, spec_FD_art_mu_noise = np.array([i[0] for i in fd_spectrum_mu_noise]), np.array([i[1] for i in fd_spectrum_mu_noise])

list_FD_numu = interp1d(Emu_FD_article, specmu_FD_article, kind='cubic') #with noise
list_FD_noisenumu = interp1d(ene_FD_art_mu_noise, spec_FD_art_mu_noise, kind='cubic') #The noise
N_FD_hist_numu= list_FD_numu(E_lisseND)  #with noise
N_FD_hist_noisenumu = list_FD_noisenumu(E_lisseND) #The noise
#(Output-with-noise data) - (Output-noise data) = Output-without-noise data
N_FD_hist_numu = np.array(N_FD_hist_numu) - np.array(N_FD_hist_noisenumu)


#One extracts the data of the flux of nu_mu in the nu_mode

Phi_nu = np.loadtxt('data/Phi_mu_nu_mode_mu4GeV.txt')
Phi_nu_mu, E = np.array([i[0] for i in Phi_nu]), np.array([i[1] for i in Phi_nu])


#Input data are smoothed by interp1D

list_ND_mc = interp1d(E, Phi_nu_mu, kind='cubic')

Phi_lisseND = list_ND_mc(E_lisseND)
Phi_lisseND = np.array(Phi_lisseND)

#Extrapolation factor is smoothed

F_nd_fd = np.loadtxt('data/F_nd_fd.txt')
ene_article, F_article = np.array([i[0] for i in F_nd_fd]), np.array([i[1]*1e-06 for i in F_nd_fd])
F = interp1d(ene_article, F_article, kind='cubic')
F_lisse = F(E_lisseND)

#Cross section data smoothing

sigma_o_E = np.loadtxt('data/sigma-E-2103-04797_4GeV.txt')
E, sigma_o_E_nue_data = np.array([i[0] for i in sigma_o_E]), np.array([i[1] for i in sigma_o_E])*1e-38
sigma_o_E_nue_interp = interp1d(E, sigma_o_E_nue_data, kind='cubic')
sigma_nue = np.array(sigma_o_E_nue_interp(E_lisseND)) * E_lisseND 
M_Ar = 40 #Molar mass of Argon in g/mol
Na = 6e23 #Avogadro cst in mol^-1
m_Ar_FD = 40e9 #masse of argon in the FD
N_nuc_per_Ar = 40 #Number of nucleon per atom of Ar
sigma_FD = m_Ar_FD * Na/M_Ar * N_nuc_per_Ar * sigma_nue

#Chi2 as a function of Dcp and theta23
"""

L_deltaCP_o_pi = np.linspace (-1,1,250) 
L_theta23 = np.linspace(0.6,1.0,250)
L_sin23 = np.sin(L_theta23)**2
L2_Int = []
LL2_N_nue = []
LL2_N_numu = []
L2_khi2 = []


for i in range (0,len(L_theta23)):
	

	LL2_N_nue = LL2_N_nue +[[]]
	LL2_N_numu = LL2_N_numu + [[]]
	L2_khi2 = L2_khi2 + [[]]
	L2_Int = L2_Int +[[]]

	for j in range (0,len(L_deltaCP_o_pi)) :
		LL2_N_nue[i] = LL2_N_nue[i] + [np.array(osc_mod.Rec2( Reco_nue,E_lisseND, Phi_lisseND * osc_mod.probability_oscillation_v2(osc_par.delta_m31, osc_par.delta_m32 ,osc_par.theta12 , osc_par.theta13, L_theta23[i], E_lisseND, 1285 ,L_deltaCP_o_pi[j]*np.pi*180/np.pi, True, True) * F_lisse * sigma_FD * Npot *Nmultiply , bias_nue))/Nmultiply]
		
		LL2_N_numu[i] = LL2_N_numu[i] + [np.array(osc_mod.Rec2_numu(Reco_numu, E_lisseND, Phi_lisseND * osc_mod.proba_disapp_0_t(osc_par.delta_m32, L_theta23[i], E_lisseND, 1285 , True) * F_lisse * sigma_FD * Npot *Nmultiply , bias_numu))/Nmultiply]
		
		
		L2_khi2[i] = L2_khi2[i] + [0]
		for k in range (0,len(E_lisseND)):
			if LL2_N_nue[i][j][k] != 0 and LL2_N_numu[i][j][k] != 0:
				L2_khi2[i][j] = L2_khi2[i][j] + ( LL2_N_nue[i][j][k] - N_FD_hist_nue[k] )**2/(4*LL2_N_nue[i][j][k]) + + ( LL2_N_numu[i][j][k] - N_FD_hist_numu[k] )**2/(4*LL2_N_numu[i][j][k]) 
			else :
				L2_khi2[i][j] = L2_khi2[i][j]

#Some statistics value

chi2min = round( L2_khi2[osc_mod.ind_min(L2_khi2)[0]][osc_mod.ind_min(L2_khi2)[1]] , 1)
chi2min_dof_str = str(round( L2_khi2[osc_mod.ind_min(L2_khi2)[0]][osc_mod.ind_min(L2_khi2)[1]] / (Nbin-2) , 2))  #Chi2min over the Number of degrees of freedom
pvalue_3sigma = osc_stat.p_value( chi2min + Dchi2_3s , Nbin - 2, 0,80) #p_value at 3sigma
pvalue_3sigma = round ( pvalue_3sigma , 3 )
p_3sigma_str = str(pvalue_3sigma) 


# 2D histogram with pcolormesh

import matplotlib.colors as colors

#Create the legend
legend = [ plt.scatter(0, 0,color = 'green', marker = '_', s = 40) , plt.scatter(0, 0,color = 'pink', marker = '_', s = 40) ,  plt.scatter(0, 0, color = 'grey', marker = '_', s = 40), plt.scatter(0, 0,color = 'pink', marker = '_', s = 40)  , plt.scatter(0, 0,color = 'red', marker = 'o', s = 40), plt.scatter(0, 0,color = 'red', marker = 'x', s = 40)]
legend_label = ['$1\sigma$', '$3\sigma$', '$5\sigma$', (f"$p(3\sigma)$={p_3sigma_str}") ,'true point',(f"$\chi^2/Ndof$={chi2min_dof_str}")]
plt.show()
plt.close()

#Plot
DCPoPI , STheta23  = np.meshgrid (L_deltaCP_o_pi , L_sin23)  #Make the grid
fig, ax = plt.subplots()
origin = 'lower'
CS2 = ax.contour(DCPoPI, STheta23 , L2_khi2, levels=[chi2min + Dchi2_1s , chi2min + Dchi2_3s ,chi2min + Dchi2_5s], colors=['green', 'pink', 'grey'],norm='log', origin=origin) #Make the contour at 1sigma , 3sigma and 5 sgima 
pc = ax.pcolormesh( DCPoPI, STheta23 , L2_khi2, norm='log')
fig.colorbar(pc,ax=ax,label='$\chi^2_{combined}$')
plt.title('$\epsilon = 50 kT.MW.yrs$')
plt.xlabel('$\delta_{CP} / \pi $ ')
plt.ylabel('$sin(\Theta_{23})^2$ ')
plt.plot(0 ,np.sin(0.86)**2,'o',color='red',label='true point')  #Plot the true point
plt.plot( L_deltaCP_o_pi[osc_mod.ind_min(L2_khi2)[1]], L_sin23[osc_mod.ind_min(L2_khi2)[0]] ,'x',color='red',label=(f"$\chi^2/Ndof$={chi2min_dof_str}"))
plt.legend(legend, legend_label, loc ='lower right')
plt.show()
"""

#Chi2 as a function of Stheta23 and Dm23				

L_delta_m32 = np.linspace (2.1e-03,2.9e-03,250) 
L_theta23 = np.linspace(0.4,1.1,250)
L_sin23 = np.sin(L_theta23)**2
L2_Int = []
LL2_N_nue = []
LL2_N_numu = []
L2_khi2 = []


for i in range (0,len(L_theta23)):
	

	LL2_N_nue = LL2_N_nue +[[]]
	LL2_N_numu = LL2_N_numu + [[]]
	L2_khi2 = L2_khi2 + [[]]
	L2_Int = L2_Int +[[]]

	for j in range (0,len(L_delta_m32)) :
		LL2_N_nue[i] = LL2_N_nue[i] + [np.array(osc_mod.Rec2( Reco_nue,E_lisseND, Phi_lisseND * osc_mod.probability_oscillation_v2(osc_par.delta_m31, osc_par.delta_m32 ,osc_par.theta12 , osc_par.theta13, L_theta23[i], E_lisseND, 1285 ,0*np.pi*180/np.pi, True, True) * F_lisse * sigma_FD * Npot *Nmultiply , bias_nue))/Nmultiply]
		
		LL2_N_numu[i] = LL2_N_numu[i] + [np.array(osc_mod.Rec2_numu(Reco_numu, E_lisseND, Phi_lisseND * osc_mod.proba_disapp_0_t(L_delta_m32[j], L_theta23[i], E_lisseND, 1285 , True) * F_lisse * sigma_FD * Npot *Nmultiply , bias_numu))/Nmultiply]
		
		
		L2_khi2[i] = L2_khi2[i] + [0]
		for k in range (0,len(E_lisseND)):
			if LL2_N_nue[i][j][k] != 0 and LL2_N_numu[i][j][k] != 0:
				L2_khi2[i][j] = L2_khi2[i][j] + ( LL2_N_nue[i][j][k] - N_FD_hist_nue[k] )**2/(4*LL2_N_nue[i][j][k]) + + ( LL2_N_numu[i][j][k] - N_FD_hist_numu[k] )**2/(4*LL2_N_numu[i][j][k]) 
			else :
				L2_khi2[i][j] = L2_khi2[i][j]

# Some statistics values 
chi2min = round( L2_khi2[osc_mod.ind_min(L2_khi2)[0]][osc_mod.ind_min(L2_khi2)[1]] , 1)
chi2min_dof_str = str(round( L2_khi2[osc_mod.ind_min(L2_khi2)[0]][osc_mod.ind_min(L2_khi2)[1]] / (Nbin-2) , 2))  #Chi2min over the Number of degrees of freedom
pvalue_3sigma = osc_stat.p_value( chi2min + Dchi2_3s , Nbin - 2, 0,80) #p_value at 3sigma
pvalue_3sigma = round ( pvalue_3sigma , 3 )
p_3sigma_str = str(pvalue_3sigma)

#Create the legend
legend = [ plt.scatter(0, 0,color = 'green', marker = '_', s = 40) , plt.scatter(0, 0,color = 'pink', marker = '_', s = 40) ,  plt.scatter(0, 0, color = 'grey', marker = '_', s = 40), plt.scatter(0, 0,color = 'pink', marker = '_', s = 40)  , plt.scatter(0, 0,color = 'red', marker = 'o', s = 40), plt.scatter(0, 0,color = 'red', marker = 'x', s = 40)]
legend_label = ['$1\sigma$', '$3\sigma$', '$5\sigma$', (f"$p(3\sigma)$={p_3sigma_str}") ,'true point',(f"$\chi^2/Ndof$={chi2min_dof_str}")]
plt.show()
plt.close()


# 2D histogram with pcolormesh
import matplotlib.colors as colors
DM2_23 , STheta23 = np.meshgrid ( L_delta_m32*1000 , L_sin23 )
fig, ax = plt.subplots()
origin = 'lower'
CS2 = ax.contour(DM2_23 , STheta23 , L2_khi2, levels=[chi2min + Dchi2_1s , chi2min + Dchi2_3s ,chi2min + Dchi2_5s], colors=['green', 'pink', 'grey'],norm='log', origin=origin) #Make the contour at 1sigma , 3sigma and 5 sgima 

pc = ax.pcolormesh(DM2_23 , STheta23, L2_khi2 , norm='log')
fig.colorbar(pc,ax=ax,label='$\chi^2_{combined}$')
plt.xlabel('$\Delta$$m^2_{23}$ ')
plt.ylabel('$sin(\Theta_{23})^2$ ')
plt.scatter(2.45*pow(10,-3)*1000, np.sin(np.radians(49.8))**2 , s=10, color='red' )
plt.title('$\epsilon = 50 kT.MW.yrs$')
plt.plot(L_delta_m32[osc_mod.ind_min(L2_khi2)[1]] *1000, L_sin23[osc_mod.ind_min(L2_khi2)[0]] ,'x',color='red',label=(f"$\chi^2/Ndof$={chi2min_dof_str}"))
plt.legend(legend, legend_label, loc ='lower right')
plt.show()

