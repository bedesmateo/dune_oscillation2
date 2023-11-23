"""
Created on 6th of april 2023
Matéo Bédès

Last edited 6th of april 2023
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

#                           REFERENCES

#Input data source : Article 2002.03005.pdf
#Output data source : Article 2109.01304v1.pdf
#Extrapolation factor source : Article  2006.16043.pdf
#Cross section source :  Article  2103.04797.pdf

Nbin=25  #Has to be the same in lib_oscil.py
Nmultiply=35
Emax = 4 #Gev has to be the same in lib_oscil.py
Reco = 15.9
bias= -2.38

#   I )                      INPUT DATA SMOOTHING

#The input data are extracted and they are smoothed. The output data are also extracted



far_detector_spectrum_mu = np.loadtxt('data/out_numu.txt') #One extracts also the output
Emu_FD_article, specmu_FD_article = np.array([i[0] for i in far_detector_spectrum_mu]), np.array([i[1] for i in far_detector_spectrum_mu])*1e03

fd_spectrum_mu_noise = np.loadtxt('data/out_numu_noise.txt') #Noise on the energy distribution of the Nvmu
ene_FD_art_mu_noise, spec_FD_art_mu_noise = np.array([i[0] for i in fd_spectrum_mu_noise]), np.array([i[1] for i in fd_spectrum_mu_noise])



#One extracts the data of the flux of nu_mu in the nu_mode

Phi_nu = np.loadtxt('data/Phi_mu_nu_mode_mu.txt')
Phi_nu_mu, E = np.array([i[0] for i in Phi_nu]) , np.array([i[1] for i in Phi_nu])

plt.plot(E, Phi_nu_mu,'x', color = 'r',label='$\phi_{\u03BD_{\mu}}^{ND}$') #One plot the data of the nu_mu flux
plt.title('Input data smoothing')

plt.yscale('log')
plt.xlabel('E GeV')
plt.xlim(0.25,Emax)
plt.ylabel('$\phi_{\u03BD_{\mu}}^{ND}$    ${\u03BD_{\mu}}$/$cm^{2}$/POT/1 GeV')
plt.legend()

#Input data are smoothed by interp1D

list_ND_mc = interp1d(E, Phi_nu_mu, kind='cubic')
#One stocks the smoothed input data
E_lisseND = np.linspace(0.5, Emax,Nbin)
Phi_lisseND = list_ND_mc(E_lisseND)
Phi_lisseND = np.array(Phi_lisseND)

plt.bar(E_lisseND,Phi_lisseND,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label='$\phi_{\u03BD_{\mu}}^{ND}$')
plt.yscale('log')
plt.legend()
plt.show()
plt.close()

#  II )                    DISAPP PROBABILITY x SMOOTHED INPUT DATA WITHOUT THE NEAR TO FAR FLUX EXTRAPOLATION


#If cos(theta13) ~= 1
theta23 = np.radians(49.8)
#theta23 = np.radians(51.6)

plt.subplot(311)
plt.title('Effect of P$_{\u03BD_{\mu} -> \u03BD_{\mu}}^0$')

plt.bar(E_lisseND,Phi_lisseND,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label='$\phi_{\u03BD_{\mu}}^{ND}$') #Plot Input data
plt.xlabel('E GeV')
plt.ylabel('$\phi_{\u03BD_{\mu}}^{ND}$ ${\u03BD_{\mu}}$/$cm^{2}$/POT/1 GeV')
plt.legend()

Ex = np.linspace(0.5,Emax,1000)

plt.subplot(312)

plt.plot( Ex, osc_mod.proba_disapp_0 ( theta23, Ex , 1285 , True ) ,label='P$_{\u03BD_{\mu} -> \u03BD_{\mu}}^0$')
plt.xlabel('E GeV')
plt.ylabel('P$_{\u03BD_{\mu} -> \u03BD_{\mu}}^0$')
plt.legend()

plt.subplot(313)


Phi_nu_mu_FD0 = Phi_lisseND * osc_mod.proba_disapp_0 ( theta23, E_lisseND , 1285 , True ) 
plt.bar(E_lisseND, Phi_nu_mu_FD0 ,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label='$\phi_{\u03BD_{\mu}}^{ND}$ * P$_{\u03BD_{\mu}->\u03BD_{\mu}}^0$')

plt.xlabel('E GeV')
plt.ylabel('$\phi_{\u03BD_{\mu}}^{ND}$xP$_{\u03BD_{\mu}->\u03BD_{\mu}}^0$ ${\u03BD_{\mu}}$/$cm^{2}$/POT/1 GeV')


plt.legend()

plt.show()
plt.close()

#If cos(theta13) =! 1


theta13 = np.radians(8.59)

plt.subplot(311)
plt.title('Effect of P$_{\u03BD_{\mu} -> \u03BD_{\mu}}^1$')

plt.bar(E_lisseND,Phi_lisseND,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label='$\phi_{\u03BD_{\mu}}^{ND}$') #Plot Input data
plt.xlabel('E GeV')
plt.ylabel('$\phi_{\u03BD_{\mu}}^{ND}$    ${\u03BD_{\mu}}$/$cm^{2}$/POT/1 GeV')
plt.legend()

Ex = np.linspace(0.5,Emax,1000)

plt.subplot(312)

plt.plot( Ex, osc_mod.proba_disapp_1 (theta13, theta23, Ex , 1285 , True ) ,label='P$_{\u03BD_{\mu} -> \u03BD_{\mu}}^1$')
plt.xlabel('E GeV')
plt.ylabel('P$_{\u03BD_{\mu} -> \u03BD_{\mu}}^1$')
plt.legend()

plt.subplot(313)


Phi_nu_mu_FD1 = Phi_lisseND * osc_mod.proba_disapp_1 ( theta13,theta23, E_lisseND , 1285 , True ) 
plt.bar(E_lisseND, Phi_nu_mu_FD1 ,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label='$\phi_{\u03BD_{\mu}}^{ND}$ * P$_{\u03BD_{\mu} -> \u03BD_{\mu}}^1$')

plt.xlabel('E GeV')
plt.ylabel('$\phi_{\u03BD_{\mu}}^{ND}$ x P$_{\u03BD_{\mu} -> \u03BD_{\mu}}^1$   ${\u03BD_{\mu}}$/$cm^{2}$/POT/1 GeV')


plt.legend()

plt.show()
plt.close()

#  III )                   NEAR TO FAR FLUX EXTRAPOLATION: DATA EXTRACTION 


F_nd_fd = np.loadtxt('data/F_nd_fd.txt')
ene_article, F_article = np.array([i[0] for i in F_nd_fd]), np.array([i[1]*1e-06 for i in F_nd_fd])

plt.plot(ene_article,F_article,'x',color="red",label='$F_{FD/ND}$')

F = interp1d(ene_article, F_article, kind='cubic')
F_lisse = F(E_lisseND)

#ene_lisse_F = osc_mod.model_flux(ene_article, F_article,1e06,0.3) 

#F_lisse=np.histogram(ene_lisse_F,Nbin)[0]
#F_lisse=np.array(F_lisse)/max(np.array(F_lisse))*2.9e-07

plt.bar(E_lisseND,F_lisse,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label='$F_{FD/ND}$')
plt.xlabel('E GeV')
plt.ylabel('$\phi_{\u03BD_{\mu}}^{FD}$ / $\phi_{\u03BD_{\mu}}^{ND}$')
plt.xlim(0.25,Emax)
#print(E_lisseND)

plt.title('Near to far flux extrapolation factor data smoothing')
plt.legend()
plt.show()
plt.close()



# IV )                   FLUX OF nu_mu AT THE FD (WITH THE NEAR TO FAR FLUX EXTRAPOLATION ; WITHOUT THE ENERGY RECONSTRUCTION FACTOR)


#Plot of the effect of the near to far flux extrapolation

plt.subplot(221) # Phi_nu_mu_ND
plt.title('Effect of the near to far flux extrapolation factor')

Npot=1.1e21 #Proton on target

plt.bar(E_lisseND, Phi_lisseND ,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin , label = '$\phi_{\u03BD_{\mu}}^{ND}$')
plt.xlabel('E GeV')
plt.ylabel('$\phi_{\u03BD_{\mu}}^{ND}$   ${\u03BD_{\mu}}$/$cm^{2}$/POT/1GeV')
plt.legend()



plt.subplot(223) # Ph_nu_mu_ND x F
plt.bar(E_lisseND, Phi_lisseND * F_lisse,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin , label = '$\phi_{\u03BD_{\mu}}^{ND}$ x F')
plt.xlabel('E GeV')
plt.ylabel('$\phi_{\u03BD_{\mu}}^{ND}$ x F')
plt.legend()


plt.subplot(222) # Ph_nu_mu_ND x Papp

plt.bar(E_lisseND,Phi_nu_mu_FD0,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin ,label='$\phi_{\u03BD_{\mu}}^{ND}$ x P$_{\u03BD_{\mu} -> \u03BD_{\mu}}$') #One has now nu_e/cm^2 per 1GeV
plt.xlabel('E GeV')
plt.ylabel('$\phi_{\u03BD_{\mu}}^{ND}$xP$_{\u03BD_{\mu}->\u03BD_e}$xF   ${\u03BD_e}$/$cm^{2}$/POT/1GeV')
plt.legend()

plt.subplot(224) # Ph_nu_mu_ND x Papp x F

Phi_nu_mu_det_FD0 = Phi_nu_mu_FD0 * F_lisse
plt.bar(E_lisseND,Phi_nu_mu_det_FD0,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin ,label='$\phi_{\u03BD_{\mu}}^{ND}$ x P$_{\u03BD_{\mu} -> \u03BD_{\mu}}^0$ x F') #One has now nu_e/cm^2 per 1GeV
plt.xlabel('E GeV')
plt.ylabel('$\phi_{\u03BD_{\mu}}^{ND}$xP$_{\u03BD_{\mu}->\u03BD_{\mu}}^0$xF')
plt.legend()


plt.legend()
plt.show()
plt.close()


#If cos(theta13) != 1
Phi_nu_mu_det_FD1 = Phi_nu_mu_FD1 * F_lisse
plt.bar(E_lisseND,Phi_nu_mu_det_FD1,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin ,label='$\phi_{\u03BD_{\mu}}^{ND}$ x P$_{\u03BD_{\mu} -> \u03BD_{\mu}}^1$ x F')
plt.xlabel('E GeV')
plt.ylabel('$\phi_{\u03BD_{\mu}}^{ND}$xP$_{\u03BD_{\mu}->\u03BD_{\mu}}^1$xF')
plt.legend()
plt.show()
plt.close()


# V )			NUMBER OF nu_mu AT THE FD (WITH THE NEAR TO FAR FLUX EXTRAPOLATION ; WITHOUT THE ENERGY RECONSTRUCTION FACTOR)


#Flux x cross section 

sigma = 1e-38 #cm^2 per nucleon
"""
sigma_o_E = np.loadtxt('data/sigma-E-2103-04797_7GeV.txt')
E, sigma_o_E_nue_data = np.array([i[0] for i in sigma_o_E]), np.array([i[1] for i in sigma_o_E])*1e-38
sigma_o_E_nue_interp = interp1d(E, sigma_o_E_nue_data, kind='cubic')
sigma_nue = np.array(sigma_o_E_nue_interp(E_lisseND)) * E_lisseND 
"""

sigma_o_E = np.loadtxt('data/xsection_cc.txt')
log10E, sigma_o_E_nue_data = np.array([i[0] for i in sigma_o_E]), np.array([i[2] for i in sigma_o_E])*1e-38
E = pow( 10 , log10E)
sigma_o_E_nue_interp = interp1d(E, sigma_o_E_nue_data, kind='cubic')
sigma_nue = np.array(sigma_o_E_nue_interp(E_lisseND)) * E_lisseND


M_Ar = 40 #Molar mass of Argon in g/mol
Na = 6e23 #Avogadro cst in mol^-1
m_Ar_FD = 40e9 #masse of argon in the FD
N_nuc_per_Ar = 40 #Number of nucleon per atom of Ar

sigma_FD = m_Ar_FD * Na/M_Ar * N_nuc_per_Ar * sigma_nue #cross section of the reaction ( nu_e + n -> e- + p ) in the FD (in cm^2), same for nu_mu than for nu_e

plt.subplot(311) #Smoothing of the cross section
plt.title('Comparison of $\phi_{\u03BD_{\mu}}^{ND}$ * P$_{\u03BD_{\mu} -> \u03BD_{\mu}}$ * F * $\sigma_{\u03BD_{\mu}}^{Ar}$   vs   N$_{\u03BD_{\mu}}^{FD}$')
plt.plot(E, m_Ar_FD * Na/M_Ar * N_nuc_per_Ar * sigma_o_E_nue_data * E, 'x', label='$\sigma_{\u03BD_e}^{Ar}$  data')
plt.bar(E_lisseND, sigma_FD,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label='$\sigma_{\u03BD_{\mu}}^{Ar}$  smoothed')
plt.legend()
plt.xlabel('E GeV')
plt.ylabel('$\sigma_{\u03BD_{\mu}}^{Ar}$   cm^2')


plt.subplot(312)

N_nu_mu_det0 = Phi_nu_mu_det_FD0* sigma_FD
plt.bar(E_lisseND,N_nu_mu_det0 * Npot,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label='$\phi_{\u03BD_{\mu}}^{ND}$ * P$_{\u03BD_{\mu} -> \u03BD_{\mu}^0}$ * F * $\sigma_{\u03BD_{\mu}}^{Ar}$') #One has now nu_mu per 1GeV

plt.xlabel('E GeV')
plt.ylabel('N$_{\u03BD_{\mu}}^{FD}$   ${\u03BD_{\mu}}$/1GeV (x1.1*10$^{21}$POT)')

print( N_nu_mu_det0 * Npot )


#DATA OF Nnu_mu at the FD

plt.plot(Emu_FD_article, specmu_FD_article,'x',color='red', label='N$_{\u03BD_{\mu}}^{FD-data}$')

plt.legend()


#If cos(theta13) =! 1
plt.subplot(313)

N_nu_mu_det1 = Phi_nu_mu_det_FD1* sigma_FD
plt.bar(E_lisseND,N_nu_mu_det1 * Npot,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label='$\phi_{\u03BD_{\mu}}^{ND}$ * P$_{\u03BD_{\mu} -> \u03BD_{\mu}^1}$ * F * $\sigma_{\u03BD_{\mu}}^{Ar}$') #One has now nu_mu per 1GeV

plt.xlabel('E GeV')
plt.ylabel('N$_{\u03BD_{\mu}}^{FD}$   ${\u03BD_{\mu}}$/1GeV (x1.1*10$^{21}$POT)')


#DATA OF Nnu_mu at the FD

plt.plot(Emu_FD_article, specmu_FD_article,'x',color='red', label='N$_{\u03BD_{\mu}}^{FD-data}$')

plt.legend()

plt.show()
plt.close()



# VI)			RECONSTRUCTED ENERGY




#print (E_lisseND)
#print(	N_nu_mu_det0 * Npot)

# To implement the reconstructed-energy factor, the statistic has to be important. The numer of nu_e was untill now in number of nu_e per POT, it is now in numbre of nu_e (*1.1*10e21) and is even multiply per Nmultiply (temporarily, after having implemented the reconstructed-energy factor, it will be divides by Nmultiply)

"""	
plt.subplot(211)
"""
plt.bar(E_lisseND , np.array(osc_mod.Rec2_numu(Reco, E_lisseND, N_nu_mu_det0 * Npot * Nmultiply , 0)) / Nmultiply,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label='$\phi_{\u03BD_{\mu}}^{ND}$ * P$_{\u03BD_{\mu} -> \u03BD_{\mu}}$ * F * $\sigma_{\u03BD_{\mu}}^{Ar}$ *T with Rec2')	

plt.plot(Emu_FD_article, np.array(specmu_FD_article),'x',color='red', label='N$_{\u03BD_{\mu}}^{FD-data}$')
plt.legend()
plt.xlabel('E GeV')
plt.ylabel('N$_{\u03BD_{\mu}}^{FD}$   ${\u03BD_{\mu}}$/1GeV (x1.1*10$^{21}$POT) with Rec2')
plt.show()
plt.close()




# VII ) The output data are smoothed to calculate the integral of output as the same way as Phi_nu_u_ND * Papp * F * sigma_FD *T 



plt.subplot(211)
plt.plot( ene_FD_art_mu_noise, spec_FD_art_mu_noise , 'x' , color='blue' , label='noise')
plt.plot(Emu_FD_article, specmu_FD_article,'x',color='red', label='N$_{\u03BD_{\mu}}^{FD-data}$ with noise') #Output data
plt.title('Output data smoothing')
plt.xlabel('E GeV')
plt.ylabel('N$_{\u03BD_{\mu}}^{FD}$   ${\u03BD_{\mu}}$/1GeV (x1.1*10$^{21}$POT)')
plt.legend()

plt.subplot(212)
#One has a list of smoothed energy data list with 1e06 output data (It is not those data that will be plot in reality)
list_FD = interp1d(Emu_FD_article, specmu_FD_article, kind='cubic') #with noise
list_FD_noise = interp1d(ene_FD_art_mu_noise, spec_FD_art_mu_noise, kind='cubic') #The noise

N_FD_hist= list_FD(E_lisseND)  #with noise
N_FD_hist_noise = list_FD_noise(E_lisseND) #The noise

#(Output-with-noise data) - (Output-noise data) = Output-without-noise data
N_FD_hist = np.array(N_FD_hist) - np.array(N_FD_hist_noise)


plt.bar(E_lisseND, N_FD_hist , color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin , label = 'N$_{\u03BD_{\mu}}^{FD-data}$ smoothed without noise') #Output smoothed data
plt.legend()
plt.show()
plt.close()



# VIII ) Energy distribution at the FD: calculated from the input VS output  FINAL PLOT



DN_nu_mu_sat = np.sqrt(N_FD_hist) #Statistic errors in srqt (N) (Poisson's law)
DN_nu_mu_sys = DN_nu_mu_sat #One takes the systematic errors equals to the statisticq one
DN_nu_mu = DN_nu_mu_sat + DN_nu_mu_sys

theta23_1 = 0.69
theta23_2= 0.78
theta23_3= 0.87

delta_m23_1= 2e-03
delta_m23_2= 2.45e-03
delta_m23_3= 3e-03

plt.subplot(331) # theta23_1 , delta_m23_1

N1=np.array(osc_mod.Rec2_numu(Reco, E_lisseND, Phi_lisseND * osc_mod.proba_disapp_0_t ( delta_m23_1 ,theta23_1, E_lisseND , 1285 , True ) * F_lisse * sigma_FD* Npot * Nmultiply , bias)) / Nmultiply
chi2_1 = osc_mod.chi2( N1 , N_FD_hist , E_lisseND)
chi2_str_1=str(round(chi2_1,2))
plt.bar(E_lisseND , N1,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin ,label=(f"$\chi^2$={chi2_str_1}"))
plt.plot(E_lisseND, N_FD_hist,'x',color='red', label='N$_{\u03BD_{\mu}}^{FD-data}$')
plt.errorbar(E_lisseND, N_FD_hist, DN_nu_mu, fmt='none',ecolor ='blue')
plt.ylabel('$\Delta$$m^2_{23} = 2e-03$ ')
plt.legend()

plt.subplot(332) # theta23_2 , delta_m23_1
plt.title('N$_{\u03BD_{\mu}}^{FD}$   ${\u03BD_{\mu}}$/1GeV (x1.1*10$^{21}$POT)')
N2=np.array(osc_mod.Rec2_numu(Reco, E_lisseND, Phi_lisseND * osc_mod.proba_disapp_0_t ( delta_m23_1 ,theta23_2, E_lisseND , 1285 , True ) * F_lisse * sigma_FD* Npot * Nmultiply , bias)) / Nmultiply
chi2_2 = osc_mod.chi2( N2 , N_FD_hist , E_lisseND)
chi2_str_2=str(round(chi2_2,2))
plt.bar(E_lisseND , N2 ,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label=(f"$\chi^2$={chi2_str_2}"))
plt.plot(E_lisseND, N_FD_hist,'x',color='red', label='N$_{\u03BD_{\mu}}^{FD-data}$')
plt.errorbar(E_lisseND, N_FD_hist, DN_nu_mu, fmt='none',ecolor ='blue')
plt.legend()

plt.subplot(333) # theta23_3 , delta_m23_1

N3=np.array(osc_mod.Rec2_numu(Reco, E_lisseND, Phi_lisseND * osc_mod.proba_disapp_0_t ( delta_m23_1 ,theta23_3, E_lisseND , 1285 , True ) * F_lisse * sigma_FD* Npot * Nmultiply , bias)) / Nmultiply
chi2_3 = osc_mod.chi2( N3 , N_FD_hist , E_lisseND)
chi2_str_3=str(round(chi2_3,2))
plt.bar(E_lisseND , N3 ,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label=(f"$\chi^2$={chi2_str_3}"))
plt.plot(E_lisseND, N_FD_hist,'x',color='red', label='N$_{\u03BD_{\mu}}^{FD-data}$')
plt.errorbar(E_lisseND, N_FD_hist, DN_nu_mu, fmt='none',ecolor ='blue')
plt.legend()

plt.subplot(334) # theta23_1 , delta_m23_2

N4=np.array(osc_mod.Rec2_numu(Reco, E_lisseND, Phi_lisseND * osc_mod.proba_disapp_0_t ( delta_m23_2 ,theta23_1, E_lisseND , 1285 , True ) * F_lisse * sigma_FD* Npot * Nmultiply , bias)) / Nmultiply
chi2_4 = osc_mod.chi2( N4 , N_FD_hist , E_lisseND)
chi2_str_4=str(round(chi2_4,2))
plt.bar(E_lisseND , N4 ,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label=(f"$\chi^2$={chi2_str_4}"))
plt.plot(E_lisseND, N_FD_hist,'x',color='red', label='N$_{\u03BD_{\mu}}^{FD-data}$')
plt.errorbar(E_lisseND, N_FD_hist, DN_nu_mu, fmt='none',ecolor ='blue')
plt.ylabel('$\Delta$$m^2_{23} = 2.45e-03$ (min/true value) ')
plt.legend()

plt.subplot(335) # theta23_2 , delta_m23_2

N5=np.array(osc_mod.Rec2_numu(Reco, E_lisseND, Phi_lisseND * osc_mod.proba_disapp_0_t ( delta_m23_2 ,theta23_2, E_lisseND , 1285 , True ) * F_lisse * sigma_FD* Npot * Nmultiply , bias)) / Nmultiply
chi2_5 = osc_mod.chi2( N5 , N_FD_hist , E_lisseND)
chi2_str_5=str(round(chi2_5,2))
plt.bar(E_lisseND , N5 ,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label=(f"$\chi^2$={chi2_str_5}"))
plt.plot(E_lisseND, N_FD_hist,'x',color='red', label='N$_{\u03BD_{\mu}}^{FD-data}$')
plt.errorbar(E_lisseND, N_FD_hist, DN_nu_mu, fmt='none',ecolor ='blue')
plt.legend()

plt.subplot(336) # theta23_3 , delta_m23_2

N6=np.array(osc_mod.Rec2_numu(Reco, E_lisseND, Phi_lisseND * osc_mod.proba_disapp_0_t ( delta_m23_2 ,theta23_3, E_lisseND , 1285 , True ) * F_lisse * sigma_FD* Npot * Nmultiply , bias)) / Nmultiply
chi2_6 = osc_mod.chi2( N6 , N_FD_hist , E_lisseND)
chi2_str_6=str(round(chi2_6,2))
plt.bar(E_lisseND , N6 ,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label=(f"$\chi^2$={chi2_str_6}"))
plt.plot(E_lisseND, N_FD_hist,'x',color='red', label='N$_{\u03BD_{\mu}}^{FD-data}$')
plt.errorbar(E_lisseND, N_FD_hist, DN_nu_mu, fmt='none',ecolor ='blue')
plt.legend()

plt.subplot(337) # theta23_1 , delta_m23_3

N7=np.array(osc_mod.Rec2_numu(Reco, E_lisseND, Phi_lisseND * osc_mod.proba_disapp_0_t ( delta_m23_3 ,theta23_1, E_lisseND , 1285 , True ) * F_lisse * sigma_FD* Npot * Nmultiply , bias)) / Nmultiply
chi2_7 = osc_mod.chi2( N7 , N_FD_hist , E_lisseND)
chi2_str_7=str(round(chi2_7,2))
plt.bar(E_lisseND , N7 ,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label=(f"$\chi^2$={chi2_str_7}"))
plt.plot(E_lisseND, N_FD_hist,'x',color='red', label='N$_{\u03BD_{\mu}}^{FD-data}$')
plt.errorbar(E_lisseND, N_FD_hist, DN_nu_mu, fmt='none',ecolor ='blue')
plt.legend()
plt.xlabel('$\Theta_{23} = 0.69 rad$  (min value)')
plt.ylabel('$\Delta$$m^2_{23} = 3e-03$ ')

plt.subplot(338) # theta23_2 , delta_m23_3

N8=np.array(osc_mod.Rec2_numu(Reco, E_lisseND, Phi_lisseND * osc_mod.proba_disapp_0_t ( delta_m23_3 ,theta23_2, E_lisseND , 1285 , True ) * F_lisse * sigma_FD* Npot * Nmultiply ,bias)) / Nmultiply
chi2_8 = osc_mod.chi2( N8 , N_FD_hist , E_lisseND)
chi2_str_8=str(round(chi2_8,2))
plt.bar(E_lisseND , N8 ,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label=(f"$\chi^2$={chi2_str_8}"))
plt.plot(E_lisseND, N_FD_hist,'x',color='red', label='N$_{\u03BD_{\mu}}^{FD-data}$')
plt.errorbar(E_lisseND, N_FD_hist, DN_nu_mu, fmt='none',ecolor ='blue')
plt.legend()
plt.xlabel('$\Theta_{23} = 78 rad$')

plt.subplot(339) # theta23_3 , delta_m23_3

N9=np.array(osc_mod.Rec2_numu(Reco, E_lisseND, Phi_lisseND * osc_mod.proba_disapp_0_t ( delta_m23_3 ,theta23_3, E_lisseND , 1285 , True ) * F_lisse * sigma_FD* Npot * Nmultiply , bias)) / Nmultiply
chi2_9 = osc_mod.chi2( N9 , N_FD_hist , E_lisseND)
chi2_str_9=str(round(chi2_9,2))
plt.bar(E_lisseND , N9 ,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label=(f"$\chi^2$={chi2_str_9}"))
plt.plot(E_lisseND, N_FD_hist,'x',color='red', label='N$_{\u03BD_{\mu}}^{FD-data}$')
plt.errorbar(E_lisseND, N_FD_hist, DN_nu_mu, fmt='none',ecolor ='blue')
plt.legend()
plt.xlabel('$\Theta_{23} = 0.87 rad$  (true value)')

plt.show()
plt.close()


#Plot of the true values spectrum
plt.bar(E_lisseND , N6 ,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label=(f"$\chi^2$={chi2_str_6}"))
plt.plot(E_lisseND, N_FD_hist,'x',color='red', label='N$_{\u03BD_{\mu}}^{FD-data}$')
plt.errorbar(E_lisseND, N_FD_hist, DN_nu_mu, fmt='none',ecolor ='blue')
plt.legend()
plt.title('N$_{\u03BD_{\mu}}^{FD}$  for $\Theta_{23} = 0.87 rad$ and $\Delta$$m^2_{23} = 2.45e-03 eV^{2}$')
plt.xlabel('E (GeV)')
plt.ylabel('N$_{\u03BD_{\mu}}$')
plt.show()
plt.close()


print(np.array(osc_mod.Rec2_numu( Reco,E_lisseND, N_nu_mu_det0 * Npot * Nmultiply , bias)) / Nmultiply )




#  IX)   KHI2 STUDY 



# IX) a) khi 2 function of sin^2(theta23) (The renormalisation factor is constant)


L_theta23 = np.linspace(-3.14,3.14,360)
L_sin23squar = np.sin(L_theta23)**2
L_N_nue = []
khi2=[]

for i in range (0,len(L_theta23)):

	L_N_nue = L_N_nue + [ np.array(osc_mod.Rec2_numu(Reco, E_lisseND, Phi_lisseND * osc_mod.proba_disapp_0_t ( osc_par.delta_m32,L_theta23[i], E_lisseND , 1285 , True )  * F_lisse * sigma_FD * Npot *Nmultiply, bias))/Nmultiply ]

	khi2 = khi2 + [0]
	for j in range(0,len(E_lisseND)):
		if L_N_nue[i][j] != 0:
			khi2[i] = khi2[i] + ( L_N_nue[i][j] - N_FD_hist[j] )**2/(4*L_N_nue[i][j])
		else :
			khi2[i] = khi2[i]


plt.plot( L_theta23 , khi2 , label='$\chi^2$ with P$_{\u03BD_{\mu} -> \u03BD_{\mu}}^0$')
plt.xlabel('$\Theta_{23}$')
plt.ylabel('$\chi^2$')
plt.yscale('log')
plt.legend()
plt.show()
plt.close()

plt.plot( L_sin23squar , khi2 , label='$\chi^2$ with P$_{\u03BD_{\mu} -> \u03BD_{\mu}}^0$')
plt.xlabel('$sin(\Theta_{23})^2$')
plt.ylabel('$\chi^2$')
plt.yscale('log')
plt.legend()
plt.show()


# IX) b) chi2 function of Dm32

L_delta_m32 = np.linspace (1e-03,4e-03,50) 
L_N_nue = []
khi2=[]

for i in range (0,len(L_delta_m32)):

	L_N_nue = L_N_nue + [ np.array(osc_mod.Rec2_numu(Reco, E_lisseND, Phi_lisseND * osc_mod.proba_disapp_0_t ( L_delta_m32[i], osc_par.theta23 , E_lisseND , 1285 , True )  * F_lisse * sigma_FD * Npot *Nmultiply , bias))/Nmultiply ]

	khi2 = khi2 + [0]
	for j in range(0,len(E_lisseND)):
		if L_N_nue[i][j] != 0:
			khi2[i] = khi2[i] + ( L_N_nue[i][j] - N_FD_hist[j] )**2/(4*L_N_nue[i][j])
		else :
			khi2[i] = khi2[i]


plt.plot( L_delta_m32*1000 , khi2 , label='$\chi^2$ with P$_{\u03BD_{\mu} -> \u03BD_{\mu}}^0$')
plt.xlabel('$\Delta$$m^2_{23}$ / $10^{-3} ev^2$')
plt.ylabel('$\chi^2$')
plt.yscale('log')
plt.legend()
plt.show()
plt.close()

#  IX) c) khi2 function of sin2(theta23) and Dm23^2
"""
L_delta_m32 = np.linspace (2.32e-03,2.52e-03,300) 
L_theta23 = np.linspace(0,np.pi/2,300)
L_sin23 = np.sin(L_theta23)**2
LL2_N_nue = []
L2_khi2 = []


for i in range (0,len(L_theta23)):
	

	LL2_N_nue = LL2_N_nue +[[]]
	L2_khi2 = L2_khi2 + [[]]

	for j in range (0,len(L_delta_m32)) :
		LL2_N_nue[i] = LL2_N_nue[i] + [np.array(osc_mod.Rec2_numu(Reco, E_lisseND, Phi_lisseND * osc_mod.proba_disapp_0_t(L_delta_m32[j], L_theta23[i], E_lisseND, 1285 , True) * F_lisse * sigma_FD * Npot *Nmultiply , bias))/Nmultiply]
				
		L2_khi2[i] = L2_khi2[i] + [0]
		for k in range (0,len(E_lisseND)):
			if LL2_N_nue[i][j][k] != 0:
				L2_khi2[i][j] = L2_khi2[i][j] + ( LL2_N_nue[i][j][k] - N_FD_hist[k] )**2/(4*LL2_N_nue[i][j][k])
			else :
				L2_khi2[i][j] = L2_khi2[i][j]

# Some statistics values 
#For 2 parameters fitted
Dchi2_1s = 2.3 
Dchi2_3s = 11.8
Dchi2_5s = 29.8 
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
fig.colorbar(pc,ax=ax,label='$\chi^2$')
plt.xlabel('$\Delta$$m^2_{23}$ ')
plt.ylabel('$sin(\Theta_{23})^2$ ')
plt.scatter(2.45*pow(10,-3)*1000 , np.sin(np.radians(49.8))**2 , s=10, color='red' )
plt.title('$\epsilon = 50 kT.MW.yrs$')
plt.plot(L_delta_m32[osc_mod.ind_min(L2_khi2)[1]] *1000, L_sin23[osc_mod.ind_min(L2_khi2)[0]] ,'x',color='red',label=(f"$\chi^2/Ndof$={chi2min_dof_str}"))
plt.legend(legend, legend_label, loc ='lower right')
plt.show()

"""

"""
# IX) d) Chi2 function of the Accuracy of the reconstruction of the energy


L_accuracy = np.linspace(10,24,70)
L_N_nue = []

khi2=[]

for i in range (0,len(L_accuracy)):

	L_N_nue = L_N_nue + [ np.array(osc_mod.Rec2_numu(L_accuracy[i], E_lisseND, Phi_lisseND * osc_mod.proba_disapp_0_t ( osc_par.delta_m32, osc_par.theta23 , E_lisseND , 1285 , True )  * F_lisse * sigma_FD * Npot *Nmultiply , 0))/Nmultiply ]

	
	khi2 = khi2 + [0]
	for j in range(0,len(E_lisseND)):
		if L_N_nue[i][j] != 0:
			khi2[i] = khi2[i] + ( L_N_nue[i][j] - N_FD_hist[j] )**2/(4*L_N_nue[i][j])
		else :
			khi2[i] = khi2[i]
			
plt.plot( L_accuracy , khi2 , label='$\chi^2$')
plt.xlabel('accuracy in %')
plt.ylabel('$\chi^2$')
#plt.yscale('log')
plt.legend()
plt.show()
plt.close()

# IX) e) Chi2 function of the bias of the reconstruction of the energy

Lbias = np.linspace(-15,10,50)
L_Ene_bias = []
L_N_nue = []

khi2=[]

for i in range (0,len(Lbias)):

	L_N_nue = L_N_nue + [ np.array(osc_mod.Rec2_numu(Reco, E_lisseND, Phi_lisseND * osc_mod.proba_disapp_0_t ( osc_par.delta_m32, osc_par.theta23 , E_lisseND , 1285 , True )  * F_lisse * sigma_FD * Npot *Nmultiply , Lbias[i]))/Nmultiply ]

	
	khi2 = khi2 + [0]
	for j in range(0,len(E_lisseND)):
		if L_N_nue[i][j] != 0:
			khi2[i] = khi2[i] + ( L_N_nue[i][j] - N_FD_hist[j] )**2/(4*L_N_nue[i][j])
		else :
			khi2[i] = khi2[i]
			
plt.plot( Lbias , khi2 , label='$\chi^2$')
plt.xlabel('bias in %')
plt.ylabel('$\chi^2$')
#plt.yscale('log')
plt.legend()
plt.show()
plt.close()
"""
# IX) f) Chi2 function of the bias and the accuracy of the reconstruction of the energy

#True values for theta23 and deltam23

L_bias = np.linspace(-15,10,100)
L_Reco= np.linspace(8,25,100)
LL2_N_nue = []
L2_khi2 = []


for i in range (0,len(L_bias)):
	

	LL2_N_nue = LL2_N_nue +[[]]
	L2_khi2 = L2_khi2 + [[]]

	for j in range (0,len(L_Reco)) :
		LL2_N_nue[i] = LL2_N_nue[i] + [np.array(osc_mod.Rec2_numu(L_Reco[j], E_lisseND, Phi_lisseND * osc_mod.proba_disapp_0_t( osc_par.delta_m32, osc_par.theta23, E_lisseND, 1285 , True) * F_lisse * sigma_FD * Npot *Nmultiply , L_bias[i]))/Nmultiply]
				
		L2_khi2[i] = L2_khi2[i] + [0]
		for k in range (0,len(E_lisseND)):
			if LL2_N_nue[i][j][k] != 0:
				L2_khi2[i][j] = L2_khi2[i][j] + ( LL2_N_nue[i][j][k] - N_FD_hist[k] )**2/(4*LL2_N_nue[i][j][k])
			else :
				L2_khi2[i][j] = L2_khi2[i][j]




# 2D histogram with pcolormesh

import matplotlib.colors as colors

Reco , BIAS = np.meshgrid ( L_Reco , L_bias )
fig, ax = plt.subplots()
pc = ax.pcolormesh(Reco , BIAS , L2_khi2,norm=colors.LogNorm(vmin=np.array(L2_khi2).min(), vmax=np.array(L2_khi2).max()))
fig.colorbar(pc,ax=ax,label='$\chi^2$')
plt.ylabel(' Bias in %')
plt.xlabel(' Reco in %')

chi2min_str = str(round( L2_khi2[osc_mod.ind_min(L2_khi2)[0]][osc_mod.ind_min(L2_khi2)[1]] , 2))
plt.plot(L_Reco[osc_mod.ind_min(L2_khi2)[1]] , L_bias[osc_mod.ind_min(L2_khi2)[0]] ,'x',color='red',label=(f"$\chi^2$={chi2min_str}"))
plt.legend()
plt.show()
plt.close()


N6=np.array(osc_mod.Rec2_numu(17.1, E_lisseND, Phi_lisseND * osc_mod.proba_disapp_0_t ( 0.00241 ,-0.6969, E_lisseND , 1285 , True ) * F_lisse * sigma_FD* Npot * Nmultiply , -4.9)) / Nmultiply
chi2_6 = osc_mod.chi2( N6 , N_FD_hist , E_lisseND)
chi2_str_6=str(round(chi2_6,2))
plt.bar(E_lisseND , N6 ,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label=(f"$\chi^2$={chi2_str_6}"))
plt.plot(E_lisseND, N_FD_hist,'x',color='red', label='N$_{\u03BD_{\mu}}^{FD-data}$')
plt.title('Bias=-3.8 % ; Reco = 16.3 %')
plt.errorbar(E_lisseND, N_FD_hist, DN_nu_mu, fmt='none',ecolor ='blue')
plt.legend()
plt.xlabel('$\Theta_{23} = 0.6969 rad$ ')
plt.ylabel('$\Delta$$m^2_{23} = 2.41e-03$ ')
plt.show()
plt.close()
