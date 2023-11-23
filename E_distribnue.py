"""
Created on 31th of March 2023
Matéo Bédès

Last edited 31th of March 2023
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
Nmultiply=30
Acc=20.77
bias = 2.31
Emax=7  #Has to be the same in lib_oscil.py


#   I )                      INPUT DATA SMOOTHING

#The input and output data are extracted and they are smoothed. The output data are also extracted


far_detector_spectrum = np.loadtxt('data/output7GeV.txt') #One extracts the output (WITH THE NOISE)	
ene_FD_article, spec_FD_article = np.array([i[0] for i in far_detector_spectrum]), np.array([i[1] for i in far_detector_spectrum])

fd_spectrum_noise = np.loadtxt('data/output_noise7GeV.txt') #Noise on the energy distribution of the Nve
ene_FD_art_noise, spec_FD_art_noise = np.array([i[0] for i in fd_spectrum_noise]), np.array([i[1] for i in fd_spectrum_noise])




#One extracts the data of the flux of nu_mu in the nu_mode

Phi_nu = np.loadtxt('data/Phi_mu_nu_mode_mu7GeV.txt')
Phi_nu_mu, E = np.array([i[0] for i in Phi_nu]), np.array([i[1] for i in Phi_nu])

plt.plot(E, Phi_nu_mu,'x', color = 'r',label='$\phi_{\u03BD_{\mu}}^{ND}$') #One plot the data of the nu_mu flux
plt.title('Input data smoothing')

plt.yscale('log')
plt.xlabel('E GeV')
plt.ylabel('$\phi_{\u03BD_{\mu}}^{ND}$    ${\u03BD_{\mu}}$/$cm^{2}$/POT/1 GeV')
plt.legend()


#list_ND_mc = osc_mod.model_flux(E, Phi_nu_mu,1e06, 0.3)

 #One has a list of smoothed energy data list with 1e06 input data (It is not those data that will be plot in reality)

"""
N=[]
for i in range(len(spec_ND_article)):
	N=N+[spec_ND_article[i]/2.7e06*12000]
#plt.plot(ene_ND_article, N,'x',label='données input brut')
"""

#plt.hist(list_ND_mc,200,label='données input lissées') 

#Input data are smoothed by interp1D

list_ND_mc = interp1d(E, Phi_nu_mu, kind='cubic')
#One stocks the smoothed input data
E_lisseND = np.linspace(0.5, Emax,Nbin)
Phi_lisseND = list_ND_mc(E_lisseND)
Phi_lisseND = np.array(Phi_lisseND)

#E_lisseND=np.histogram(list_ND_mc,Nbin)[1]
#E_lisseND=E_lisseND[:-1] #On retire la dernière donnée
#Phi_lisseND=np.histogram(list_ND_mc,Nbin)[0] #ene_lisse et N_lisse sont les données lissée!
#Phi_lisseND=np.array(Phi_lisseND)/max(np.array(Phi_lisseND))*max(Phi_nu_mu)  #Histogramme renormalisé sur les données de l'input



plt.bar(E_lisseND,Phi_lisseND,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label='$\phi_{\u03BD_{\mu}}^{ND}$')
plt.yscale('log')
plt.legend()
plt.show()
plt.close()


#  II )                    OSCILLATION PROBABILITY x SMOOTHED INPUT DATA WITHOUT THE NEAR TO FAR FLUX EXTRAPOLATION

theta12 = np.radians(33.82)
theta13 = np.radians(8.59)
theta23 = np.radians(49.8)

plt.subplot(311)
plt.title('Effect of P$_{\u03BD_{\mu} -> \u03BD_e}$')

plt.bar(E_lisseND,Phi_lisseND,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label='$\phi_{\u03BD_{\mu}}^{ND}$') #Plot Input data
plt.xlabel('E GeV')
plt.ylabel('$\phi_{\u03BD_{\mu}}^{ND}$    ${\u03BD_{\mu}}$/$cm^{2}$/POT/1 GeV')
plt.legend()

Ex = np.linspace(0.5,Emax,1000)
plt.subplot(312)
plt.plot( Ex, osc_mod.probability_oscillation( theta12 , theta13 , theta23, Ex, 1285 ,0, True, True),label='P$_{\u03BD_{\mu} -> \u03BD_e}$')
plt.xlabel('E GeV')
plt.ylabel('P$_{\u03BD_{\mu} -> \u03BD_e}$')
plt.legend()



plt.subplot(313)

Phi_mue_FD = Phi_lisseND * osc_mod.probability_oscillation( theta12 , theta13 , theta23,E_lisseND, 1285 ,0, True, True)
plt.bar(E_lisseND,Phi_mue_FD,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label='$\phi_{\u03BD_{\mu}}^{ND}$ x P$_{\u03BD_{\mu} -> \u03BD_e}$') #Plot du spectre au FD th

plt.xlabel('E GeV')
plt.ylabel('$\phi_{\u03BD_{\mu}}^{ND}$ x P$_{\u03BD_{\mu} -> \u03BD_e}$   ${\u03BD_e}$/$cm^{2}$/POT/1 GeV')


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



# IV )                   FLUX OF mu_e AT THE FD (WITH THE NEAR TO FAR FLUX EXTRAPOLATION ; WITHOUT THE ENERGY RECONSTRUCTION FACTOR)

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

plt.bar(E_lisseND,Phi_mue_FD,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin ,label='$\phi_{\u03BD_{\mu}}^{ND}$ x P$_{\u03BD_{\mu} -> \u03BD_e}$') #One has now nu_e/cm^2 per 1GeV
plt.xlabel('E GeV')
plt.ylabel('$\phi_{\u03BD_{\mu}}^{ND}$xP$_{\u03BD_{\mu}->\u03BD_e}$xF   ${\u03BD_e}$/$cm^{2}$/POT/1GeV')
plt.legend()

plt.subplot(224) # Ph_nu_mu_ND x Papp x F

Phi_mue_det_FD = Phi_mue_FD * F_lisse
plt.bar(E_lisseND,Phi_mue_det_FD,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin ,label='$\phi_{\u03BD_{\mu}}^{ND}$ x P$_{\u03BD_{\mu} -> \u03BD_e}$ x F') #One has now nu_e/cm^2 per 1GeV
plt.xlabel('E GeV')
plt.ylabel('$\phi_{\u03BD_{\mu}}^{ND}$xP$_{\u03BD_{\mu}->\u03BD_e}$xF')
plt.legend()


plt.legend()
plt.show()
plt.close()


# V )			NUMBER OF mu_e AT THE FD (WITH THE NEAR TO FAR FLUX EXTRAPOLATION ; WITHOUT THE ENERGY RECONSTRUCTION FACTOR)

#Flux x cross section 

sigma = 1e-38 #cm^2 per nucleon

"""
sigma_o_E = np.loadtxt('data/sigma-E-2103-04797_7GeV.txt')
E, sigma_o_E_nue_data = np.array([i[0] for i in sigma_o_E]), np.array([i[1] for i in sigma_o_E])*1e-38
sigma_o_E_nue_interp = interp1d(E, sigma_o_E_nue_data, kind='cubic')
sigma_nue = np.array(sigma_o_E_nue_interp(E_lisseND)) * E_lisseND 
"""

sigma_o_E = np.loadtxt('data/xsection_cc.txt')
log10E, sigma_o_E_nue_data = np.array([i[0] for i in sigma_o_E]), np.array([i[1] for i in sigma_o_E])*1e-38
E = pow( 10 , log10E)
sigma_o_E_nue_interp = interp1d(E, sigma_o_E_nue_data, kind='cubic')
sigma_nue = np.array(sigma_o_E_nue_interp(E_lisseND)) * E_lisseND


M_Ar = 40 #Molar mass of Argon in g/mol
Na = 6e23 #Avogadro cst in mol^-1
m_Ar_FD = 40e9 #masse of argon in the FD
N_nuc_per_Ar = 40 #Number of nucleon per atom of Ar

sigma_FD = m_Ar_FD * Na/M_Ar * N_nuc_per_Ar * sigma_nue #cross section of the reaction ( nu_e + n -> e- + p ) in the FD (in cm^2)

plt.subplot(211)
plt.title('Comparison of $\phi_{\u03BD_{\mu}}^{ND}$ * P$_{\u03BD_{\mu} -> \u03BD_e}$ * F * $\sigma_{\u03BD_e}^{Ar}$   vs   N$_{\u03BD_e}^{FD}$')
plt.plot(E, m_Ar_FD * Na/M_Ar * N_nuc_per_Ar * sigma_o_E_nue_data * E, 'x', label='$\sigma_{\u03BD_e}^{Ar}$  data')
plt.bar(E_lisseND, sigma_FD,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label='$\sigma_{\u03BD_e}^{Ar}$  smoothed')
plt.legend()
plt.xlabel('E GeV')
plt.ylabel('$\sigma_{\u03BD_e}^{Ar}$   cm^2')


plt.subplot(212)


N_mue_det = Phi_mue_det_FD * sigma_FD
plt.bar(E_lisseND,N_mue_det * Npot,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label='$\phi_{\u03BD_{\mu}}^{ND}$ * P$_{\u03BD_{\mu} -> \u03BD_e}$ * F * $\sigma_{\u03BD_e}^{Ar}$') #One has now nu_e per 1GeV

plt.xlabel('E GeV')
plt.ylabel('N$_{\u03BD_e}^{FD}$   ${\u03BD_e}$/1GeV (x1.1*10$^{21}$POT)')


plt.plot(ene_FD_article, np.array(spec_FD_article),'x',color='red', label='N$_{\u03BD_e}^{FD-data}$')
plt.legend()

plt.show()
plt.close()


#print(len(E_lisseND))



# VI)			RECONSTRUCTED ENERGY




# To implement the reconstructed-energy factor, the statistic has to be important. The numer of nu_e was untill now in number of nu_e per POT, it is now in numbre of nu_e (*1.1*10e21) and is even multiply per Nmultiply (temporarily, after having implemented the reconstructed-energy factor, it will be divides by Nmultiply)
	

plt.bar(E_lisseND , np.array(osc_mod.Rec2(Acc, E_lisseND, N_mue_det * Npot * Nmultiply , bias)) / Nmultiply,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label='$\phi_{\u03BD_{\mu}}^{ND}$ * P$_{\u03BD_{\mu} -> \u03BD_e}$ * F * $\sigma_{\u03BD_e}^{Ar}$ *T with Rec2')	

plt.plot(ene_FD_article, np.array(spec_FD_article)/max(np.array(spec_FD_article))*max(np.array(osc_mod.Rec2(Acc, E_lisseND, N_mue_det * Npot , bias))),'x',color='red', label='N$_{\u03BD_e}^{FD-data}$')
#plt.plot(ene_FD_article, np.array(spec_FD_article) ,'x',color='red', label='N$_{\u03BD_e}^{FD-data}$')
plt.ylabel('N$_{\u03BD_e}^{FD}$   ${\u03BD_e}$/1GeV (x1.1*10$^{21}$POT)')
plt.xlabel('E GeV')
plt.legend()


plt.title('Normalization of N$_{\u03BD_e}^{FD-th}$ by the max ratio')

##plt.plot(ene_FD_article, np.array(spec_FD_article)/max(np.array(spec_FD_article))*max(np.array(osc_mod.Rec2(Acc, E_lisseND, N_mue_det * Npot * Nmultiply )) / Nmultiply),'x',color='red', label='N$_{\u03BD_e}^{FD-data}$')


plt.show()
plt.close()



# VII ) The output data are smoothed to calculate the integral of output as the same way as Phi_nu_u_ND * Papp * F * sigma_FD *T 

plt.subplot(211)
plt.plot(ene_FD_article, np.array(spec_FD_article),'x',color='red', label='N$_{\u03BD_e}^{FD-data}$ with noise') #Output data
plt.plot(ene_FD_art_noise, spec_FD_art_noise, 'x', color='blue', label='Background') #Noise data
plt.title('Output data smoothing')

plt.ylabel('N$_{\u03BD_e}^{FD}$   ${\u03BD_e}$/1GeV (x1.1*10$^{21}$POT)')
plt.legend()

#Output-with-noise data smoothing
list_FD = interp1d(ene_FD_article, spec_FD_article, kind='cubic')
N_FD_hist= list_FD(E_lisseND)

#Output-noise data smoothing
list_FD_noise = interp1d(ene_FD_art_noise, spec_FD_art_noise, kind='cubic')
N_FD_hist_noise = list_FD_noise(E_lisseND)

#(Output-with-noise data) - (Output-noise data) = Output-without-noise data
N_FD_hist = np.array(N_FD_hist) - np.array(N_FD_hist_noise)

plt.subplot(212)
#print(osc_mod.cumulative(N_FD_hist))
#print(osc_mod.cumulative(N_FD_hist)[-1])
plt.bar(E_lisseND, N_FD_hist , color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin , label = 'N$_{\u03BD_e}^{FD-data}$ smoothed without background') #Output smoothed data
plt.xlabel('E GeV')
plt.legend()
plt.show()
plt.close()



# VIII ) Energy distribution at the FD: calculated from the input VS output  FINAL PLOT : For differents values of DeltaCP and theta23

DN_nu_e_sat = np.sqrt(N_FD_hist) #Statistic errors in srqt (N) (Poisson's law)
DN_nu_e_sys = DN_nu_e_sat #One takes the systematic errors equals to the statisticq one
DN_nu_e = DN_nu_e_sat + DN_nu_e_sys

#Generate a random nu_e energy distribution spectrum at the FD  --->  One pseudo experience 
Nexp1=[]
for i in range(0,len(DN_nu_e)):
	Nexp1.append(rand.uniform(N_FD_hist[i]-DN_nu_e[i]/2,N_FD_hist[i]+DN_nu_e[i]/2))
plt.bar(E_lisseND , Nexp1 ,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label='random N$_{\u03BD_e}^{FD-exp}$')
plt.errorbar(E_lisseND, Nexp1, DN_nu_e, fmt='none',ecolor ='blue')
plt.xlabel('E GeV')
plt.ylabel('N$_{\u03BD_e}^{FD}$   ${\u03BD_e}$/1GeV (x1.1*10$^{21}$POT)')
plt.legend()
plt.title('Generating one pseudo experience')
plt.show()
plt.close()

theta23_1=0.80
theta23_2=0.87
theta23_3=0.87
deltaCP_1=-0.45
deltaCP_2=0
deltaCP_3=-0.35



plt.subplot(321) # deltaCP_1 ;theta23_1 

N_1 = np.array(osc_mod.Rec2(Acc, E_lisseND, Phi_lisseND * osc_mod.probability_oscillation( theta12 , theta13 , theta23_1 ,E_lisseND, 1285 ,deltaCP_1*180/np.pi, True, True) * F_lisse * sigma_FD * Npot *Nmultiply ,bias))/Nmultiply
#chi2_1 = osc_mod.chi2( N_1 , N_FD_hist , E_lisseND)
chi2_1 = osc_mod.chi2( N_1 , Nexp1 , E_lisseND)
chi2_str_1=str(round(chi2_1,2))
plt.bar(E_lisseND , N_1 ,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label=(f"$\chi^2$={chi2_str_1}"))
#plt.plot(E_lisseND, N_FD_hist ,'x',color='red', label='N$_{\u03BD_e}^{FD-data}$')
#plt.errorbar(E_lisseND, N_FD_hist, DN_nu_e, fmt='none',ecolor ='blue')
plt.plot(E_lisseND, Nexp1 ,'x',color='red', label='N$_{\u03BD_e}^{FD-data}$')
plt.errorbar(E_lisseND,Nexp1 , DN_nu_e, fmt='none',ecolor ='blue')
plt.ylabel('$\delta_{CP}= -0.45 $')
plt.legend()

plt.subplot(322) # deltaCP_1 ;theta23_2 
plt.title('N$_{\u03BD_e}^{FD}$   ${\u03BD_e}$/1GeV (x1.1*10$^{21}$POT)')
N_2 = np.array(osc_mod.Rec2(Acc, E_lisseND, Phi_lisseND * osc_mod.probability_oscillation( theta12 , theta13 , theta23_2,E_lisseND, 1285 ,deltaCP_1*180/np.pi, True, True) * F_lisse * sigma_FD * Npot *Nmultiply ,bias))/Nmultiply
#chi2_2 = osc_mod.chi2( N_2 , N_FD_hist , E_lisseND)
chi2_2 = osc_mod.chi2( N_2 , Nexp1 , E_lisseND)
chi2_str_2=str(round(chi2_2,2))
plt.bar( E_lisseND , N_2,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label=(f"$\chi^2$={chi2_str_2}"))
#plt.plot(E_lisseND, N_FD_hist ,'x',color='red', label='N$_{\u03BD_e}^{FD-data}$')
#plt.errorbar(E_lisseND, N_FD_hist, DN_nu_e, fmt='none',ecolor ='blue')
plt.plot(E_lisseND, Nexp1 ,'x',color='red', label='N$_{\u03BD_e}^{FD-data}$')
plt.errorbar(E_lisseND, Nexp1, DN_nu_e, fmt='none',ecolor ='blue')
plt.legend()

"""
plt.subplot(323) # deltaCP_1 ;theta23_3

N_3 = np.array(osc_mod.Rec2(Acc, E_lisseND, Phi_lisseND * osc_mod.probability_oscillation( theta12 , theta13 , theta23_3 ,E_lisseND, 1285 ,deltaCP_1*180/np.pi, True, True) * F_lisse * sigma_FD * Npot *Nmultiply ,bias))/Nmultiply
chi2_3 = osc_mod.chi2( N_3 , N_FD_hist , E_lisseND)
chi2_str_3=str(round(chi2_3,2))
plt.bar( E_lisseND , N_3,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label=(f"$\chi^2$={chi2_str_3}"))
plt.plot(E_lisseND, N_FD_hist ,'x',color='red', label='N$_{\u03BD_e}^{FD-data}$')
plt.errorbar(E_lisseND, N_FD_hist, DN_nu_e, fmt='none',ecolor ='blue')
plt.legend()
"""
plt.subplot(323) # deltaCP_2 ;theta23_1 

N_4 = np.array(osc_mod.Rec2(Acc, E_lisseND, Phi_lisseND * osc_mod.probability_oscillation( theta12 , theta13 , theta23_1 ,E_lisseND, 1285 ,deltaCP_2*180/np.pi, True, True) * F_lisse * sigma_FD * Npot *Nmultiply ,bias))/Nmultiply
#chi2_4 = osc_mod.chi2( N_4 , N_FD_hist , E_lisseND)
chi2_4 = osc_mod.chi2( N_4 , Nexp1 , E_lisseND)
chi2_str_4=str(round(chi2_4,2))
plt.bar( E_lisseND , N_4,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label=(f"$\chi^2$={chi2_str_4}"))
#plt.plot(E_lisseND, N_FD_hist ,'x',color='red', label='N$_{\u03BD_e}^{FD-data}$')
#plt.errorbar(E_lisseND, N_FD_hist, DN_nu_e, fmt='none',ecolor ='blue')
plt.plot(E_lisseND, Nexp1 ,'x',color='red', label='N$_{\u03BD_e}^{FD-data}$')
plt.errorbar(E_lisseND, Nexp1, DN_nu_e, fmt='none',ecolor ='blue')
plt.legend()
plt.ylabel('$\delta_{CP}= 0 $  (true value)')

plt.subplot(324) # deltaCP_2 ;theta23_2 

N_5 = np.array(osc_mod.Rec2(Acc, E_lisseND, Phi_lisseND * osc_mod.probability_oscillation( theta12 , theta13 , theta23_2 ,E_lisseND, 1285 ,deltaCP_2*180/np.pi, True, True) * F_lisse * sigma_FD * Npot *Nmultiply ,bias))/Nmultiply
#chi2_5 = osc_mod.chi2( N_5 , N_FD_hist , E_lisseND)
chi2_5 = osc_mod.chi2( N_5 , Nexp1 , E_lisseND)
chi2_str_5=str(round(chi2_5,2))
plt.bar( E_lisseND , N_5,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label=(f"$\chi^2$={chi2_str_5}"))
#plt.plot(E_lisseND, N_FD_hist ,'x',color='red', label='N$_{\u03BD_e}^{FD-data}$')
#plt.errorbar(E_lisseND, N_FD_hist, DN_nu_e, fmt='none',ecolor ='blue')
plt.plot(E_lisseND, Nexp1 ,'x',color='red', label='N$_{\u03BD_e}^{FD-data}$')
plt.errorbar(E_lisseND,Nexp1 , DN_nu_e, fmt='none',ecolor ='blue')
plt.legend()
"""
plt.subplot(326) # deltaCP_2 ;theta23_3

N_6 = np.array(osc_mod.Rec2(Acc, E_lisseND, Phi_lisseND * osc_mod.probability_oscillation( theta12 , theta13 , theta23_3 ,E_lisseND, 1285 ,deltaCP_2*180/np.pi, True, True) * F_lisse * sigma_FD * Npot *Nmultiply ,bias))/Nmultiply
chi2_6 = osc_mod.chi2( N_6 , N_FD_hist , E_lisseND)
chi2_str_6=str(round(chi2_6,2))
plt.bar( E_lisseND , N_6,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label=(f"$\chi^2$={chi2_str_6}"))
plt.plot(E_lisseND, N_FD_hist ,'x',color='red', label='N$_{\u03BD_e}^{FD-data}$')
plt.errorbar(E_lisseND, N_FD_hist, DN_nu_e, fmt='none',ecolor ='blue')
plt.legend()
"""
plt.subplot(325) # deltaCP_3 ;theta23_1 

N_7 = np.array(osc_mod.Rec2(Acc, E_lisseND, Phi_lisseND * osc_mod.probability_oscillation( theta12 , theta13 , theta23_1 ,E_lisseND, 1285 ,deltaCP_3*180/np.pi, True, True) * F_lisse * sigma_FD * Npot *Nmultiply ,bias))/Nmultiply
#chi2_7 = osc_mod.chi2( N_7 , N_FD_hist , E_lisseND)
chi2_7 = osc_mod.chi2( N_7 , Nexp1 , E_lisseND)
chi2_str_7=str(round(chi2_7,2))
plt.bar( E_lisseND , N_7,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label=(f"$\chi^2$={chi2_str_7}"))
#plt.plot(E_lisseND, N_FD_hist ,'x',color='red', label='N$_{\u03BD_e}^{FD-data}$')
#plt.errorbar(E_lisseND, N_FD_hist, DN_nu_e, fmt='none',ecolor ='blue')
plt.plot(E_lisseND, Nexp1 ,'x',color='red', label='N$_{\u03BD_e}^{FD-data}$')
plt.errorbar(E_lisseND, Nexp1, DN_nu_e, fmt='none',ecolor ='blue')
plt.xlabel('$\Theta_{23} = 0.80 rad$')
plt.ylabel('$\delta_{CP}= -0.35 (min value)$')
plt.legend()

plt.subplot(326) # deltaCP_3 ;theta23_2 

N_8 = np.array(osc_mod.Rec2(Acc, E_lisseND, Phi_lisseND * osc_mod.probability_oscillation( theta12 , theta13 , theta23_2 ,E_lisseND, 1285 ,deltaCP_3*180/np.pi, True, True) * F_lisse * sigma_FD * Npot *Nmultiply ,bias))/Nmultiply
#chi2_8 = osc_mod.chi2( N_8 , N_FD_hist , E_lisseND)
chi2_8 = osc_mod.chi2( N_8 , Nexp1 , E_lisseND)
chi2_str_8=str(round(chi2_8,2))
plt.bar( E_lisseND , N_8,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label=(f"$\chi^2$={chi2_str_8}"))
#plt.plot(E_lisseND, N_FD_hist ,'x',color='red', label='N$_{\u03BD_e}^{FD-data}$')
#plt.errorbar(E_lisseND, N_FD_hist, DN_nu_e, fmt='none',ecolor ='blue')
plt.plot(E_lisseND, Nexp1 ,'x',color='red', label='N$_{\u03BD_e}^{FD-data}$')
plt.errorbar(E_lisseND, Nexp1 , DN_nu_e, fmt='none',ecolor ='blue')
plt.xlabel('$\Theta_{23} = 0.87 (true/min value) rad$')
plt.legend()
"""
plt.subplot(329) # deltaCP_3 ;theta23_3 

N_9 = np.array(osc_mod.Rec2(Acc, E_lisseND, Phi_lisseND * osc_mod.probability_oscillation( theta12 , theta13 , theta23_3 ,E_lisseND, 1285 ,deltaCP_3*180/np.pi, True, True) * F_lisse * sigma_FD * Npot *Nmultiply ,bias))/Nmultiply
chi2_9 = osc_mod.chi2( N_9 , N_FD_hist , E_lisseND)
chi2_str_9=str(round(chi2_9,2))
plt.bar( E_lisseND , N_9,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label=(f"$\chi^2$={chi2_str_9}"))
plt.plot(E_lisseND, N_FD_hist ,'x',color='red', label='N$_{\u03BD_e}^{FD-data}$')
plt.errorbar(E_lisseND, N_FD_hist, DN_nu_e, fmt='none',ecolor ='blue')
plt.xlabel('$\Theta_{23} = 0.87 rad$  (min value)')
plt.legend()
"""
plt.show()
plt.close()

#Plot of the true values energy distributon spectrum

plt.bar( E_lisseND , N_5,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label=(f"$\chi^2$={chi2_str_5}"))
#plt.plot(E_lisseND, N_FD_hist ,'x',color='red', label='N$_{\u03BD_e}^{FD-data}$')
#plt.errorbar(E_lisseND, N_FD_hist, DN_nu_e, fmt='none',ecolor ='blue')
plt.plot(E_lisseND, Nexp1 ,'x',color='red', label='N$_{\u03BD_e}^{FD-data}$')
plt.errorbar(E_lisseND, Nexp1, DN_nu_e, fmt='none',ecolor ='blue')
plt.title('N$_{\u03BD_e}^{FD}$   for $\Theta_{23}$ = 0.87 rad and $\delta_{CP}$= 0')
plt.xlabel('E (GeV)')
plt.ylabel('N$_{\u03BD_e}$    ${\u03BD_e}$/1GeV (x1.1*10$^{21}$POT) ')
plt.legend()
plt.show()
plt.close()



# IX ) Khi square
# IX) a) Khi square function of delta CP


delta_CP = np.linspace(-np.pi,np.pi,36)
L_Papp = []
L_N_nue = []
L_Int = []
khi2=[]

for i in range (0,len(delta_CP)):
	L_Papp = L_Papp + [osc_mod.probability_oscillation(theta12 , theta13 , theta23, E_lisseND, 1285 ,delta_CP[i]*180/np.pi, True, True)]
	
	
	L_N_nue = L_N_nue + [ np.array(osc_mod.Rec2(Acc, E_lisseND, Phi_lisseND * osc_mod.probability_oscillation( theta12 , theta13 , theta23,E_lisseND, 1285 ,delta_CP[i]*180/np.pi, True, True) * F_lisse * sigma_FD * Npot *Nmultiply ,bias))/Nmultiply ]
	
	L_N_nue[i] = L_N_nue[i] 	
	khi2 = khi2 + [0]
	for j in range(0,len(E_lisseND)):
		if L_N_nue[i][j] != 0:
		#	khi2[i] = khi2[i] + ( L_N_nue[i][j] - N_FD_hist[j] )**2/(4*L_N_nue[i][j])
                        khi2[i] = khi2[i] + ( L_N_nue[i][j] - Nexp1[j] )**2/(4*L_N_nue[i][j])
		else :
			khi2[i] = khi2[i]

#plt.bar(E_lisseND , L_N_nue[75] ,color='grey',edgecolor = 'black',width=(7-0.5)/Nbin,label='$\phi_{\u03BD_{\mu}}^{ND}$ * P$_{\u03BD_{\mu} -> \u03BD_e}$ * F * $\sigma_{\u03BD_e}^{Ar}$ *T')
plt.plot( delta_CP , khi2 ,label='$\chi^2$')
plt.ylabel('$\chi^2$')
plt.yscale('log')
plt.xlabel('$\delta_{CP}$  radians')
plt.legend()


#plt.plot(ene_FD_article, np.array(spec_FD_article) ,'x',color='red', label='N$_{\u03BD_e}^{FD-data}$')
plt.show()	
plt.close()



# IX) b) khi 2 function of sin^2(theta23) (The renormalisation factor is constant)

L_theta23 = np.linspace(-3.14,3.14,72)
L_sin23squar = np.sin(L_theta23)**2
L_Papp = []
L_N_nue = []
L_Int = []
khi2=[]

for i in range (0,len(L_theta23)):
	L_Papp = L_Papp + [osc_mod.probability_oscillation(theta12 , theta13 , L_theta23[i], E_lisseND, 1285 ,0*180/np.pi, True, True)]
	
	
	L_N_nue = L_N_nue + [ np.array(osc_mod.Rec2(Acc, E_lisseND, Phi_lisseND * osc_mod.probability_oscillation( theta12 , theta13 , L_theta23[i],E_lisseND, 1285 ,0*180/np.pi, True, True) * F_lisse * sigma_FD * Npot *Nmultiply ,bias))/Nmultiply ]
	
	
	L_N_nue[i] = L_N_nue[i] 
	khi2 = khi2 + [0]
	for j in range(0,len(E_lisseND)):
		if L_N_nue[i][j] != 0:
		#	khi2[i] = khi2[i] + ( L_N_nue[i][j] - N_FD_hist[j] )**2/(4*L_N_nue[i][j])
			khi2[i] = khi2[i] + ( L_N_nue[i][j] - Nexp1[j] )**2/(4*L_N_nue[i][j])
		else :
			khi2[i] = khi2[i]



plt.plot( L_theta23 , khi2 , label='$\chi^2$')
plt.xlabel('$\Theta_{23}$')
plt.ylabel('$\chi^2$')
plt.yscale('log')
plt.legend()
plt.show()
plt.close()

plt.plot( L_sin23squar , khi2 , label='$\chi^2$')
plt.xlabel('$sin(\Theta_{23})^2$')
plt.ylabel('$\chi^2$')
plt.yscale('log')
plt.legend()
plt.show()




#plt.plot(L_theta23,L_sin23squar)
plt.show()

plt.plot( L_sin23squar , np.sqrt(khi2-min(khi2)),label='$\delta$$\chi^2$')
plt.ylabel('$\sqrt{\Delta\chi^2}$')
plt.xlabel('$sin(\Theta_{23})^2$')
plt.legend()

plt.show()

# IX) c) Chi2 as a function of sin2(theta23) and deltaCP


L_deltaCP_o_pi = np.linspace (-1,1,100) 
L_theta23 = np.linspace(0.3,1.3,100)
L_sin23 = np.sin(L_theta23)**2
L2_Int = []
LL2_N_nue = []
L2_khi2 = []


for i in range (0,len(L_theta23)):
	

	LL2_N_nue = LL2_N_nue +[[]]
	L2_khi2 = L2_khi2 + [[]]
	L2_Int = L2_Int +[[]]

	for j in range (0,len(L_deltaCP_o_pi)) :
		LL2_N_nue[i] = LL2_N_nue[i] + [np.array(osc_mod.Rec2( Acc,E_lisseND, Phi_lisseND * osc_mod.probability_oscillation_v2(osc_par.delta_m31, osc_par.delta_m32 ,theta12 , theta13, L_theta23[i], E_lisseND, 1285 ,L_deltaCP_o_pi[j]*np.pi*180/np.pi, True, True) * F_lisse * sigma_FD * Npot *Nmultiply , bias))/Nmultiply]
	
		LL2_N_nue[i][j] =  LL2_N_nue[i][j] 	
		L2_khi2[i] = L2_khi2[i] + [0]
		for k in range (0,len(E_lisseND)):
			if LL2_N_nue[i][j][k] != 0:
			#	L2_khi2[i][j] = L2_khi2[i][j] + ( LL2_N_nue[i][j][k] - N_FD_hist[k] )**2/(4*LL2_N_nue[i][j][k])
	                        L2_khi2[i][j] = L2_khi2[i][j] + ( LL2_N_nue[i][j][k] - Nexp1[k] )**2/(4*LL2_N_nue[i][j][k])
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
DCPoPI , STheta23  = np.meshgrid (L_deltaCP_o_pi , L_sin23)
fig, ax = plt.subplots()
origin = 'lower'
CS2 = ax.contour(DCPoPI , STheta23 , L2_khi2, levels=[chi2min + Dchi2_1s , chi2min + Dchi2_3s ,chi2min + Dchi2_5s], colors=['green', 'pink', 'grey'],norm='log', origin=origin) #Make the contour at 1sigma , 3sigma and 5 sgima 
pc = ax.pcolormesh( DCPoPI, STheta23 , L2_khi2, norm='log')
fig.colorbar(pc,ax=ax,label='$\chi^2$')
plt.title('$\epsilon = 50 kT.MW.yrs$')
plt.xlabel('$\delta_{CP} / \pi $ ')
plt.ylabel('$sin(\Theta_{23})^2$ ')
plt.plot(0 ,np.sin(0.86)**2,'o',color='red')
chi2min_str = str(round( L2_khi2[osc_mod.ind_min(L2_khi2)[0]][osc_mod.ind_min(L2_khi2)[1]] , 2))
chi2min_dof_str = str(round( L2_khi2[osc_mod.ind_min(L2_khi2)[0]][osc_mod.ind_min(L2_khi2)[1]] / (Nbin-2) , 2))
chi2min = round( L2_khi2[osc_mod.ind_min(L2_khi2)[0]][osc_mod.ind_min(L2_khi2)[1]] , 1)
plt.plot( L_deltaCP_o_pi[osc_mod.ind_min(L2_khi2)[1]], L_sin23[osc_mod.ind_min(L2_khi2)[0]] ,'x',color='red',label=(f"$\chi^2/Ndof$={chi2min_dof_str}"))
plt.legend(legend, legend_label, loc ='lower right')
plt.show()
#print(L2_khi2)


# IX) d) Chi2 function of the accuracy of the reconstruction of the energy
"""

L_accuracy = np.linspace(14,30,70)
L_N_nue = []

khi2=[]

for i in range (0,len(L_accuracy)):

	L_N_nue = L_N_nue + [ np.array(osc_mod.Rec2(L_accuracy[i], E_lisseND, Phi_lisseND * osc_mod.probability_oscillation( osc_par.theta12 , osc_par.theta13 , osc_par.theta23,E_lisseND, 1285 ,0*180/np.pi, True, True) * F_lisse * sigma_FD * Npot *Nmultiply ,bias))/Nmultiply ]
	
	khi2 = khi2 + [0]
	for j in range(0,len(E_lisseND)):
		if L_N_nue[i][j] != 0:
		#	khi2[i] = khi2[i] + ( L_N_nue[i][j] - N_FD_hist[j] )**2/(4*L_N_nue[i][j])
		        khi2[i] = khi2[i] + ( L_N_nue[i][j] - Nexp1[j] )**2/(4*L_N_nue[i][j])
		else :
			khi2[i] = khi2[i]
			
plt.plot( L_accuracy , khi2 , label='$\chi^2$')
plt.xlabel('accuracy in %')
plt.ylabel('$\chi^2$')
#plt.yscale('log')
plt.legend()
plt.show()
plt.close()
"""
# IX) d)2) distribution function of the accuracy of the reconstruction of the energy

"""
Acc_1 = 20 
deltaCP_min_acc1 = -0.61
theta23_min_acc1 = 0.82
Acc_2 = 21
deltaCP_min_acc2 = -0.61
theta23_min_acc2 = 0.82
Acc_3 = 22 
deltaCP_min_acc3 = -0.52 
theta23_min_acc3 = 0.83

plt.subplot(321)

plt.title('min values')
N_1 = np.array(osc_mod.Rec2(Acc_1, E_lisseND, Phi_lisseND * osc_mod.probability_oscillation( theta12 , theta13 , theta23_min_acc1 ,E_lisseND, 1285 ,deltaCP_min_acc1*180/np.pi, True, True) * F_lisse * sigma_FD * Npot *Nmultiply ,bias))/Nmultiply
#chi2_1 = osc_mod.chi2( N_1 , N_FD_hist , E_lisseND)
chi2_1 = osc_mod.chi2( N_1 , Nexp1 , E_lisseND)
chi2_str_1=str(round(chi2_1,2))
plt.bar(E_lisseND , N_1 ,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label=(f"$\chi^2$={chi2_str_1}"))
#plt.plot(E_lisseND, N_FD_hist ,'x',color='red', label='N$_{\u03BD_e}^{FD-data}$')
#plt.errorbar(E_lisseND, N_FD_hist, DN_nu_e, fmt='none',ecolor ='blue')
plt.plot(E_lisseND, Nexp1 ,'x',color='red', label='N$_{\u03BD_e}^{FD-data}$')
plt.errorbar(E_lisseND, Nexp1, DN_nu_e, fmt='none',ecolor ='blue')
plt.ylabel('Accuracy = 20%')
#plt.xlabel('$\Theta23 = 0.86 $')
plt.xlabel('$\delta_{CP}= -0.61 ; \Theta_{23} = 0.86  $')
plt.legend()

plt.subplot(322)

plt.title('true values : $\delta_{CP}= 0 ; \Theta_{23} = 0.86  $')
N_2 = np.array(osc_mod.Rec2(Acc_1, E_lisseND, Phi_lisseND * osc_mod.probability_oscillation( theta12 , theta13 , theta23 ,E_lisseND, 1285 ,0*180/np.pi, True, True) * F_lisse * sigma_FD * Npot *Nmultiply ,bias))/Nmultiply
#chi2_2 = osc_mod.chi2( N_2 , N_FD_hist , E_lisseND)
chi2_2 = osc_mod.chi2( N_2 , Nexp1 , E_lisseND)
chi2_str_2=str(round(chi2_2,2))
plt.bar( E_lisseND , N_2,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label=(f"$\chi^2$={chi2_str_2}"))
#plt.plot(E_lisseND, N_FD_hist ,'x',color='red', label='N$_{\u03BD_e}^{FD-data}$')
#plt.errorbar(E_lisseND, N_FD_hist, DN_nu_e, fmt='none',ecolor ='blue')
plt.plot(E_lisseND, Nexp1 ,'x',color='red', label='N$_{\u03BD_e}^{FD-data}$')
plt.errorbar(E_lisseND, Nexp1, DN_nu_e, fmt='none',ecolor ='blue')
plt.legend()

plt.subplot(323)


N_3 = np.array(osc_mod.Rec2(Acc_2, E_lisseND, Phi_lisseND * osc_mod.probability_oscillation( theta12 , theta13 , theta23_min_acc2 ,E_lisseND, 1285 ,deltaCP_min_acc2*180/np.pi, True, True) * F_lisse * sigma_FD * Npot *Nmultiply ,bias))/Nmultiply
#chi2_3 = osc_mod.chi2( N_3 , N_FD_hist , E_lisseND)
chi2_3 = osc_mod.chi2( N_3 , Nexp1 , E_lisseND)
chi2_str_3=str(round(chi2_3,2))
plt.bar(E_lisseND , N_3 ,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label=(f"$\chi^2$={chi2_str_3}"))
#plt.plot(E_lisseND, N_FD_hist ,'x',color='red', label='N$_{\u03BD_e}^{FD-data}$')
#plt.errorbar(E_lisseND, N_FD_hist, DN_nu_e, fmt='none',ecolor ='blue')
plt.plot(E_lisseND, Nexp1 ,'x',color='red', label='N$_{\u03BD_e}^{FD-data}$')
plt.errorbar(E_lisseND, Nexp1, DN_nu_e, fmt='none',ecolor ='blue')
#plt.ylabel('$\delta_{CP}= -0.61 $  (min value for acc=0.2)')
#plt.xlabel('$\Theta23 = 0.86 $')
plt.ylabel('Accuracy = 21%')
plt.xlabel('$\delta_{CP}= -0.61 ; \Theta_{23} = 0.86  $')
plt.legend()

plt.subplot(324)


N_4 = np.array(osc_mod.Rec2(Acc_2, E_lisseND, Phi_lisseND * osc_mod.probability_oscillation( theta12 , theta13 , theta23 ,E_lisseND, 1285 ,0*180/np.pi, True, True) * F_lisse * sigma_FD * Npot *Nmultiply ,bias))/Nmultiply
#chi2_4 = osc_mod.chi2( N_4 , N_FD_hist , E_lisseND)
chi2_4 = osc_mod.chi2( N_4 , Nexp1 , E_lisseND)
chi2_str_4=str(round(chi2_4,2))
plt.bar( E_lisseND , N_4,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label=(f"$\chi^2$={chi2_str_4}"))
#plt.plot(E_lisseND, N_FD_hist ,'x',color='red', label='N$_{\u03BD_e}^{FD-data}$')
#plt.errorbar(E_lisseND, N_FD_hist, DN_nu_e, fmt='none',ecolor ='blue')
plt.plot(E_lisseND, Nexp1 ,'x',color='red', label='N$_{\u03BD_e}^{FD-data}$')
plt.errorbar(E_lisseND, Nexp1, DN_nu_e, fmt='none',ecolor ='blue')
plt.legend()

plt.subplot(325)


N_5 = np.array(osc_mod.Rec2(Acc_3, E_lisseND, Phi_lisseND * osc_mod.probability_oscillation( theta12 , theta13 , theta23_min_acc3 ,E_lisseND, 1285 ,deltaCP_min_acc3*180/np.pi, True, True) * F_lisse * sigma_FD * Npot *Nmultiply ,bias))/Nmultiply
#chi2_5 = osc_mod.chi2( N_5 , N_FD_hist , E_lisseND)
chi2_5 = osc_mod.chi2( N_5 , Nexp1 , E_lisseND)
chi2_str_5=str(round(chi2_5,2))
plt.bar(E_lisseND , N_5 ,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label=(f"$\chi^2$={chi2_str_5}"))
#plt.plot(E_lisseND, N_FD_hist ,'x',color='red', label='N$_{\u03BD_e}^{FD-data}$')
#plt.errorbar(E_lisseND, N_FD_hist, DN_nu_e, fmt='none',ecolor ='blue')
plt.plot(E_lisseND, Nexp1 ,'x',color='red', label='N$_{\u03BD_e}^{FD-data}$')
plt.errorbar(E_lisseND, Nexp1, DN_nu_e, fmt='none',ecolor ='blue')
#plt.ylabel('$\delta_{CP}= -0.52 $  (min value for acc=0.2)')
#plt.xlabel('$\Theta23 = 0.83 $')
plt.ylabel('Accuracy = 22%')
plt.xlabel('$\delta_{CP}= -0.52 ; \Theta_{23} = 0.83  $')
plt.legend()

plt.subplot(326)


N_6 = np.array(osc_mod.Rec2(Acc_3, E_lisseND, Phi_lisseND * osc_mod.probability_oscillation( theta12 , theta13 , theta23 ,E_lisseND, 1285 ,0*180/np.pi, True, True) * F_lisse * sigma_FD * Npot *Nmultiply ,bias))/Nmultiply
#chi2_6 = osc_mod.chi2( N_6 , N_FD_hist , E_lisseND)
chi2_6 = osc_mod.chi2( N_6 , Nexp1 , E_lisseND)
chi2_str_6=str(round(chi2_6,2))
plt.bar( E_lisseND , N_6,color='grey',edgecolor = 'black',width=(Emax-0.5)/Nbin,label=(f"$\chi^2$={chi2_str_6}"))
#plt.plot(E_lisseND, N_FD_hist ,'x',color='red', label='N$_{\u03BD_e}^{FD-data}$')
#plt.errorbar(E_lisseND, N_FD_hist, DN_nu_e, fmt='none',ecolor ='blue')
plt.plot(E_lisseND, Nexp1 ,'x',color='red', label='N$_{\u03BD_e}^{FD-data}$')
plt.errorbar(E_lisseND, Nexp1, DN_nu_e, fmt='none',ecolor ='blue')
plt.legend()
plt.show()
"""
# IX) e) Chi2 function of the bias of the reconstruction of the energy
"""
Lbias = np.linspace(-10,15,50)
L_Ene_bias = []
L_N_nue = []

khi2=[]

for i in range (0,len(Lbias)):

	L_N_nue = L_N_nue + [ np.array(osc_mod.Rec2(Acc, E_lisseND, Phi_lisseND * osc_mod.probability_oscillation ( theta12 , theta13 , theta23 ,E_lisseND, 1285 ,0*180/np.pi, True, True )  * F_lisse * sigma_FD * Npot *Nmultiply , Lbias[i]))/Nmultiply ]

	
	khi2 = khi2 + [0]
	for j in range(0,len(E_lisseND)):
		if L_N_nue[i][j] != 0:
#			khi2[i] = khi2[i] + ( L_N_nue[i][j] - N_FD_hist[j] )**2/(4*L_N_nue[i][j])
                        khi2[i] = khi2[i] + ( L_N_nue[i][j] - Nexp1[j] )**2/(4*L_N_nue[i][j])
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
"""
#True values for theta23 and deltam23

L_bias = np.linspace(-10,10,100)
L_Reco= np.linspace(15,27,100)
LL2_N_nue = []
L2_khi2 = []


for i in range (0,len(L_bias)):
	

	LL2_N_nue = LL2_N_nue +[[]]
	L2_khi2 = L2_khi2 + [[]]

	for j in range (0,len(L_Reco)) :
		LL2_N_nue[i] = LL2_N_nue[i] + [np.array(osc_mod.Rec2(L_Reco[j], E_lisseND, Phi_lisseND * osc_mod.probability_oscillation ( theta12 , theta13 , theta23 ,E_lisseND, 1285 ,0*180/np.pi, True, True )  * F_lisse * sigma_FD * Npot *Nmultiply , L_bias[i]))/Nmultiply]
				
		L2_khi2[i] = L2_khi2[i] + [0]
		for k in range (0,len(E_lisseND)):
			if LL2_N_nue[i][j][k] != 0:
#				L2_khi2[i][j] = L2_khi2[i][j] + ( LL2_N_nue[i][j][k] - N_FD_hist[k] )**2/(4*LL2_N_nue[i][j][k])
                                L2_khi2[i][j] = L2_khi2[i][j] + ( LL2_N_nue[i][j][k] - Nexp1[k] )**2/(4*LL2_N_nue[i][j][k])
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
"""
