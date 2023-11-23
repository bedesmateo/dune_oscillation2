# -*- coding: utf-8 -*-
"""
Created on 4th of Nov 2022
Thibaut Houdy

Last edited 31th of March 2023
Matéo Bédès

"""
import matplotlib.pyplot as plt
import numpy as np
import random as rand
from random import gauss
import params as osc_par

Nbin=25
Emax_numu = 4
Emax=7

#cumulative function for measured energy spectrum 
def cumulative(y):
    cumul, f_cumul = 0, []
    f_cumul.append(0)
    for yy in y:
        cumul += yy
        f_cumul.append(cumul)

    return f_cumul







#generate a gaussian energy spectrum for near or far detector

def mc_function(x_ene, f_y, Ntrials, sigma): 
    new_spectrum =[]
    for i in range(int(Ntrials)):
        t_rand = rand.uniform(0,1)
        max_bin = max(np.where(t_rand>f_y)[0])
        try : 
            true_ene = x_ene[max_bin]
            rand_energie = rand.gauss(true_ene,sigma)
            if rand_energie > 0.5 and rand_energie < Emax :
                new_spectrum.append(rand_energie)
        except : 
            continue
    return new_spectrum



#simulate neutrino spectrum using measured energy spectrum
def model_flux(x,y,N,sigma):
    y_cumul = cumulative(y)
    y_cumul = y_cumul/max(y_cumul)
    flux = np.array(mc_function(x,y_cumul,N,sigma))
    return flux
    

#  calculate oscillation probability (Normal Ordering)  [E]-GeV, [L]-km, deltaCP & a
"""
def probability_oscillation(E=[2.5], L=1000, dCP=osc_par.delta_CP, b_normal_hierarchy=True, b_neutrino=True  ):   

    c1 = np.sin(osc_par.theta23)**2*np.sin(2*osc_par.theta13)**2
    c2 = np.sin(2*osc_par.theta23)*np.sin(2*osc_par.theta13)*np.sin(2*osc_par.theta12)
    c3 = np.cos(osc_par.theta23)**2*np.sin(2*osc_par.theta12)**2

    #osc_par.normal_hierarchy(b_normal_hierarchy)
    #osc_par.neutrinos(b_neutrino)

    if b_normal_hierarchy:
        delta_m3l = osc_par.delta_m31
    else:
        delta_m3l = osc_par.delta_m32

    if b_neutrino :
        a = osc_par.a_matter
    else:
        a=-1*osc_par.a_matter
        dCP=-dCP

    def D3l(E):
        #rint(E)
        return(1.267*delta_m3l*L/E)

    def D21(E):
        return(1.267*osc_par.delta_m21*L/E)

    if dCP<-180 or dCP>180 : 
        dCP = osc_par.delta_CP
   
    Proba = c1 * np.sin(D3l(E)-a*L)**2/(D3l(E)-a*L)**2 * D3l(E)**2 \
                  + c2 *np.sin(D3l(E)-a*L)/(D3l(E)-a*L)*D3l(E)*np.sin(a*L)/a/L \
                    * D21(E)*np.cos(D3l(E)+np.radians(dCP)) + c3*D21(E)**2*np.sin(a*L)**2/(a*L)**2
    
    return Proba
"""
    
    
def probability_oscillation(theta12 , theta13 , theta23,E=[2.5], L=1000, dCP=osc_par.delta_CP, b_normal_hierarchy=True, b_neutrino=True):   

    c1 = np.sin(theta23)**2*np.sin(2*theta13)**2
    c2 = np.sin(2*theta23)*np.sin(2*theta13)*np.sin(2*theta12)
    c3 = np.cos(theta23)**2*np.sin(2*theta12)**2

    #osc_par.normal_hierarchy(b_normal_hierarchy)
    #osc_par.neutrinos(b_neutrino)

    if b_normal_hierarchy:
        delta_m3l = osc_par.delta_m31
    else:
        delta_m3l = osc_par.delta_m32

    if b_neutrino :
        a = osc_par.a_matter
    else:
        a=-1*osc_par.a_matter
        dCP=-dCP

    def D3l(E):
        #rint(E)
        return(1.267*delta_m3l*L/E)

    def D21(E):
        return(1.267*osc_par.delta_m21*L/E)

    if dCP<-180 or dCP>180 : 
        dCP = osc_par.delta_CP
   
    Proba = c1 * np.sin(D3l(E)-a*L)**2/(D3l(E)-a*L)**2 * D3l(E)**2 \
                  + c2 *np.sin(D3l(E)-a*L)/(D3l(E)-a*L)*D3l(E)*np.sin(a*L)/a/L \
                    * D21(E)*np.cos(D3l(E)+np.radians(dCP)) + c3*D21(E)**2*np.sin(a*L)**2/(a*L)**2
    
    return Proba
    
    
def probability_oscillation_v2(delta_m31, delta_m32, theta12 , theta13 , theta23,E=[2.5], L=1000, dCP=osc_par.delta_CP, b_normal_hierarchy=True, b_neutrino=True):   

    c1 = np.sin(theta23)**2*np.sin(2*theta13)**2
    c2 = np.sin(2*theta23)*np.sin(2*theta13)*np.sin(2*theta12)
    c3 = np.cos(theta23)**2*np.sin(2*theta12)**2

    #osc_par.normal_hierarchy(b_normal_hierarchy)
    #osc_par.neutrinos(b_neutrino)

    if b_normal_hierarchy:
        delta_m3l = delta_m31
    else:
        delta_m3l = delta_m32

    if b_neutrino :
        a = osc_par.a_matter
    else:
        a=-1*osc_par.a_matter
        dCP=-dCP

    def D3l(E):
        #rint(E)
        return(1.267*delta_m3l*L/E)

    def D21(E):
        return(1.267*osc_par.delta_m21*L/E)

    if dCP<-180 or dCP>180 : 
        dCP = osc_par.delta_CP
   
    Proba = c1 * np.sin(D3l(E)-a*L)**2/(D3l(E)-a*L)**2 * D3l(E)**2 \
                  + c2 *np.sin(D3l(E)-a*L)/(D3l(E)-a*L)*D3l(E)*np.sin(a*L)/a/L \
                    * D21(E)*np.cos(D3l(E)+np.radians(dCP)) + c3*D21(E)**2*np.sin(a*L)**2/(a*L)**2
    
    return Proba    



#If cos(theta13) ~= 1
def proba_disapp_0 ( theta23, E , L , b_normal_hierarchy=True ):

	if b_normal_hierarchy:
		delta_matm = osc_par.delta_m32
	else : 
		delta_matm = osc_par.delta_m31
	def Datm(E):
		return(1.267*delta_matm*L/E)
		
	Proba = 1 - (np.sin(2*theta23)**2) * (np.sin(Datm(E)))**2
	
	return Proba

def proba_disapp_0_t ( delta_m32 ,theta23, E , L , b_normal_hierarchy=True ):

	if b_normal_hierarchy:
		delta_matm = delta_m32
	else : 
		delta_matm = osc_par.delta_m31
	def Datm(E):
		return(1.267*delta_matm*L/E)
		
	Proba = 1 - (np.sin(2*theta23)**2) * (np.sin(Datm(E)))**2
	
	return Proba	

theta23 = np.radians(49.8)	
Ex = np.linspace(0.5,7,1000)
plt.plot( Ex, proba_disapp_0 ( theta23, Ex , 1285 , True ) ,label='P$_{\u03BD_{\mu} -> \u03BD_{\mu}}$')	
plt.xlabel('E GeV')
plt.ylabel('P$_{\u03BD_{\mu} -> \u03BD_{\mu}}$')
plt.legend()
plt.show()
plt.close()
	
	

#Else
def proba_disapp_1 (  theta13 ,  theta23, E , L , b_normal_hierarchy=True ):

	if b_normal_hierarchy:
		delta_matm = osc_par.delta_m32
	else : 
		delta_matm = osc_par.delta_m31
	def Datm(E):
		return(1.267*delta_matm*L/E)
		
	Proba =  1 - 4 * np.cos(theta13)**2 * np.sin(theta23)**2 * ( 1 - np.cos(theta13)**2 * np.sin(theta23)**2 ) * (np.sin(Datm(E)))**2 
	
	return Proba
	
	
	


#One defines a first function that has as argumets two scalar. 1 energy value and 1 Nevents value corresponding to this energy. T_rec returnes a Gaussian distribution of eprance E and standard deviation 0.2*E.
	
#Final reconstruction function, the good one !
def G ( x , x0 , sigma ):
	return 1 / ( sigma * np.sqrt(2*np.pi) ) * np.exp( -1/2 * ((x-x0)/sigma)**2  )


E_lisseND = np.linspace(0.5, Emax,Nbin)
def Rec2 (accuracy,E,N,bias):
	O=[]
	for k in range ( 0 , len(E)):
		O=O+[0]
		for i in range (0, len(E)):
			O[k] = O[k] + (G( E_lisseND + bias/100*E[i] , E[i] , accuracy/100*E[i] ) * N[i]* (Emax-0.5) / Nbin)[k]
	return O

E_lisseND_numu = np.linspace(0.5, Emax_numu,Nbin)	
def Rec2_numu (accuracy,E,N,bias):
	O=[]
	for k in range ( 0 , len(E)):
		O=O+[0]
		for i in range (0, len(E)):
			O[k] = O[k] + (G( E_lisseND_numu + bias/100*E[i], E[i] , accuracy/100*E[i] ) * N[i]* (Emax_numu-0.5) / Nbin)[k]
	return O



def chi2 (Nmod , Ndata, Liste_biningE):
	chi2 = 0
	for i in range(0,len(Liste_biningE)):
		chi2 = chi2 + ( Nmod[i] - Ndata[i] )**2 / (4*Nmod[i])
	return chi2
	
def ind_min ( chi ):
	M=[]
	I=[]
	for i in range (0,len(chi)):
		M=M+[min(chi[i])]
		I=I+[chi[i].index(min(chi[i]))]
	return [ M.index(min(M)) , I[M.index(min(M))] ]


def chi2min_1sigma (chi2):
	L=[]
	A=[]
	B=[]
	for i in range (0,len(chi2)):
		for j in range (0,len(chi2[i])):
			if chi2[i][j]== round( chi2[ind_min(chi2)[0]][ind_min(chi2)[1]] , 0) + 2:
				L=L+[[i,j]]
			else :
				L=L
	return L


def chi2min_2sigma (chi2):
	L=[]
	A=[]
	B=[]
	for i in range (0,len(chi2)):
		for j in range (0,len(chi2[i])):
			if chi2[i][j]== round( chi2[ind_min(chi2)[0]][ind_min(chi2)[1]] , 0) + 6:
				L=L+[[i,j]]
			else :
				L=L
	return L

def chi2min_3sigma (chi2):
	L=[]
	A=[]
	B=[]
	for i in range (0,len(chi2)):
		for j in range (0,len(chi2[i])):
			if chi2[i][j]== round( chi2[ind_min(chi2)[0]][ind_min(chi2)[1]] , 0) + 12:
				L=L+[[i,j]]
			else :
				L=L
	return L




	
