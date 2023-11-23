import matplotlib.pyplot as plt
import numpy as np
import random as rand
from random import gauss

def gamma (n):
	g = 1
	fact1 =1
	fact2= 1
	
	if (2*n)%2 == 0 :
		for i in range (1,int(n)):
			g=g*i
		return g
	elif (2*n)%2 == 1 :
		for i in range (1,int(2*n-1)+1):
			fact1=fact1*i
		for i in range (1,int(n-1/2)+1):
			fact2=fact2*i
		if n==1/2 :
			return np.sqrt(np.pi)
		else :
			return fact1 / ( pow(2 , 2*n-1) * fact2) * np.sqrt(np.pi)

#Probability density function of chi2 for a given degrees of freedom

def fchi2 (x,k):
	return pow(  1/2 , k/2 ) / gamma(k/2) * pow( x , k/2-1 ) * np.exp(-x/2)


def f5 (x) :
	return pow(1/2,5/2) / (3*np.sqrt(np.pi)/4) * pow(x,(5/2)-1) * np.exp(-x/2)


#Return the critical chi2 for a given risk, and for a given degrees of freedom

def chi_crit( risque , k ,chi2a,chi2b):
	Int=0
	Chi_c=0
	x=np.linspace(chi2a,chi2b,10000) #the integral of the probability density (fchi2) has to be 1 on the interval [chi2a,chi2b]
	for i in range (0,len(x)):
		if Int < 1-risque :
			Int = Int + fchi2(x[i],k)*(chi2b-chi2a)/10000
			Chi_c = x[i]	
		else :
			Int = Int
	return Chi_c


#Return the p_value of an adjustment for a given value of chi2, and for a given degrees of freedom

def p_value (chi2 , k ,chi2a,chi2b):
	Int=0
	x=np.linspace(chi2a,chi2b,10000) #the integral of the probability density (fchi2) has to be 1 on the interval [chi2a,chi2b]
	ind_chi2 = int(chi2/(chi2b-chi2a)*10000)
	for i in range (ind_chi2,len(x)):
		Int = Int + fchi2(x[i],k)*(chi2b-chi2a)/10000
		Int = round(Int,20)
	return Int

print(p_value (25 , 1 ,0,110) )
print(chi_crit(0.05,2,0,100))
