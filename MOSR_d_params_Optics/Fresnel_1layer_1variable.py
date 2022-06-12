# -*- coding: utf-8 -*-
"""
Created on Thu May 19 16:56:59 2022

@author: mrsri
"""
import numpy as np
import matplotlib
import scipy
import sympy
import matplotlib.pyplot as plt
import pandas as pd
"Clearing all the previous varaibles's values"
#from IPython import get_ipython
#get_ipython().magic('reset -sf')

"Définition des fonctions"
def angle(ni,ne,thetain):
    exit_angle=np.arcsin((ni/ne)*np.sin(thetain))
    return exit_angle

def phaseshift(d,ne,e,exit_angle):
    shift=(1j*2*np.pi*d*ne*np.cos(exit_angle))
    return shift

def refle_coeff(ni,ne,in_angle,exit_angle):
    if (pol=='s'):
        coeff_refl=((ni*np.cos(in_angle)-ne*np.cos(exit_angle))/(ni*np.cos(in_angle)+ne*np.cos(exit_angle)))
    else:
        coeff_refl=((ne*np.cos(in_angle)-ni*np.cos(exit_angle))/(ne*np.cos(in_angle)+ni*np.cos(exit_angle)))
    return coeff_refl
def trans_coeff(ni,ne,in_angle,exit_angle):
    if (pol=='s'):
        coeff_trans=((2*ni*np.cos(in_angle))/(ni*np.cos(in_angle)+ne*np.cos(exit_angle)))
    else:
        coeff_trans=((2*ni*np.cos(in_angle))/(ne*np.cos(in_angle)+ni*np.cos(exit_angle)))
    return coeff_trans

"Parameters simulation"
pol='s'
variable='wavelength'
n1 = 1.0
n2= 1.55
n3=1.0
theta=np.deg2rad(30)
wavelength=630
e1=114.96
d=e1/wavelength
"d est pour avoir une unique variable de e imposé et la longueur d'onde varie"

def calculate_R_coeff(): 
    theta_out1=angle(n1, n2, theta)
    theta_out2=angle(n2,n3,theta_out1)
        
    wave=np.linspace(200,1000,1000)
    d=e1/wave
        
    r12=refle_coeff(n1, n2, theta,theta_out1)
    r23=refle_coeff(n2, n3, theta_out1, theta_out2)
        
    t12=trans_coeff(n1,n2,theta, theta_out1)
    t23=trans_coeff(n2, n3,theta_out1, theta_out2)
        
    depha1=phaseshift(d, n2, e1,theta_out1)
        
    r123=(r12+r23*np.exp(2*depha1))/(1+r12*r23*np.exp(2*depha1))
    t123=(t12*t23*np.exp(depha1))/(1+r12*r23*np.exp(2*depha1))
        
    R=np.abs(r123)**2
    T=(np.abs(t123)**2)*(n3*theta_out2/n1*theta)
    A=1-R-T
    
    return {
        'R':R,
        'T':T,
        'A':A
    }



# wave=np.linspace(200,1000,1000)
# plt.figure(1)
# plt.plot(wave,calculate_R_coeff().get('R'),'r',label='Reflectance of structure')
# plt.plot(wave,calculate_R_coeff().get('T'),'b',label='Transmittance of structure')
# plt.plot(wave,calculate_R_coeff().get('A'),'g',label='Absorbance of structure')
# plt.xlabel('Wavelength[nm]')
# plt.ylabel('Efficiencies')
# plt.legend()
# plt.show()