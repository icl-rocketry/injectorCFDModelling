#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 20 23:30:55 2021

@author: piotr
"""

from scipy.optimize import brentq
import math
import numpy as np
from matplotlib import pyplot as plt

def h_l(T):
    b_1 = -200.
    b_2 = 116.043
    b_3 = -917.225
    b_4 = 794.779
    b_5 = -589.587
    T_c = 309.57
    T = (T+273.15)/T_c
    return b_1 + b_2*(1-T)**(1/3) + b_3*(1-T)**(2/3) + b_4*(1-T) + b_5*(1-T)**(4/3)

def h_g(T):
    b_1 = -200.
    b_2 = 440.055
    b_3 = -459.701
    b_4 = 434.081
    b_5 = -485.338
    T_c = 309.57
    T = (T+273.15)/T_c
    return b_1 + b_2*(1-T)**(1/3) + b_3*(1-T)**(2/3) + b_4*(1-T) + b_5*(1-T)**(4/3)

def rho_l(T):
    rho_c = 452.
    b_2 = 1.72328
    b_3 = -0.83950
    b_4 =  0.51060
    b_5 = -0.10412
    T_c = 309.57
    T = (T+273.15)/T_c
    return rho_c * math.exp(b_2*(1-T)**(1/3) + b_3*(1-T)**(2/3) + b_4*(1-T) + b_5*(1-T)**(4/3))
    
def rho_g(T):
    rho_c = 452.
    b_2 = -1.00900
    b_3 = -6.28792
    b_4 = 7.50332
    b_5 = -7.90463
    b_6 = 0.629427
    T_c = 309.57
    T = (T+273.15)/T_c
    return rho_c * math.exp(b_2*(1/T-1)**(1/3) + b_3*(1/T-1)**(2/3) + b_4*(1/T-1) + b_5*(1/T-1)**(4/3) + b_6*(1/T-1)**(5/3))

def func(T):
    return (h_l(10)+((0.09/(rho_l(10)))**2)/2)-((0.09**2)/2*(1/rho_g(T))**2 + h_g(T) + h_l(T) - h_l(10))

def p(T):
    p_c = 7251e3
    b_2 = -6.71893
    b_3 = 1.35966
    b_4 =  -1.3779
    b_5 = -4.051
    T_c = 309.57
    T = (T)/T_c
    return p_c * math.exp(1./T*(b_2*(1.-T) + b_3*(1.-T)**(3./2.) + b_4*(1.-T)**(5./2.) + b_5*(1.-T)**(5.)))

def cp_r(T):
    b_1 = 132.632
    b_2 =  0.052187
    b_3 = -0.364923
    b_4 = -1.20233
    b_5 = 0.536141
    T_c = 309.57
    T = (T)/T_c
    return 1e3*b_1*(1. + b_2*(1.-T)**(-2./3.) + b_3*(1.-T)**(-1./3.) + b_4*(1.-T)**(1./3.) + b_5*(1.-T)**(2./3.))

def z(T,p):
    W  = 44.013
    NA = 6.0221417930e+23
    k = 1.38065e-23
    RR = NA*k*1e3
    P_c = 7251e3
    T_c = 309.57
    V_c = 97e-3
    omega = 0.16
    
    T_r = T/T_c
    a = 0.45724*(RR*T_c)**2./P_c
    b = 0.0778*RR*T_c/P_c
    kappa = 0.37464 + 1.54226*omega - 0.26992*omega**2.
    alpha = (1. + kappa*(1. - math.sqrt(T_r)))**2.
     
    A = a*alpha*p/(RR*T)**2.
    B = b*p/(RR*T)
    
    a2 = B -1.
    a1 = A -2.*B-3*B**2.
    a0 = -A*B + B**2. + b**3.
    
    Q = (3.*a1 - a2*a2)/9.
    Rl = (9.*a2*a1 -27.*a0 -2.*a2*a2*a2)/54.
    
    Q3 = Q*Q*Q
    D = Q3 + Rl*Rl
    
    root = -1
    
    if D <= 0:
        th = math.acos(Rl/math.sqrt(-Q3))
        qm = 2*math.sqrt(-Q)
        r1 = qm*math.cos(th/3.) - a2/3.
        r2 = qm*math.cos((th + 2.*math.pi)/3.) - a2/3.
        r3 = qm*math.cos((th + 4**math.pi)/3.) - a2/3.
        root = max([r1,r2,r3])
    else:
        D05 = math.sqrt(D)
        S = (Rl+D05)**(1./3.)
        Tl = 0
        if D05 > Rl:
            Tl = - (abs(Rl-D05))**(1./3.)
        else:
            Tl = (Rl-D05)**(1./3.)
        root = S+Tl-a2/3.
    return root

def cp_p(T,p):
    W  = 44.013
    NA = 6.0221417930e+23
    k = 1.38065e-23
    RR = NA*k*1e3
    P_c = 7251e3
    T_c = 309.57
    V_c = 97e-3
    omega = 0.16
    
    T_r = T/T_c
    a = 0.45724*(RR*T_c)**2./P_c
    b = 0.0778*RR*T_c/P_c
    kappa = 0.37464 + 1.54226*omega - 0.26992*omega**2.
    alpha = (1. + kappa*(1. - math.sqrt(T_r)))**2.
    
    A = a*alpha*p/(RR*T)**2.
    B = b*p/(RR*T)
    
    Z=z(T,p)
    
    ap = kappa*a*(kappa/T_c-(1.+kappa)/math.sqrt(T*T_c))
    app = kappa*a*(1. + kappa)/(2.*math.sqrt(T**3. * T_c))
    
    M = (Z**2. +2*Z*B - B**2.)/(Z-B)
    N = ap*B/(b*RR)
    
    root2 = math.sqrt(2.)
    return (app*(T/(2.*root2*b))*math.log((Z+(root2+1.)*B)/(Z-(root2-1.)*B)) + RR*(M-N)**2. / (M**2. -2.*A*(Z+B)) - RR)/W

def cp_i(T):
    b_1 = -0.169903
    b_2 =  0.099053
    b_3 =  1.20822
    b_4 = -0.248324
    T_c = 309.57
    T = (T)/T_c
    return 1e3*(b_1 + b_2*T**(-1./2.) + b_3*T**(1./2.) + b_4*T)
"""
T = np.linspace(182.33,300.57,500)
CP_R = np.linspace(182.33,300.57,500)
CP_P = np.linspace(182.33,300.57,500)

for i in range(0,500):
    CP_P[i] = cp_p(CP_P[i],p(CP_P[i]))+cp_i(CP_P[i])
    CP_R[i] = cp_r(CP_R[i])
    
plt.plot(T,CP_P)
plt.plot(T,CP_R)
#plt.plot(T,(CP_R-CP_P)/CP_P)
"""
atm = 1.01325e5

print(z(36.46+273.15,70*atm))