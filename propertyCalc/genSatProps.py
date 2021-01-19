from CoolProp.CoolProp import PropsSI
from math import exp, pow

def display_menu(P,rhov,rhol,Vv,Vl,sigma):
    print("Temp: {:.2f} K".format(T))
    print("Pressure: {:.4e} Pa".format(P))
    print("Surface Tension: {:.4e}".format(sigma))
    print("-------------------------------------------------------------------")
    print("State\t\tDensity (kg/m^3)\tKinematic Viscosity (m^2/s)")
    print("-------------------------------------------------------------------")
    print("vapour\t\t{:.4e}\t\t{:.4e}".format(rhov,Vv/rhov))
    print("liquid\t\t{:.4e}\t\t{:.4e}".format(rhol,Vl/rhol))
    print("-------------------------------------------------------------------")

def visc_sat_liq(T):
    Tc = 273.15 + 36.4
    b1 = 1.6089
    b2 = 2.0439
    b3 = 5.24
    b4 = 0.0293423
    theta = (Tc-b3)/(T-b3)
    eta = b4*exp(b1*pow(theta-1,1/3) + b2*pow(theta-1,4/3))
    return eta

def visc_sat_vap(T):
    Tc = 273.15 + 36.4
    b1 = 3.3281
    b2 = -1.18237
    b3 = -0.055155
    Tr = T/Tc
    eta = exp(b1 + b2*pow(1/Tr - 1, 1/3) + b3*pow(1/Tr - 1, 4/3))
    return eta

def surface_tension(T):
    Tc = 273.15 + 36.4
    b1 = 69.31
    b2 = 1.19346
    b3 = 0
    Tr = T/Tc
    sigma = b1*pow((1-Tr),b2)*(1+b3*(1-Tr))
    sigma_SI = sigma/1e5 *1e2
    return sigma_SI
# for using coolprop:
# T = temperature
# P = pressure
# H = enthalpy
# C = specific heat
# D = density
# V = viscosity
# Q = quality
# U = internal energy
# S = entropy
# L = thermal conductivity
print("enter temps in °C")
while 1:
    fluid = 'NitrousOxide' #input("fluid: ")
    T = 273.15 + float(input("temp: "))   
    try:
        Pv = PropsSI('P','T',T,'Q',1,fluid)
        rhov = PropsSI('D','T',T,'Q',1,fluid)
        rhol = PropsSI('D','T',T,'Q',0,fluid)
        #if fluid == 'NitrousOxide':
        Vv = visc_sat_vap(T)
        Vl = visc_sat_liq(T)
        sigma = surface_tension(T)
        #else:
        #    Vv = PropsSI('V','T',T,'Q',1,fluid)
        #    Vl = PropsSI('V','T',T,'Q',0,fluid)
        display_menu(Pv,rhov,rhol,Vv,Vl,sigma)
    except ValueError:
        print("outside of saturation curve")
# references:
# https://webbook.nist.gov/cgi/fluid.cgi?TLow=+216.592&THigh=304.1282&TInc=1&Applet=on&Digits=5&ID=C124389&Action=Load&Type=SatP&TUnit=K&PUnit=bar&DUnit=kg%2Fm3&HUnit=kJ%2Fmol&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm&RefState=DEF
# http://edge.rit.edu/edge/P07106/public/Nox.pdf
# https://cdnsciencepub.com/doi/pdf/10.1139/v64-439#:~:text=The%20only%20previous%20work%20reported,and%2026.3%20dyn%20cm-l.