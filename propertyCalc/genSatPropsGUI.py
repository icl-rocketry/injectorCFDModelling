from CoolProp.CoolProp import PropsSI
from math import exp, pow
import PySimpleGUI as sg

def visc_sat_liq(T):
    Tc = 273.15 + 36.4
    b1 = 1.6089
    b2 = 2.0439
    b3 = 5.24
    b4 = 0.0293423
    theta = (Tc-b3)/(T-b3)
    eta = b4*exp(b1*pow(theta-1,1/3) + b2*pow(theta-1,4/3))*pow(10, -3)
    return eta

def visc_sat_vap(T): # takes T in degC
    Tc = 273.15 + 36.42
    b1 = 3.3281
    b2 = -1.18237
    b3 = -0.055155
    Tr = T/Tc
    eta = exp(b1 + b2*pow(1/Tr - 1, 1/3) + b3*pow(1/Tr - 1, 4/3))*pow(10, -6)
    return eta

def surface_tension(T): 
    Tc = 273.15 + 36.42
    b1 = 69.31
    b2 = 1.19346
    b3 = 0
    Tr = T/Tc
    sigma = b1*pow((1-Tr),b2)*(1+b3*(1-Tr))*pow(10, -3)
    return sigma

def update_disp(P,rhov,rhol,Vv,Vl,sigma,dHl,dHv):
    window["PRESSURE"].update(value="{:.4e}".format(P))
    window["TEMPERATURE"].update(value="{:.2f}".format(T))
    window["SIGMA"].update(value="{:.4e}".format(sigma))
    window["vapDensity"].update(value="{:.4e}".format(rhov))
    window["liqDensity"].update(value="{:.4e}".format(rhol))
    window["vapVisc"].update(value="{:.4e}".format(Vv))
    window["liqVisc"].update(value="{:.4e}".format(Vl))
    return

def calc_props(T,fluid):
    Pv = PropsSI('P','T',T,'Q',1,fluid)
    rhov = PropsSI('D','T',T,'Q',1,fluid)
    rhol = PropsSI('D','T',T,'Q',0,fluid)
    dHv = PropsSI("H","T",T,"Q",1,fluid)
    dHl = PropsSI("H","T",T,"Q",0,fluid)
    if fluid == 'NitrousOxide':
        Vv = visc_sat_vap(T)/rhov
        Vl = visc_sat_liq(T)/rhol
        sigma = surface_tension(T)
    else:
        Vv = PropsSI('V','T',T,'Q',1,fluid)/rhov
        Vl = PropsSI('V','T',T,'Q',0,fluid)/rhol
        sigma = PropsSI("surface_tension","T",T,"Q",0,fluid)
    update_disp(Pv,rhov,rhol,Vv,Vl,sigma,dHv,dHl)
    return

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

#sg.theme('Dark')

headings = ["State","Density (kg/m^3)","Kinematic Viscosity (m^2/s)"]
validFluids = ["NitrousOxide","CarbonDioxide","Water"]
layout = [ 
            [sg.Text("Choose fluid, or type valid fluid in")],
            [sg.InputCombo(validFluids,key="FLUID")],
            [sg.Text("Enter temperatures in °C")],
            [sg.Input(size=(12,1),key="TEMP_INPUT")],
            [sg.Button('Calculate')],
            [sg.Text("Temp: ",size=(12,1)),sg.Input("",key="TEMPERATURE", size=(10,1)),sg.Text("K")],
            [sg.Text("Pressure: ",size=(12,1)),sg.Input("",key="PRESSURE",size=(10,1)),sg.Text("Pa")],
            [sg.Text("Surface Tension: ",size=(12,1)),sg.Input("",key="SIGMA",size=(10,1)),sg.Text("N/m")],
            [sg.Text("State",size=(7,1)),sg.Text("Density (kg/m³)",size=(12,1)),sg.Text("Kinematic Viscosity (m²/s)",size=(20,1))],
            [sg.Text("vapour",size=(7,1)),sg.Input("",key="vapDensity",size=(10,1)),sg.Text("",size=(2,1)),sg.Input("",key="vapVisc",size=(10,1))
            ,sg.Text("",size=(9,1))],
            [sg.Text("liquid",size=(7,1)),sg.Input("",key="liqDensity",size=(10,1)),sg.Text("",size=(2,1)),sg.Input("",key="liqVisc",size=(10,1))
            ,sg.Text("",size=(9,1))]
         ]

window = sg.Window('saturation properties calculator', layout, grab_anywhere=True)

while 1:
    event, inputs = window.read()
    if event == sg.WIN_CLOSED:     # If user closed window with X or if user clicked "Exit" button then exit
        break
    try:
        T = 273.15 + float(inputs["TEMP_INPUT"])
        fluid = inputs["FLUID"]

        if event == 'Calculate':  
            Tcrit = PropsSI(fluid,"Tcrit")
            Ttriple = PropsSI(fluid,"Ttriple")
            if T > Tcrit:
                T = Tcrit
                window["TEMP_INPUT"].update(value="{:.2f}".format(T - 273.15))
            elif T < Ttriple:
                T = Ttriple
                window["TEMP_INPUT"].update(value="{:.2f}".format(T - 273.15))

            try:
               calc_props(T,fluid) 
            except ValueError:
                print("broken")
        
    except ValueError:
        print("invalid value of T")

window.close()
# pyinstaller -wF genSatPropsGUI.py
# references:
# https://webbook.nist.gov/cgi/fluid.cgi?TLow=+216.592&THigh=304.1282&TInc=1&Applet=on&Digits=5&ID=C124389&Action=Load&Type=SatP&TUnit=K&PUnit=bar&DUnit=kg%2Fm3&HUnit=kJ%2Fmol&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm&RefState=DEF
# http://edge.rit.edu/edge/P07106/public/Nox.pdf
# https://cdnsciencepub.com/doi/pdf/10.1139/v64-439#:~:text=The%20only%20previous%20work%20reported,and%2026.3%20dyn%20cm-l.