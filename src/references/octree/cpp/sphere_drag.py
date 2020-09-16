#!/bin/env python

# Reference:
# http://pages.mtu.edu/~fmorriso/DataCorrelationForSphereDrag2016.pdf

import math

vel = 30
Re=vel*1/0.2

cd = 24.0/Re*(1.0+0.125*math.pow(Re,0.72))
print("Lapple formula Drag coefficient Cd of sphere for Re=",Re," equals: ",cd)
#cd=24./Re*(1+0.1935*(math.pow(Re,0.6305)))
#print("Simple Drag coefficient Cd of sphere for Re=",Re," equals: ",cd)
cd=24./Re+(2.6*Re/5.)/(1+math.pow(Re/5.,1.52))+(0.411*math.pow(Re/(2.63*100000),-7.94))/(1+math.pow((Re)/(2.63*100000),-8.0))+(0.25*Re/1000000)/(1+Re/1000000)
print("Drag coefficient Cd of sphere for Re=",Re," equals: ",cd)

# https://www.sciencedirect.com/sdfe/pdf/download/eid/1-s2.0-0032591086800122/first-page-pdf
cd_2 = (math.pow(10,-0.7133+0.6305*math.log10(Re))+1)*24./Re
print("Drag (other formula) coefficient Cd of sphere for Re=",Re," equals: ",cd_2)

Ma=0.088
gam=1.4
TP=300.
TG=300.
ccd=24./(Re+Ma*math.sqrt(gam/2.)*(4.33+((3.65-1.53*TP/TG)/(1+0.353*TP/TG))*math.exp(-0.247*Re/(Ma*math.sqrt(gam/2.)))))+math.exp(-0.5*Ma/(math.sqrt(Re)))*((4.5+0.38*(0.03*Re+0.48*math.sqrt(Re)))/(1+0.03*Re+0.48*math.sqrt(Re))+0.1*Ma*Ma+0.2*math.pow(Ma,8))+0.6*Ma*math.sqrt(gam/2.)*(1-math.exp(-Ma/Re))
print("Compressible drag coefficient Cd of sphere for Re=",Re," Ma=",Ma," gam=",gam," TP=",TP, "TG=",TG," equals: ",ccd)



rho=1.225
radius=0.5
A=3.141592*radius*radius
drag = cd*0.5*rho*vel*vel*A
print("Drag for the given parameters equals: ",drag)

drag = ccd*0.5*rho*vel*vel*A
print("Compressible drag for the given parameters equals: ",drag)
