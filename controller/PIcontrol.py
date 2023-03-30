# Explicit feedback for climate modeling
# Copyright (C) 2020  Ben Kravitz
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# Sample control parameters file
#
# Written by Ben Kravitz (bkravitz@iu.edu or ben.kravitz.work@gmail.com)
# Last updated 11 July 2019
#
# This script provides information about the feedback algorithm.  All of this
# is user-customizable.  The other parts of the script will give you outvals,
# lats, lons, and times.  The output of this script should be a list called
# nlvals, which consists of pairs.  The first item in the pair is the name
# of the namelist value.  The second item is the value itself as a string.
#
# This script is written in native python format.  Be careful with brackets [],
# white space, and making sure everything that needs to be a string is actually
# a string by putting it in quotes ''.  All lines beginning with # are comments.

#### USER-SPECIFIED CONTROL PARAMETERS ####
#refvals=[288.13,0.76,-5.98] reference values for CESM2(WACCM6) (Tilmes et al. 2020) 2020-2039
refvals=[288.64,0.8767,-5.89] # updated to be average over years 2010-2029
new_refvals=[288.14,0.8497,-5.90] # new target values: averages during the years 2000-2019 (in which T0 is ~0.5 degrees lower than the old targets)
#new_refvals=[287.64,0.7348,-5.97] # new target values: averages during the years 1984-2003 (in which T0 is ~1.0 degrees lower than the old targets)

#refvals=[288.21,0.594,-6.006] # new version of the model (GLENS values)
kivals=[0.0183,0.0753,0.3120] #[0.0183,0.0753,0.3120]
kpvals=[0.0183,0.0753,0.3120] #[0.0183,0.0753,0.3120]
# kivals=[0.,0.,0.]
# kpvals=[0.,0.,0.]
firstyear=2035
baseyear=2030#2030
#Change to Zero
x_ramp = 0#5.0 # defines a range of years over which the feedback is ramped up


#### USER SPECIFIED CALCULATIONS ####
logfilename='ControlLog_'+runnames+'.txt'

logheader=['Timestamp','dT0','sum(dT0)','dT1','sum(dT1)','dT2','sum(dT2)','L0','L1N','L1S','L2','30S(Tg)','15S(Tg)','15N(Tg)','30N(Tg)']

firsttime=0
if os.path.exists(maindir+logfilename)==False:
    # print("This is the first file")
    firsttime=1
else:
    # print("this continues a file")
    loglines=readlog(maindir+logfilename)

import numpy as np
w=makeweights(lats,lons)

outvals = numpy.array(outvals)

T0=numpy.mean(gmean(outvals[0],w)) #288.64 #variable#numpy.mean(gmean(outvals[0],w))
T1=numpy.mean(l1mean(outvals[0],w,lats)) #0.8497 #numpy.mean(l1mean(outvals[0],w,lats))
T2=numpy.mean(l2mean(outvals[0],w,lats)) #-5.89 #numpy.mean(l2mean(outvals[0],w,lats))

# T0 = 290.38532417121905
# T1 =  0.8008744066435923
# T2 = -5.7479355939554795
# print(T0, T1, T2)

de=numpy.array([T0-refvals[0],T1-refvals[1],T2-refvals[2]]) # error terms

timestamp=firstyear
sumde=de
sumdt2=de[2]

# if firsttime==1:
#     timestamp=firstyear
#     sumde=de
#     sumdt2=de[2]
# else:
#     timestamp=int(loglines[-1][0])+1
#     sumdt0=float(loglines[-1][2])+(T0-refvals[0])
#     sumdt1=float(loglines[-1][4])+(T1-refvals[1])
#     sumdt2=float(loglines[-1][6])+(T2-refvals[2])
#     sumde=numpy.array([sumdt0,sumdt1,sumdt2])

# dt should be 1
dt=1 #timestamp-baseyear
dt2=timestamp-firstyear
# feedforward
#l0hat=0.011*dt
#l1hat=-0.005*dt
#l2hat=0.006*dt
# updated based on feedback simulation
l0hat=0.00347*dt
l1hat=-0.0*dt
l2hat=0.0*dt

############
####Changed to zero
ramp_up = 1.0
# if (dt2<x_ramp):
#     ramp_up = dt2 / x_ramp

# feedback EQUATION 5 Kravits 2017
l0kp1=(kpvals[0]*de[0]+kivals[0]*sumde[0])*ramp_up
l1kp1=(kpvals[1]*de[1]+kivals[1]*sumde[1]-0.5*l0kp1)*ramp_up
l2kp1=(kpvals[2]*de[2]+kivals[2]*sumde[2]-l0kp1)*ramp_up
############3

# all of the feeds
l0step4=l0kp1+l0hat
l1step4=l1kp1+l1hat
l2step4=l2kp1+l2hat

l0=max(l0step4,0)
l1n=min(max(l1step4,0),l0)
l1s=min(max(-l1step4,0),l0)
l2=min(max(l2step4,0),l0-l1s-l1n)
ell=numpy.array([[l0],[l1n],[l1s],[l2]])


# preventing integrator wind-up
if (l2==(l0-l1s-l1n)):
    sumdt2=sumdt2-(T2-refvals[2])
    sumde[2]=sumdt2

M=numpy.array([[0,30,30,0],[0,0,45,20],[20,45,0,0],[40,0,0,40]])
F=numpy.array([[1,1,1,1],[0,1,0,0],[0,0,1,0],[0,0,0,1]])

q=numpy.dot(numpy.dot(numpy.transpose(M),numpy.linalg.inv(F)),ell)

# sensitivity[:,la,lo] = q[:,0]

for k in range(len(q)):
    q[k]=max(q[k],0)

newline=[str(data_type),str(timestamp),str(de[0]),str(sumde[0]),str(de[1]),str(sumde[1]),str(de[2]),str(sumde[2]),str(l0),str(l1n),str(l1s),str(l2),str(q[0])[1:-1],str(q[1])[1:-1],str(q[2])[1:-1],str(q[3])[1:-1]]
if firsttime==1:
    linestowrite=[logheader,newline]
else:
    linestowrite=[]
    for k in range(len(loglines)):
        linestowrite.append(loglines[k])
    linestowrite.append(newline)

writelog(maindir+'/'+logfilename,linestowrite)

