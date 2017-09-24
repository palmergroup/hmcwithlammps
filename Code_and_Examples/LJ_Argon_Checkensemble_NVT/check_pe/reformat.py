#!/usr/bin/python
import string
import os
import sys
import numpy as numpy

nskp = 0
neq = 20000
nfooter = 0

t1=numpy.genfromtxt('../132K/thermo.dat', skip_header=nskp+neq, skip_footer=nfooter)
npts = t1[::1,0].size
t1out= open('t1.dat', 'w')
for i in range(npts):
        t1out.write("%f \n" % (4.184*t1[i,1]))
t1out.close()

t2=numpy.genfromtxt('../145K/thermo.dat', skip_header=nskp+neq, skip_footer=nfooter)
npts = t2[::1,0].size
t2out= open('t2.dat', 'w')
for i in range(npts):
        t2out.write("%f \n" % (4.184*t2[i,1]))
t2out.close()


print numpy.average(t1[:,1]),numpy.std(t1[:,1]), numpy.average(t2[:,1]),numpy.std(t2[:,1])
