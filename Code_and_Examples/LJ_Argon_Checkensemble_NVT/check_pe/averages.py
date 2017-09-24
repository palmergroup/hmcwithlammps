#!/usr/bin/python
import string
import os
import sys
import numpy as np
import timeseries




dat1 = np.loadtxt("t1.dat", delimiter=None)
dat2 = np.loadtxt("t2.dat", delimiter=None)

g1 = timeseries.statisticalInefficiency(dat1,fast=True)
g2 = timeseries.statisticalInefficiency(dat2,fast=True)

nsamp1 = dat1.size/g1
nsamp2 = dat1.size/g2

avg1 = np.mean(dat1)
avg2 = np.mean(dat2)

sig1 = np.power(np.var(dat1)/nsamp1,0.5)
sig2 = np.power(np.var(dat2)/nsamp2,0.5)

print avg1, 3.*sig1
print avg2, 3.*sig2
