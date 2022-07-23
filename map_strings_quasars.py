# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 23:33:32 2022

@author: olgso
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
#coordinates of quasar pairs from google doc
pair_RA=[214.041, 214.041, 161.78, 161.78, 185.576,185.576,161.772,164.452,164.452,170.193,170.193,218.181,216.585,291.452]
pair_DEC=[-30.6038, -30.6034, -31.6942,-31.6939, -36.1837, -36.1835,-32.7616, -31.7225,-31.7219,-32.3082,-32.3075, -30.1963,-30.1976,-41.1854]
#working with icrf .txt
#4536 sourses
icrf_RA=[0]*4536
icrf_DEC=[0]*4536
icrf_namez=['0']*4536
f=open('icrf.txt','r')
for i in range(4536):
    f.read(5)
    icrf_namez[i]=f.read(16)
    f.read(19)
    ra_1=f.read(2)
    ra_1=float(ra_1)
    f.read(1)
    ra_2=f.read(2)
    ra_2=float(ra_2)
    f.read(1)
    ra_3=f.read(11)
    ra_3=float(ra_3)
    icrf_RA[i]=15*(ra_1+ra_2/60+ra_3/3600)
    f.read(4)
    minus=f.read(1)
    dec_1=f.read(2)
    dec_1=float(dec_1)
    f.read(1)
    dec_2=f.read(2)
    dec_2=float(dec_2)
    f.read(1)
    dec_3=f.read(10)
    dec_3=float(dec_3)
    icrf_DEC[i]=dec_1+dec_2/60+dec_3/3600
    if minus=='-':
        icrf_DEC[i]=icrf_DEC[i]*(-1)
    f.readline()
f.close()
'''f=open('icrf_better.txt','w')
for i in range(4536):
    f.write(icrf_namez[i])
    f.write(' ')
    f.write(str(icrf_RA[i]))
    f.write(' ')
    f.write(str(icrf_DEC[i]))
    f.write('\n')
f.close()'''
icrf_numberz=[3590,3632,3704,3759,3894,3999,4020,4078]
icrf_string_RA=[0]*len(icrf_numberz)
icrf_string_DEC=[0]*len(icrf_numberz)
for i in range(len(icrf_numberz)):
    icrf_string_RA[i]=icrf_RA[icrf_numberz[i]]
    icrf_string_DEC[i]=icrf_DEC[icrf_numberz[i]]
#converting RA to -180 to 180 format
for i in range(4536):
    if icrf_RA[i]>180:
        icrf_RA[i]=(360-icrf_RA[i])*(-1)
for i in range(len(icrf_numberz)):
    if icrf_string_RA[i]>180:
        icrf_string_RA[i]=(360-icrf_string_RA[i])*(-1)
for i in range(len(pair_RA)):
    if pair_RA[i]>180:
        pair_RA[i]=(360-pair_RA[i])*(-1)

#converting to radians
for i in range(len(pair_RA)):
    pair_RA[i]=pair_RA[i]*(np.pi/180)
for i in range(len(pair_DEC)):
    pair_DEC[i]=pair_DEC[i]*(np.pi/180)
for i in range(len(icrf_RA)):
    icrf_RA[i]=icrf_RA[i]*(np.pi/180)
for i in range(len(icrf_DEC)):
    icrf_DEC[i]=icrf_DEC[i]*(np.pi/180)
for i in range(len(icrf_string_RA)):
    icrf_string_RA[i]=icrf_string_RA[i]*(np.pi/180)
for i in range(len(icrf_string_DEC)):
    icrf_string_DEC[i]=icrf_string_DEC[i]*(np.pi/180)
np_icrf_RA=np.array(icrf_RA)
np_icrf_DEC=np.array(icrf_DEC)
np_icrf_string_RA=np.array(icrf_string_RA)
np_icrf_string_DEC=np.array(icrf_string_DEC)
np_pair_RA=np.array(pair_RA)
np_pair_DEC=np.array(pair_DEC)


plt.figure()
plt.subplot(projection="mollweide")
plt.title("")
plt.grid(True)
marker = '.'
marker2='+'
ms, mew = 1, 0.5
ms2,mew2 =10, 1.
ms3,mew3=10,1.
plt.plot(np_pair_RA, np_pair_DEC,'ro', color='red', marker=marker2, ms=ms2, mew=mew2,label='pairs')
plt.plot(np_icrf_RA, np_icrf_DEC,'ro', color='blue', marker=marker, ms=ms, mew=mew,label='ICRF quasars')
plt.plot(np_icrf_string_RA, np_icrf_string_DEC,'ro', color='green', marker=marker2, ms=ms3, mew=mew3, label='ICRF quasars in 1 deg vivinity of a string')
plt.legend(loc='lower right',bbox_to_anchor=(0.5,-0.4, 0.5, 0.5))
plt.show()
