#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 09:16:29 2021

@author: ayastan1930
"""

import pandas as pd
import numpy as np
import scipy.integrate as sp
import matplotlib.pyplot as plt
import math
from scipy.signal import argrelextrema
import csv

'''
x[0] - Sun x position
x[1] - Sun velocity in x direction
x[2] - Sun y position
x[3] - Sun velocity in y direction
x[4] - Jupiter x position
x[5] - Jupiter velocity in x direction
x[6] - Jupiter y position
x[7] - Jupiter velocity in y direction
x[8] - Earth x position
x[9] - Earth x velocity
x[10] - Earth y position
x[11] - Earth y velocity
'''


mass_b = 5.972e24
mass_star = 1*1.989*10**30
mass_sun = 2.5*1.989*10**30

b_star_ratio = mass_b/mass_star 
b_sun_ratio = mass_b/mass_sun

g_sun = 1.98279e-29*2.5*1.989e+30
g_star = 1.98279e-29*1*1.989e+30




def fod_planet(x, t):
    sun_star_vector = ((x[0]-x[4])**2 + (x[2]-x[6])**2)**1.5
    sun_b_vector = ((x[0]-x[8])**2 + (x[2]-x[10])**2)**1.5
    star_sun_vector = ((x[4] - x[0])**2 + (x[6] - x[2])**2)**1.5
    star_b_vector = ((x[4] - x[8])**2 + (x[6] - x[10])**2)**1.5
    b_sun_vector = ((x[8]-x[0])**2 + (x[10]-x[2])**2)**1.5
    b_star_vector = ((x[8] - x[4])**2 + (x[10] - x[6])**2)**1.5

  
    
    derivs = [x[1], -(4*np.pi**2+g_star)*(1/3.5)*((x[0]-x[4])/sun_star_vector) -
             
              (4*np.pi**2)*b_sun_ratio*((x[0]-x[8])/sun_b_vector),
              
              x[3], -(4*np.pi**2+g_star)*(1/3.5)*((x[2]-x[6])/sun_star_vector) -
              
              (4*np.pi**2)*b_sun_ratio*((x[2]-x[10])/sun_b_vector),
              
              x[5], -(4*np.pi**2+g_sun)*(2.5/3.5)*(x[4]-x[0])/star_sun_vector -
              
              (4*np.pi**2)*b_star_ratio*(x[4]-x[8])/star_b_vector,
              
              x[7], -(4*np.pi**2+g_sun)*(2.5/3.5)*(x[6]-x[2])/star_sun_vector -
              
              (4*np.pi**2)*b_star_ratio*(x[6]-x[10])/star_b_vector, 
              
              x[9], -4*(np.pi**2)*((x[8]-x[0])/b_sun_vector) -
              
              4*(np.pi**2)*((x[8]-x[4])/b_star_vector), 
              
              x[11], -(4*np.pi**2)*((x[10]-x[2])/b_sun_vector) -
              
              (4*np.pi**2)*((x[10]-x[6])/b_star_vector) 
              ]
    return derivs

peri_star = 1.5
peri_b = 5

period_star = ((peri_star)**3)**(1/2)
period_b = (peri_b**3)**(1/2)


sun_pos = (mass_star*peri_star)/(mass_sun+mass_star)
star_pos = peri_star - sun_pos
b_pos = peri_b - sun_pos

speed_sun = (2*np.pi*sun_pos)/(period_star)
speed_star = (2*np.pi*star_pos)/(period_star)

speed_sun_1 = 0.9
speed_star_1 = 3.94

b_peri_speed = (2*np.pi*b_pos)/(period_b)
b_peri_speed_1 = 4.2



t = np.linspace(0.0, 40, num=50000 ) 

ic = [-sun_pos, 0, 0, -speed_sun_1,
      star_pos, 0, 0, speed_star_1,
      b_pos, 0, 0, b_peri_speed_1
      ]



soln = sp.odeint(fod_planet, ic, t)

soln_ms = soln[:, 1]*4743.18

soln_AU_year = soln[:, 1]

df = pd.DataFrame({'rad_x_vel': soln_ms, 'years': t }, columns=['rad_x_vel', 'years'])


plto1= plt.figure(1)
plt.plot(soln[:, 0], soln[:,2],'red')
plt.plot(soln[:, 4], soln[:,6], 'orange')
plt.plot(soln[:, 8], soln[:,10], 'blue')
plt.grid('k')
plt.xlabel('X position (AU)')
plt.ylabel('Y position (AU)')
plt.title('Orbital Motion')


plot2 = plt.figure(2)
plt.plot(t, soln_ms,'red')
plt.grid('k')
plt.ylabel('X velocity (m/s)')
plt.xlabel('Time (years)')
plt.title('Radial Velocity Plot')


plt.show()


#rad_vel_habitable = soln_ms

#np.savetxt("data_4_binary.csv", df, delimiter=",")

df.to_csv('Sim3.csv', index = False)

