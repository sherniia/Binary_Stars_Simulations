#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  6 20:38:56 2021

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

#mass_b = 5.972e24
mass_b = 5.972e24
#mass_c = 1.981*1.89813*10**27
#mass_d = 4.132*1.89813*10**27
mass_star = 2*1.989*10**30
mass_sun = 3*1.989*10**30

b_star_ratio = mass_b/mass_star 
b_sun_ratio = mass_b/mass_sun

g_sun = 1.98279e-29*3*1.989e+30
g_star = 1.98279e-29*2*1.989e+30



alpha = 1.1*10**-8


def fod_planet(x, t):
    sun_star_vector = ((x[0]-x[4])**2 + (x[2]-x[6])**2)**1.5
    sun_b_vector = ((x[0]-x[8])**2 + (x[2]-x[10])**2)**1.5
    #sun_b_vector = ((x[0]-x[12])**2 + (x[2]-x[14])**2)**1.5
    star_sun_vector = ((x[4] - x[0])**2 + (x[6] - x[2])**2)**1.5
    star_b_vector = ((x[4] - x[8])**2 + (x[6] - x[10])**2)**1.5
    #b_d_vector = ((x[4] - x[12])**2 + (x[6] - x[14])**2)**1.5
    b_sun_vector = ((x[8]-x[0])**2 + (x[10]-x[2])**2)**1.5
    b_star_vector = ((x[8] - x[4])**2 + (x[10] - x[6])**2)**1.5
    #c_b_vector = ((x[8] - x[12])**2 + (x[10] - x[14])**2)**1.5
    #d_sun_vector = ((x[12]-x[0])**2 + (x[14]-x[2])**2)**1.5
    #d_b_vector = ((x[12] - x[4])**2 + (x[14] - x[6])**2)**1.5
    #d_c_vector = ((x[12] - x[8])**2 + (x[14] - x[10])**2)**1.5
    
    #rel = 1+(alpha/(x[0]**2 + x[2]**2))
    
    
    derivs = [x[1], -(4*np.pi**2+g_star)*(2/5)*((x[0]-x[4])/sun_star_vector) -
             
              (4*np.pi**2)*b_sun_ratio*((x[0]-x[8])/sun_b_vector),# -
              
              #4*(np.pi**2)*d_sun_ratio*(-(x[0]-x[12])/sun_b_vector),
              
              x[3], -(4*np.pi**2+g_star)*(2/5)*((x[2]-x[6])/sun_star_vector) -
              
              (4*np.pi**2)*b_sun_ratio*((x[2]-x[10])/sun_b_vector),# -
              
              #4*(np.pi**2)*d_sun_ratio*(-(x[2]-x[14])/sun_d_vector),
              
              x[5], -(4*np.pi**2+g_sun)*(3/5)*(x[4]-x[0])/star_sun_vector -
              
              (4*np.pi**2)*b_star_ratio*(x[4]-x[8])/star_b_vector,# -
              
              #4*(np.pi**2)*d_sun_ratio*(x[4]-x[12])/b_d_vector,
              
              x[7], -(4*np.pi**2+g_sun)*(3/5)*(x[6]-x[2])/star_sun_vector -
              
              (4*np.pi**2)*b_star_ratio*(x[6]-x[10])/star_b_vector, # - does one needs to include the mass ratio of the star with the planet?
              
             # -4*(np.pi**2)*d_sun_ratio*(x[6]-x[14])/b_d_vector,
              
              x[9], -4*(np.pi**2)*((x[8]-x[0])/b_sun_vector) -
              
              4*(np.pi**2)*((x[8]-x[4])/b_star_vector), #-
              
              #4*(np.pi**2)*d_sun_ratio*(-(x[8]-x[12])/c_d_vector),
              
              x[11], -(4*np.pi**2)*((x[10]-x[2])/b_sun_vector) -
              
              (4*np.pi**2)*((x[10]-x[6])/b_star_vector) #-
              
              #4*(np.pi**2)*d_sun_ratio*(-(x[10]-x[14])/c_d_vector),
              
              #x[13], 4*(np.pi**2)*(-(x[12]-x[0])/d_sun_vector) -  
              
             # 4*(np.pi**2)*b_sun_ratio*(-(x[12]-x[4])/d_b_vector)-
              
              #4*(np.pi**2)*c_sun_ratio*(-(x[12]-x[8])/d_c_vector),
              
             # x[15], 4*(np.pi**2)*(-(x[14]-x[2])/d_sun_vector) - 
              
              #4*(np.pi**2)*b_sun_ratio*(-(x[14]-x[6])/d_b_vector)-
              
              #4*(np.pi**2)*c_sun_ratio*(-(x[14]-x[10])/d_c_vector)
              ]
    return derivs

peri_star = 3
peri_b = 9.5

period_star = ((peri_star)**3)**(1/2)
period_star_1 = (((peri_star**3)*4*np.pi**2)/(1.98279e-29*(mass_sun+mass_star)))**(1/2)

period_b_1 = (((peri_b**3)*4*np.pi**2)/(1.98279e-29*(mass_star+mass_b+mass_sun)))**(1/2)
period_b = (peri_b**3)**(1/2)


sun_pos = (mass_star*peri_star)/(mass_sun+mass_star)
star_pos = peri_star - sun_pos
#star_pos = (mass_star*peri_star)/(mass_sun+mass_star)
#sun_pos = (mass_star*peri_star)/(2*mass_star+mass_b)
#star_pos = peri_star - sun_pos
b_pos = peri_b - sun_pos

#star_peri_speed = (np.pi*2*star_pos)/period_star
speed_sun = (2*np.pi*sun_pos)/(period_star)
speed_star = (2*np.pi*star_pos)/(period_star)

#speed_sun_1 = 2.2


speed_sun_1 = 1.3
speed_star_1 = 2.6

b_peri_speed = (2*np.pi*b_pos)/(period_b)
b_peri_speed_1 = 3 # the calculated speed is 1.78



t = np.linspace(0.0, 70, num=500000 ) 
t_1 = t*(3.154*10**7)

ic = [-sun_pos, 0, 0, -speed_sun_1,
      star_pos, 0, 0, speed_star_1,
      b_pos, 0, 0, b_peri_speed_1
      #pos_c, 0, 0, c_peri_speed
      #pos_d, 0, 0, d_peri_speed
      ]



soln = sp.odeint(fod_planet, ic, t)

soln_ms = soln[:, 1]*4743.18
soln_p = soln[:, 9]*4743.18
sun_x = soln[:, 0]
sun_y = soln[:, 2]
star_x = soln[:, 4]
star_y= soln[:, 6]
planet_x = soln[:, 8]
planet_y = soln[:,10]

soln_AU_year = soln[:, 1]

df = pd.DataFrame({'rad_x_vel': soln_ms, 'years': t }, columns=['rad_x_vel', 'years'])



df1 = pd.DataFrame({'x_pos_1': sun_x, 'y_pos_1': sun_y , 'x_pos_2': star_x, 'y_pos_2':
                    star_y, 'planet_x':planet_x, 'planet_y':planet_y}, columns=['x_pos_1', 'y_pos_1', 'x_pos_2', 'y_pos_2', 'planet_x', 'planet_y'])


plto1= plt.figure(1)
plt.plot(soln[:, 0], soln[:,2],'red')
plt.plot(soln[:, 4], soln[:,6], 'orange')
#plt.plot(soln[:, 8], soln[:,10], 'blue')
#plt.plot(soln[:, 12], soln[:,14], 'cyan')
plt.grid('k')
plt.xlabel('X position (AU)')
plt.ylabel('Y position (AU)')


'''
plot2 = plt.figure(2)
plt.plot(t, soln_ms,'red')
plt.grid('k')
plt.ylabel('X velocity (m/s)')
plt.xlabel('Time (years)')

'''

plot2 = plt.figure(2)
plt.plot(t, soln_ms,'red')
plt.grid('k')
plt.ylabel('X velocity (m/s)')
plt.xlabel('Time (years)')

'''

plot2 = plt.figure(3)
plt.plot(t, soln_ms,'red')
plt.grid('k')
plt.ylabel('X velocity in m/s')
plt.xlabel('time in seconds')


plo3 = plt.figure(3)
plt.plot(t_1,soln_ms*np.cos(0.5235), '-')
plt.grid('k')
plt.ylabel('X velocity in m/s')
plt.xlabel('time in seconds')
'''
plt.show()


#rad_vel_habitable = soln_ms

#np.savetxt("data_4_binary.csv", df, delimiter=",")

df.to_csv('Sim4.csv', index = False)


