# Project Name: Prey and Predator Population Model
# Date: 06/20/2019
# Compiler type: PyCharm
# File Type: Python 3
from numpy import *
import pylab as p
from scipy.integrate import odeint
import numpy as np

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


#ignores division by zero error in calculation effective population


np.seterr(divide='ignore')


'''# Purpose: Initialization of Variables for Lokta Volterra Equation
prey = 1.1
predator = 0.4
# Prey and Predator Species names
preySpecies = "Rabbits"
predatorSpecies = "Wolves"
# Time period, how frequently the addition of populations differ
timePeriod = 0.1
# preyGrowthRate, how fast the prey pop grows
carryingCapacity= 5
preyGR = 1.3
attackRate = 0.5
# Ratio of death rate for predator populations
predatorDeathRate = 0.7
# Rate at how efficiently predator reproduce based off the prey
predatorEfficiencyRate = 1.6'''

# ODE DE equations

# the natural growth rate of species1
a = 1
# death rate of rabbits due to species2
b = 0.1
# natural death rate of predator when there are no prey
c = 1.5
# the facor describing how many predator eat the prey
d = 0.75


# Calculates growth rate
def dX_dt(X, t=0):
   # Returns the growth rate of predator and prey populations
   return array([a * X[0] + b * X[0] * X[1], c * X[1] + d * b * X[0] * X[1]])


# Population Equilibirum: equilibirum occurs when the growth rate equals 0, this will give two fixed points
X_f0 = array([0.0, 0.0])
X_f1 = array([c / (d * b), a / b])
all(dX_dt(X_f0) == zeros(2)) and all(dX_dt(X_f1) == zeros(2))  # => True


# Stability of the fixed points
# dX_dt = A_f*X where A is the Jacobian matrix evaluated at the corresponding point.
# We have to define the Jacobian matrix:
def d2X_dt2(X, t=0):
   # Return the jacobian matrix evaluated in X
   return array([[a - b * X[1], -b * X[0]], [b * d * X[1], -c + b * d * X[0]]])

A_f1 = d2X_dt2(X_f1)

# Eigenvalue are +/- square root (c*a).j:
lambda1, lambda2 = linalg.eigvals(A_f1)

# They are imaginary numbers so the predator and prey populations are periodic and the periods are given by
T_f1 = 2 * pi / abs(lambda1)

# integrate the ODE with scipy.integrate

from scipy import integrate

# time t = linspace(start,stop,num)
# Returns num(the third parameter) evenly spaced samples, calculated over the interval [start, stop].
t = linspace(0, 1000, 1000)

# initial conditions, 10 predator and 5 prey
X0 = array([10, 5])

X, infodict = integrate.odeint(dX_dt, X0, t, full_output=True)
infodict['message']

#Linear algebra transposition of 2d-array(matrix)
rabbits, wolves = X.T


# Creating a list data structure for rabbits and wolves
ne_r = list(rabbits)
ne_w = list(wolves)

#function to calculate effective population of rabbit
def effectivePopRabbits(X,t):
   summation = 0

   for i in range(0, t.size):
       summation += 1 / rabbits[i]
       X[i] = float((1 / (t[i] * float(summation))))
    #returns ne_r
   return X


#function to calculate effective population of wolves
def effectivePopWolves(X,t):
   summation = 0

   for i in range(0, t.size):
       # print("this is wolves[i]: " + str(wolves[i]))
       summation += 1 / wolves[i]
       # print("this is summation: " + str(summation))

       X[i] = float((1 / (t[i] * float(summation))))
    #returns ne_w
   return X

#Fills ne_r and new_w with effective populations for rabbits and wolves
ne_r = effectivePopRabbits(ne_r,t)
ne_w = effectivePopWolves(ne_w,t)



# # Population Scaled Mutation Rates ( 4 NEu) u= specific mutation rates, ne = effectie pop size
#
# # Pop scaled mutation rates in predator
# mutationRate_r=[]
# mutationRate_w=[]
# u=0
#
# for i in range(0,ne_r.size):
#     mutationRate_r[i]= (4*(ne_r[i]*u))
#     mutationRate_w[i] = (4*(ne_r[i]*u))
#
# # Pop scaled mutation rates in prey
# # predatorMutationRatesArray=(4 * (preyEffectiveArray[i] * u)




# Graph of population census of predator and prey
f1 = p.figure()
p.plot(t, rabbits, 'r-', label='Rabbits')
p.plot(t, wolves, 'b-', label='wolves')
p.grid()
p.legend(loc='best')
p.xlabel('time')
p.ylabel('population')
p.title('Wolves and rabbit populations')
f1.savefig('rabbits_and_wolves_1.png')
plt.show()

# Graph of
f2 = p.figure(2)
# p.plot(t,rabbits, 'r-', label='Rabbits' )
# p.plot(t, wolves, 'b-', label='wolves')
p.plot(t, ne_r, 'r--', label='Rabbits_m')
p.plot(t, ne_w, 'b--', label='wolves_m')
p.grid()
p.legend(loc='best')
p.xlabel('time')
p.ylabel('population')
p.title('Wolves and rabbit populations $(1/N_e)$')
f2.savefig('rabbits_and_wolves_1.png')
plt.show()


#Graph of mutation rates

# f3=p.figure(3)
# p.plot(t,rabbits, 'r-', label='Rabbits' )
# p.plot(t, wolves, 'b-', label='wolves')
# p.plot(t,r_hm, 'r--', label='Rabbits_m' )
# p.plot(t, w_hm, 'b--', label='wolves_m')
# p.grid()
# p.legend(loc='best')
# p.xlabel('time')
# p.ylabel('population')
# p.title('Wolves and rabbit populations $(1/N_e)$')
# f3.savefig('rabbits_and_wolves_1.png')
# plt.show()