#!usr/bin/env python3

'''
  File name: generate_montecarlo_viral_load.py
  Author:    PFA
  Date:      Fri May 14, 2021

  Using MonteCarlo integration, compute the probabilities needed in the DP for the optimal testing policy of suspect Covid-19 cases.
'''

import numpy as np
import pandas as pd

# Global Vars

T = 14  # Days to simulate
np.random.seed(2)
N = 1e5 # Number of simulations

# Methods

def generate_viral_load_curve(t_0, t_p, t_f, v_p):
  '''
  Generate viral load curve V_t given the parameters
  '''

  m1 = (v_p - 3) / (t_p - t_0)
  m2 = (6 - v_p) / (t_f - t_p)
  b1 = 3 - m1 * t_0
  b2 = v_p - m2 * t_p
  
  l1 = [m1 * t + b1 for t in range(T + 1) if t <= t_p]
  l2 = [m2 * t + b2 for t in range(T + 1) if t > t_p]
  
  return l1 + l2
  
def generate_random_variables():
  '''
  Generate a set of random variables (one draw from the Monte Carlo simulation)
  '''
  
  t_0 = np.random.uniform(low=2.5, high=3.5)
  t_p = t_0 + min(3, 0.5 + np.random.gamma(shape=1.5))
  v_p = np.random.uniform(low=7, high=11)
  symptoms_flag = np.random.randint(low=0, high=2)
  t_s = t_p + np.random.uniform(low=0, high=3) * symptoms_flag
  t_f = t_s + np.random.uniform(low=4, high=9)
  
  return t_0, t_p, v_p, t_s, t_f, symptoms_flag

# Main

def main():
  mc_matrix = []
  for it in range(int(N)):
    t_0, t_p, v_p, t_s, t_f, symptoms_flag = generate_random_variables()
    row = generate_viral_load_curve(t_0, t_p, t_f, v_p) + [t_s, symptoms_flag]
    mc_matrix.append(row)
  
  results_df = pd.DataFrame(mc_matrix)
  results_df.to_csv('viral_load_mc.csv')
  
main()