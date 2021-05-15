#!usr/bin/env python3

'''
  File name: evaluate_testing_policies_infection_days
  Author:    PFA
  Date:      Fri May 14, 2021

  Evaluate testing policies with respect to expected number of days infecting
'''

import numpy as np
import pandas as pd

from itertools import combinations, product

# Global Vars

T = 12

# Methods

def count_infecting(df, t):
  '''
  Count number of viral load curves in df that are infecting at day t
  '''
  
  c = str(t)
  aux = df.loc[df[c] >= 6]
  return aux.shape[0]

def filter_antigent_tested(df, t):
  '''
  Filter out viral load curves that are antigen-positive at day t
  '''
  
  c = str(t)
  return df.loc[df[c] < 5]
  
def filter_pcr_tested(df, t):
  '''
  Filter out viral load curves with a positive PCR result at day t
  '''
  
  c = str(t - 1)
  return df.loc[df[c] < 3]

def generate_budgeted_policies(num_antigen, num_pcr):
  '''
  Generate an array of possible testing policies with budgets for each type of test given by num_antigen and num_pcr
  '''
  
  antigen_days = combinations(range(1,  T + 1), num_antigen)
  pcr_days = combinations(range(2, T + 1), num_pcr)
  
  return product(antigen_days, pcr_days)

def evaluate_policy(df, policy):
  '''
  Evaluate expected cost of a policy.
  A policy consists in a tuple, where the first entry are the antigen result test days, and the second are the PCR result days.
  '''
  
  aux = df.copy()
  num_mc = df.shape[0]
  
  D_a = policy[0]
  D_p = policy[1]
  
  prob_inf = 0
  
  for t in range(1, T + 1):
    if t - 1 in D_a: 
      aux = filter_antigent_tested(aux, t - 1)
    if t - 1 in D_p:
      aux = filter_pcr_tested(aux, t - 1)
      
    num = count_infecting(aux, t)
    den = aux.shape[0]
    if den != 0:
      prob_inf += num / den
    
  return prob_inf

def main():
  df = pd.read_csv('viral_load_mc.csv')
  
  policies = generate_budgeted_policies(2, 1)
  policy_cost = []
  for p in policies:
    prob_inf = evaluate_policy(df, p)
    aux_df = pd.DataFrame({'antigen_days': [p[0]],
                           'pcr_days': [p[1]], #[[d - 1 for d in p[1]]],
                           'expected_infecting_days': prob_inf})
    policy_cost.append(aux_df)
  
  pd.concat(policy_cost).to_csv('results2.csv')
    
main()
  