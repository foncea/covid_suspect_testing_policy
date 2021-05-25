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

import time

# Global Vars

T = 7

# Methods

def count_infecting(df, t):
  '''
  Count number of viral load curves in df that are infecting at day t
  '''
  
  c = str(t)
  aux = df.loc[df[c] >= 6]
  return aux.shape[0]

def count_infecting_corr(df, t, s):
  '''
  Count number of viral load curves in df that are infecting at day t AND also at day s.
  Method necessary to compute the variance of the estimator (covariance across time).
  '''
  
  c_t = str(t)
  c_s = str(s)
  aux = df[(df[c_t] >= 6) & (df[c_s] >= 6)]
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

def filter_symptomatic(df, t):
  '''
  Filter out symptomatic individuals
  '''
  
  return df[(df['14'] > t) | (df['14'] < 1)] 

def generate_budgeted_policies(num_antigen, num_pcr):
  '''
  Generate an array of possible testing policies with budgets for each type of test given by num_antigen and num_pcr
  '''
  
  antigen_days = combinations(range(0,  T), num_antigen)
  pcr_days = combinations(range(1, T + 1), num_pcr)
  
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
  
  for t in range(0, T + 1):
    aux = filter_symptomatic(aux, t)
    if t - 1 in D_a: 
      aux = filter_antigent_tested(aux, t - 1)
    if t - 1 in D_p:
      aux = filter_pcr_tested(aux, t - 1)
      
    num = count_infecting(aux, t)
    den = aux.shape[0]
    if den != 0:
      prob_inf += num / num_mc #den
    
  return prob_inf  

def evaluate_policy_lag(df, policy, t_l=0):
  '''
  Evaluate expected cost of a policy when there is a lag of t_l in identifying exposure.
  A policy consists in a tuple, where the first entry are the antigen result test days, and the second are the PCR result days.
  '''
  
  aux = df.copy()
  num_mc = df.shape[0]
  
  D_a = [d + t_l for d in policy[0]]
  D_p = [d + t_l for d in policy[1]]
  
  prob_inf = []  # (1): Expected number of infecting days
  work_days = [] # (2): Expected number of non-infecting days
  
  variance_pi = 0 # (3): Variance of (1)
  variance_wd = 0 # (4): Variance of (2)
  
  # (5): x% quantile of (1)
  quantile_pi = {.95: 0,
                 .90: 0,
                 .80: 0,
                 .50: 0,
                 .20: 0,
                 .10: 0,
                 .05: 0,}
  # (6): x% quantile of (2)
  quantile_wd = {.95: 0,
                 .90: 0,
                 .80: 0,
                 .50: 0,
                 .20: 0,
                 .10: 0,
                 .05: 0}
  
  for t in range(0, 14):
    aux = filter_symptomatic(aux, t)
    if t in D_a: 
      aux = filter_antigent_tested(aux, t)
    if t in D_p:
      aux = filter_pcr_tested(aux, t)
      
    den = aux.shape[0]
    if den != 0:
      num = count_infecting(aux, t)
      prob_inf.append(num / num_mc)
      work_days.append((den - num) / num_mc)
      
      variance_pi += prob_inf[t] * (1 - prob_inf[t]) 
      variance_wd += work_days[t] * (1 - work_days[t])
      for s in range(t):
        num_tau = count_infecting_corr(aux, t, s)
        variance_pi += 2 * (num_tau / num_mc - prob_inf[t] * prob_inf[s])
        variance_wd += 2 * ((den - num) / num_mc - work_days[t] * work_days[s])
      
      for q in quantile_pi:
        quantile_pi[q] = quantile_pi[q] + (prob_inf[t] >= 1 - q)
        quantile_wd[q] = quantile_wd[q] + (work_days[t] >= 1 - q)
      
  return sum(prob_inf), sum(work_days), variance_pi, variance_wd, quantile_pi, quantile_wd

def weight_policy_symptom_day(df, policy):
  '''
  Weight the cost of each policy by the (distribution) of the day in which the individual was infected by the index, where the distribution was already computed based on the conditional distribution of infected days given symptoms on day 0
  '''
  
  inf_day_pdf = {0: 0.1669353401405535,
                  1: 0.22514759942847082,
                  2: 0.26480058806134404,
                  3: 0.21587741781643335,
                  4: 0.10776722889318635,
                  5: 0.01918233381609943}
  inf_day_pdf = {0: 0.23550498586741955,
                  1: 0.18526076665827995,
                  2: 0.1491931526909608,
                  3: 0.12145630421492379,
                  4: 0.09797953217203573,
                  5: 0.07684305227089296,
                  6: 0.057105543597348894,
                  7: 0.03934966233862648,
                  8: 0.024580317017036127,
                  9: 0.012726683172475764}
  inf_day_pdf = {d: inf_day_pdf[d] / sum(inf_day_pdf.values()) for d in inf_day_pdf}
  
  prob_inf = []
  work_days = []
  
  variance_pi = []
  variance_wd = []
  
  quantile_pi = []
  quantile_wd = []

  for d in inf_day_pdf:
    pi, wd, v_pi, v_wd, q_pi, q_wd = evaluate_policy_lag(df, policy, t_l=d)
    
    prob_inf.append(pi)
    work_days.append(wd)
    
    variance_pi.append(v_pi)
    variance_wd.append(v_wd)
    
    quantile_pi.append(q_pi)
    quantile_wd.append(q_wd)
  
  prob_inf_w = sum([prob_inf[d] * inf_day_pdf[d] for d in inf_day_pdf])
  work_days_w = sum([work_days[d] * inf_day_pdf[d] for d in inf_day_pdf])
  
  variance_pi_w = sum([variance_pi[d] * inf_day_pdf[d] for d in inf_day_pdf])
  variance_wd_w = sum([variance_wd[d] * inf_day_pdf[d] for d in inf_day_pdf])
  
  quantile_pi_w = {q: sum([quantile_pi[d][q] * inf_day_pdf[d] for d in inf_day_pdf]) for q in q_pi}
  quantile_wd_w = {q: sum([quantile_wd[d][q] * inf_day_pdf[d] for d in inf_day_pdf]) for q in q_wd}
  
  return prob_inf_w, work_days_w, variance_pi_w, variance_wd_w, quantile_pi_w, quantile_wd_w

def evaluate_budgeted_policy(df, num_antigen, num_pcr):
  '''
  Evaluate all posible policies that use exactly num_antigen antigen test and num_pcr PCR tests.
  '''
  
  policies = generate_budgeted_policies(num_antigen, num_pcr)
  policy_cost = []
  
  for p in policies:
    prob_inf, work_days, variance_pi, variance_wd, quantile_pi, quantile_wd = weight_policy_symptom_day(df, p)
    aux_df = pd.DataFrame({'num_antigen': num_antigen,
                           'num_pcr': num_pcr,
                           'antigen_days': [[d for d in p[0]]],
                           'pcr_days': [[d - 1 for d in p[1]]],
                           'expected_infecting_days': prob_inf,
                           'variance_infecting_days': variance_pi,
                           'expected_non_infecting_days': work_days,
                           'variance_non_infecting_days': variance_wd})
    aux_q_pi = {'infecting_days_' + str(d): [quantile_pi[d]] for d in quantile_pi}
    aux_q_wd = {'non_infecting_days_' + str(d): [quantile_wd[d]] for d in quantile_wd}
    aux_df = aux_df.join(pd.DataFrame(aux_q_pi)).join(pd.DataFrame(aux_q_wd))
    
    policy_cost.append(aux_df)
  
  return policy_cost
  
# Main

def main2():
  df = pd.read_csv('data/viral_load_mc.csv').iloc[:int(2e5)]
  
  policy_space_limits = {'antigen': [0, 5], 
                         'pcr': [0, 2]}
  policy_space = {'antigen': range(policy_space_limits['antigen'][0],
                                   policy_space_limits['antigen'][1] + 1),
                  'pcr': range(policy_space_limits['pcr'][0],
                               policy_space_limits['pcr'][1] + 1)}

  policy_cost = []
  for num_pcr in policy_space['pcr']:
    for num_antigen in policy_space['antigen']:
      start_time = time.time()
      
      results = evaluate_budgeted_policy(df, num_antigen, num_pcr)
      policy_cost = policy_cost + results
      
      print(str((num_antigen, num_pcr)) + ' finished in ' + str((time.time() - start_time) / 60))
            
  pd.concat(policy_cost).to_csv('data/results_antigen_25_05_1.csv')

# Execute

start = time.time()    
main2()
print((time.time() - start) / 60)