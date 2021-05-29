#!usr/bin/env python3

'''
  File name: evaluate_isolation_policy_benchmark
  Author:    PFA
  Date:      Fri 28 May 21

  Evaluate the performance of the benchmark given by the policy that immediately isolates the suspected individual after the index is discovered to be infected.
'''

import numpy as np
import pandas as pd

import time

# Methods

def count_infecting(df, t):
  '''
  Count number of viral load curves in df that are infecting at day t
  '''
  
  c = str(t)
  aux = df.loc[df[c] >= 6]
  return aux.shape[0]

def filter_symptomatic(df, t):
  '''
  Filter out symptomatic individuals
  '''
  
  return df[(df['14'] > t) | (df['14'] < 1)] 

def evaluate_isolation(df, t_l=0):
  '''
  Count the proportion of viral loads curves in df that are infecting before the index is identified.
  '''
  
  aux = df.copy()
  
  num_infecting = 0
  for t in range(t_l):
    aux = filter_symptomatic(aux, t)
    num_infecting += count_infecting(aux, t)
    
  return num_infecting / df.shape[0]

def weight_policy_index_day(df, index='symptoms'):
  '''
  Weight the cost of each policy by the (distribution) of the day in which the individual was infected by the index, where the distribution was already computed based on the conditional distribution of infected days given symptoms on day 0 (or antigen at day 0), depending on which one was used to identify the index.
  '''

  if index == 'symptoms':
    inf_day_pdf = {0: 0.1669353401405535,
                    1: 0.22514759942847082,
                    2: 0.26480058806134404,
                    3: 0.21587741781643335,
                    4: 0.10776722889318635,
                    5: 0.01918233381609943}
  if index == 'antigen':
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

  inf_days = {d: evaluate_isolation(df, t_l=d) for d in inf_day_pdf}
  exp_inf_days = sum([inf_days[d] * inf_day_pdf[d] for d in inf_days])
    
  return exp_inf_days
  
# Main
 
def main():
  df = pd.read_csv('data/viral_load_mc.csv').iloc[:int(2e5)]

  index_case = ['symptoms', 'antigen']
  for i in index_case:
    exp_inf_days = weight_policy_index_day(df, index=i)
    print(i + ': ' + str(exp_inf_days))
  
# Execute

start = time.time()    
main()
print((time.time() - start) / 60) 