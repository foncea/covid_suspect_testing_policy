#!usr/bin/env python3

'''
  File name: montecarlo_probs.py
  Author:    PFA
  Date:      Fri May 14, 2021

  Using MonteCarlo integration, compute the probabilities needed in the DP for the optimal testing policy of suspect Covid-19 cases.
'''

import numpy as np
import pandas as pd

def generate_viral_load(t_0, t_p, t_f, v_p):
  '''
  Generate viral load curve V_t given the parameters
  '''
  
  