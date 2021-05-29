library(tidyverse)
library(plotly)

# Read Data
df = read.csv('data/results_antigen_25_05_1.csv')

# Plot 1: facet_grid with all points
df %>% mutate(num_antigen = paste('num_antigen:', num_antigen),
              num_pcr = paste('num_pcr:', num_pcr)) %>%
  ggplot(aes(x=expected_infecting_days,
             y=expected_non_infecting_days,
             color=interaction(num_antigen,	num_pcr))) +
  geom_point() +
  facet_grid(num_antigen	~ num_pcr, scales='free') +
  theme(legend.title = element_blank(),
        legend.position = "none")

# Plot 2: Only best points per budget
df %>% arrange(expected_infecting_days, expected_non_infecting_days) %>%
  group_by(num_pcr, num_antigen) %>%
  summarize(expected_infecting_days = first(expected_infecting_days),
            expected_non_infecting_days = first(expected_non_infecting_days)) %>%
  filter(num_pcr + num_antigen > 0) %>% 
  mutate(num_antigen = paste('num_antigen:', num_antigen),
         num_pcr = paste('num_pcr:', num_pcr)) %>% 
  ggplot(aes(x=expected_infecting_days,
             y=expected_non_infecting_days,
             color=interaction(num_antigen,	num_pcr))) +
  geom_point()

# Plot 3: Interactive with best points
df_aux = df %>% 
  arrange(expected_infecting_days, expected_non_infecting_days) %>%
  mutate(Policy_Budget = paste('Antigen: ', num_antigen, '- PCR: ', num_pcr),
         num_a = paste('Antigen:', num_antigen),
         num_p = paste('PCR:', num_pcr)) %>%
  group_by(num_pcr, num_antigen, Policy_Budget, num_a, num_p) %>%
  summarize(expected_infecting_days = round(first(expected_infecting_days), 2),
            expected_non_infecting_days = round(first(expected_non_infecting_days), 2),
            p_A = first(antigen_days),
            p_P = first(pcr_days)) 

df_aux %>% 
  plot_ly(text = ~paste('(A, P) = (', num_antigen, ', ', num_pcr, ')', sep=''),
          x = ~expected_infecting_days,
          y = ~expected_non_infecting_days,
          color = ~Policy_Budget,
          mode = 'line+markers',
          marker = list(size=8))

df_aux %>% 
  plot_ly(text = ~paste('(A, P) = (', num_antigen, ', ', num_pcr, ')', sep=''),
          x = ~expected_infecting_days,
          y = ~expected_non_infecting_days,
          color = ~num_a,
          symbol = ~num_p,
          symbols = c('0', 'circle', 'cross'),
          mode = 'markers',
          marker = list(size=12))

df_aux %>% 
  group_by(num_p) %>% 
  do(p = plot_ly(., text = ~paste('num_A = ', num_antigen, ', num_P = ', num_pcr, '\npi_A = ', p_A, '\npi_P = ', p_P),
          x = ~expected_infecting_days,
          y = ~expected_non_infecting_days,
          color = ~num_a,
          mode = 'markers',
          marker = list(size=8))) %>% 
  subplot(nrows = 3, shareX=TRUE)

# Plot 4: Frontier for all budgets together
pareto_df = read.csv('data/pareto_points_antigen_25_05_1.csv')

pareto_df %>% 
  mutate(Policy_Budget = paste('Antigen: ', num_antigen, '- PCR: ', num_pcr),
         num_a = paste('Antigen:', num_antigen),
         num_p = paste('PCR:', num_pcr)) %>%
  plot_ly(text = ~paste('num_A = ', num_antigen, ', num_P = ', num_pcr, '\n',
                        'pi_A = ', antigen_days, '\n',
                        'pi_P = ', pcr_days),
          x = ~expected_infecting_days,
          y = ~expected_non_infecting_days,
          color = ~Policy_Budget,
          mode = 'lines+markers',
          marker = list(size=8),
          line = list(shape='hv')
  )

pareto_df %>% 
  mutate(Policy_Budget = paste('Antigen: ', num_antigen, '- PCR: ', num_pcr),
         num_a = paste('Antigen:', num_antigen),
         num_p = paste('PCR:', num_pcr)) %>%
  plot_ly(text = ~paste('num_A = ', num_antigen, ', num_P = ', num_pcr, '\n',
                        'pi_A = ', antigen_days, '\n',
                        'pi_P = ', pcr_days),
          x = ~expected_infecting_days,
          y = ~expected_non_infecting_days,
          color = ~num_a,
          symbol = ~num_p,
          symbols = c('0', 'circle', 'cross'),
          mode = 'lines+markers',
          marker = list(size=8),
          line = list(shape='hv')
          )

# Plot 5: Percentile interval per policy

df %>% 
  arrange(expected_infecting_days, -expected_non_infecting_days) %>% 
  group_by(num_antigen, num_pcr) %>% 
  summarize(expected_infecting_days = first(expected_infecting_days),
            q80 = first(infecting_days_0.8),
            q20 = first(infecting_days_0.2),
            antigen_days = first(antigen_days),
            pcr_days = first(pcr_days)) %>% 
  mutate(budget_policy = paste('(', num_antigen, ', ', num_pcr, ')', sep=''),
         lower_error = expected_infecting_days - q20,
         upper_error = q80 - expected_infecting_days) %>% 
  ggplot(aes(x=budget_policy, y=expected_infecting_days)) +
  geom_point() +
  geom_errorbar(aes(ymin=q20, ymax=q80)) +
  xlab('(# Antigen, # PCR)')

aux = df %>% 
  arrange(expected_infecting_days, -expected_non_infecting_days) %>% 
  group_by(num_antigen, num_pcr) %>% 
  summarize(expected_infecting_days = first(round(expected_infecting_days, 2)),
            q80 = first(round(infecting_days_0.8, 2)),
            q20 = first(round(infecting_days_0.2, 2)),
            antigen_days = first(antigen_days),
            pcr_days = first(pcr_days)) %>% 
  mutate(budget_policy = paste('(', num_antigen, ', ', num_pcr, ')', sep=''),
         lower_error = expected_infecting_days - q20,
         upper_error = q80 - expected_infecting_days) %>% 
  plot_ly(type = 'scatter',
          text = ~paste('A: ', antigen_days, ', P: ', pcr_days),
          x = ~budget_policy,
          y = ~expected_infecting_days,
          error_y = ~list(symmetric=FALSE,
                          arrayminus= lower_error,
                          array= upper_error,
                          width = 10)
  )
aux  
