library(tidyverse)
library(plotly)

# Read Data
df = read.csv('data/results_17_05_4.csv')

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
df_aux = df %>% arrange(expected_infecting_days, expected_non_infecting_days) %>%
  group_by(num_pcr, num_antigen) %>%
  summarize(expected_infecting_days = first(expected_infecting_days),
            expected_non_infecting_days = first(expected_non_infecting_days),
            p_A = first(antigen_days),
            p_P = first(pcr_days)) %>%
  filter(num_pcr + num_antigen > 0)

df_aux %>% plot_ly(text = ~paste('num_A = ', num_antigen, ', num_P = ', num_pcr, '\npi_A = ', p_A, '\npi_P = ', p_P),
                         x = ~expected_infecting_days,
                         y = ~expected_non_infecting_days,
                         color = ~interaction(num_antigen, num_pcr))

# Plot 4: Frontier for all budgets together
pareto_df = read.csv('data/pareto_points_18_05_1.csv')

pareto_df %>% 
  mutate(Policy_Budget = paste('Antigen: ', num_antigen, '- PCR: ', num_pcr)) %>%
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
          line = list(shape='vh')
          )
