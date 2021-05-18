library(tidyverse)

df = read.csv('results_17_05_4.csv')

df %>% mutate(num_antigen = paste('num_antigen:', num_antigen),
              num_pcr = paste('num_pcr:', num_pcr)) %>%
        ggplot(aes(x=expected_infecting_days,
                  y=expected_non_infecting_days,
                  color=interaction(num_antigen,	num_pcr))) +
        geom_point() +
        facet_grid(num_antigen	~ num_pcr, scales='free') +
        theme(legend.title = element_blank(),
              legend.position = "none")
 
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
     
library(plotly)
df_aux = df %>% arrange(expected_infecting_days, expected_non_infecting_days) %>%
        group_by(num_pcr, num_antigen) %>%
        summarize(expected_infecting_days = first(expected_infecting_days),
                  expected_non_infecting_days = first(expected_non_infecting_days),
                  p_A = first(antigen_days),
                  p_P = first(pcr_days)) %>%
        filter(num_pcr + num_antigen > 0)

fig = df_aux %>% plot_ly(text = ~paste('num_A = ', num_antigen, ', num_P = ', num_pcr, '\npi_A = ', p_A, '\npi_P = ', p_P),
                         x = ~expected_infecting_days,
                         y = ~expected_non_infecting_days,
                         color = ~interaction(num_antigen, num_pcr))

fig
