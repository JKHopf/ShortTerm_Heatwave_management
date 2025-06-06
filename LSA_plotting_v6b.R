# Plotting results of local sensitivity analysis
#
# Jess Hopf Nov 2024
#


library(tidyverse)

# load in data
LSA_readin <- read_csv('LocalSA_outputs_v6b.csv', col_select=Para:Inc_persist)

LSA <- LSA_readin %>% 
  mutate(Paravalue = ordered(Paravalue, levels = c("L","B","H"), labels = c('Low','Baseline','High')),
         Scenario = ordered(Scenario, levels = c('Fishery closure','Urchin removal','Kelp seeding')),
         Para = ordered(Para, levels = c('urchin.PH',"kelp.rS",'kelp.c','kelp.g','urchin.RUstdv')))

# plotting

ggplot(LSA) +
  geom_col(aes(y = Inc_persist, x = MngtTime, fill = Scenario), position = 'dodge') +
  facet_grid(Paravalue~Para) +
  scale_fill_manual(values = c('#1f70b3','#684fa1','#218943'))+
  theme_bw() +
  ylab('Increase in probability of persistance (%)') +
  xlab('Start of management action (year relative to heatwave)')
  





