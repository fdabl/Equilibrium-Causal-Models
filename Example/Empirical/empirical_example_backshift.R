library(dplyr)
library(tidyr)
library(qgraph)
library(backShift)

# Available from Blanken et al. (2019) upon request
df <- read.csv('network_intervention_data_fixed.csv', sep = ';') %>% 
  # Rename pre and post
  rename_with(~sub('^pre_', 'week0_', .x), starts_with('pre_')) %>% 
  rename_with(~sub('^post_', 'week9_', .x), starts_with('post_'))

# RCT: 52 received CBTI, 52 no treatment
# T0-T1: Before treatment
# T2-T6: Treatment
# T7-T9: After treatment
df_sel_ISI <- df %>% select(proefpersoon, conditie, matches('^week\\d+_(ISI\\d+)$'))
df_sel_PHQ <- df %>% select(proefpersoon, conditie, matches('^week\\d+_(PHQ\\d+)$'))

df_ISI <- pivot_longer(
  df_sel_ISI, cols = starts_with('week'), names_to = c('time', '.value'),
  names_pattern = '(week\\d+)_(ISI\\d+)', values_to = 'ISI_score'
)

df_PHQ <- pivot_longer(
  df_sel_PHQ, cols = starts_with('week'), names_to = c('time', '.value'),
  names_pattern = '(week\\d+)_(PHQ\\d+)', values_to = 'PHQ_score'
) %>% 
  select(
    # Not used in the Blanken et al. paper
    -PHQ3, # Trouble falling asleep
    -PHQ10 # Impairment in daily life
  )

df_long <- full_join(df_ISI, df_PHQ, by = c('proefpersoon', 'conditie', 'time'))# %>% na.omit
df_long$time <- unname(sapply(df_long$time, function(x) as.numeric(strsplit(x, '')[[1]][5])))

df_missing <- df_long %>% 
  filter(rowSums(is.na(.)) >= 1)

# Around 12.2% of observations are missing
nrow(df_missing) / nrow(df_long)

to_numeric <- function(x) as.numeric(as.character(x))

# Let's average scores across pre, treatment, and post periods 
df_avg <- df_long %>%
  mutate(
    time = as.character(time),
    time_period = case_when(
      time %in% c('0', '1') ~ 1,
      time %in% c('2', '3', '4', '5', '6') ~ 2,
      time %in% c('7', '8', '9') ~ 3
    )
  ) %>% 
  group_by(proefpersoon, conditie, time_period) %>% 
  summarize(
    ISI1 = mean(ISI1, na.rm = TRUE),
    ISI2 = mean(ISI2, na.rm = TRUE),
    ISI3 = mean(ISI3, na.rm = TRUE),
    ISI4 = mean(ISI4, na.rm = TRUE),
    ISI5 = mean(ISI5, na.rm = TRUE),
    ISI6 = mean(ISI6, na.rm = TRUE),
    ISI7 = mean(ISI7, na.rm = TRUE),
    PHQ1 = mean(PHQ1, na.rm = TRUE),
    PHQ2 = mean(PHQ2, na.rm = TRUE),
    PHQ4 = mean(PHQ4, na.rm = TRUE),
    PHQ5 = mean(PHQ5, na.rm = TRUE),
    PHQ6 = mean(PHQ6, na.rm = TRUE),
    PHQ7 = mean(PHQ7, na.rm = TRUE),
    PHQ8 = mean(PHQ8, na.rm = TRUE),
    PHQ9 = mean(PHQ9, na.rm = TRUE)
  ) %>% ungroup()

df_avg_missing <- df_avg %>% 
  filter(rowSums(is.na(.)) >= 1)

# Only 6% are now missing 
# This is because some NAs are at certain time points but not others, and so averaging yields a value
nrow(df_avg_missing) / nrow(df_avg)

# There are three people for which the condition variable itself is missing — remove those
df_avg <- filter(df_avg, !is.na(conditie))

## Average across all weeks for the control group
## Average pre- and post-treatment weeks for treatment group
df_sel <- bind_rows(
  df_avg %>% 
    # Average treatment group as discussed
    filter(conditie == 1 & time_period != 2) %>% 
    group_by(proefpersoon, time_period) %>% 
    summarize(across(where(is.numeric), function(x) mean(x, na.rm = TRUE))),
  
  # Join with averaged control group as discussed
  df_avg %>% 
    filter(conditie == 0) %>% 
    group_by(proefpersoon) %>% 
    summarize(across(where(is.numeric), function(x) mean(x, na.rm = TRUE))),
) %>% ungroup %>% 
  # Variable names in line with Blanken et al. (2019)
  rename(
    DIS = ISI1,
    DMS = ISI2,
    EMA = ISI3,
    Dissat = ISI4,
    NIQoL = ISI5,
    Worry = ISI6,
    IDF = ISI7,
    LoI = PHQ1,
    `Dep Mood` = PHQ2,
    Fatigue = PHQ4,
    Appet = PHQ5,
    Worth = PHQ6,
    Con = PHQ7,
    `Psych Mot` = PHQ8,
    Sui = PHQ9
  )

df_sel_missing <- df_sel %>% 
  filter(rowSums(is.na(.)) >= 1)

# 3.8% of observations are missing (6 observations in total)
nrow(df_sel_missing) / nrow(df_sel)

# We remove those 6 observations
df_sel <- df_sel %>% 
  filter(rowSums(is.na(.)) == 0)

varmapping <- c(
  'DIS: Difficulty initiating sleep',
  'DMS: Difficulty maintaining sleep',
  'EMA: Early morning awakening',
  'Dissat: Dissatisfaction with sleep',
  'NIQoL: Noticeability of impaired quality of life',
  'Worry: Worry about sleep',
  'IDF: Interference with daily functioning',
  'LoI: Loss of interest',
  'DpM: Feeling down or depressed',
  'Fatigue: Fatigue',
  'Appet: Poor appetite or overeating',
  'Worth: Feelings of worthlessness',
  'Con: Concentration problems',
  'Psych Mot: Psychomotor agitation or retardation',
  'Sui: Suicidal thoughts'
)

X <- as.matrix(df_sel %>% ungroup() %>% select(-proefpersoon, -conditie, -time_period))
ind <- round(df_sel$time_period)

# Include stability selection with 500 iterations
res <- backShift(X, ind, threshold = 0.75, ev = 1, nsim = 500, verbose = TRUE)
# res$Ahat is point estimate
# res$AhatAdjacency give percentages of times an edge has been retained

cols <- RColorBrewer::brewer.pal(3, 'Set1')
pdf('Figures/ECM_empirical_example2.pdf', width = 20, height = 10)
par(mfrow = c(1, 2))

g1 <- qgraph(
  res$Ahat, layout = 'spring', edge.color = ifelse(res$Ahat < 0, cols[1], cols[2]),
  label.cex = 1.25, label.scale = FALSE, label.scale.equal = FALSE,
  labels = colnames(res$Ahat), legend = FALSE, vsize = 10, asize = 4
)

title(main = 'Estimates before stability selection', line = 2, font.main = 1, cex.main = 2)

par(mar = c(5, 1, 4, 2) + 0.1)
plot.new()
legend(
  x = 'left', y = 0.95,
  legend = varmapping,
  title = '',
  bty = 'n',
  cex = 1.75,
  y.intersp = 1.2, 
  xpd = NA
)
dev.off()