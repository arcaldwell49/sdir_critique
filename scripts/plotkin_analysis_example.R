library(tidyverse)
library(here)
library(flextable)
library(gtsummary)
library(emmeans)
# Import Plotkin Data -----
source(
  here("scripts",
  "data_import.R")
)
# Import Gadbury functions ----
source(
  here("functions",
       "gadbury_no_covariate_functions.R")
)

# reccode ----
# 
df = df %>%
  mutate(
    group = factor(
      group,
      levels = 0:1,
      labels = c("LOAD",
                 "REPS")
    )
  )

# Summary statistics
df %>%
  select(group, x1rm_pre,
         x1rm_post,
         x1rm_change) %>%
  tbl_summary(by = group,
              label = list(
                x1rm_pre = "Pre",
                x1rm_post = 'Post',
                x1rm_change = "Change"
              ),
              statistic = list(all_continuous() ~ "{mean} ({sd})"),
              digits = list(everything() ~ c(1)))

# Model

model = lm(x1rm_change ~ x1rm_pre + group,
           data = df) 

ATE = pairs(emmeans(model, ~ group))

df %>%
  group_by(group) %>%
  summarize(pre_mean = mean(x1rm_pre,
                            na.rm = TRUE),
            pre_sd = sd(x1rm_pre,
                       na.rm = TRUE),
            post_mean = mean(x1rm_post,
                             na.rm = TRUE),
            post_sd = sd(x1rm_post,
                         na.rm = TRUE),
            change_mean = mean(x1rm_change,
                               na.rm = TRUE),
            change_sd = sd(x1rm_change,
                          na.rm = TRUE),
            .groups = 'drop')

var_d_mle <- sensitivity_analysis(subset(df, group == "REPS")$x1rm_change,
                                  subset(df, group == "LOAD")$x1rm_change,
                                  method = "mle")
var_d_mle %>%
  ggplot(aes(x=rho,y=sigma_d,
             ymin = sigma_d_lower,
             ymax = sigma_d_upper)) +
  geom_hline(yintercept = 0, color="darkred", alpha = .8) +
  geom_ribbon(fill = "grey",
              alpha = .2) +
  geom_line(linetype=4, #color = "red",
            linewidth = 1.1) +
  ggprism::theme_prism()

var_d_mle %>%
  ggplot(aes(x=rho,y=p_minus,
             ymin = p_minus_lower,
             ymax = p_minus_upper)) +
  geom_hline(yintercept = 0, color="darkred", alpha = .8) +
  geom_ribbon(fill = "grey",
              alpha = .2) +
  geom_line(linetype=4, #color = "red",
            linewidth = 1.1) +
  ggprism::theme_prism()
