library(tidyverse)
library(here)
library(flextable)
library(gtsummary)
library(emmeans)
library(skedastic)
library(ggdist)
library(distributional)
library(patchwork)
# Import Plotkin Data -----
# # Data origin
# https://doi.org/10.7717/peerj.14142
source(
  here("scripts",
  "data_import.R")
)

# Import other functions -----
source(
  here("functions",
       "var_functions.R")
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
sum_table = df %>%
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

# Model -----

model = lm(x1rm_change ~ x1rm_pre + sex + group,
           data = df) 
 

## EMMEANS -----
# Get average treatment effect from estimated marginal means
# NOTE: differs from paper which used bootstrap with BCa CI methods
ATE = pairs(emmeans(model, ~ group)) %>% confint(level = .9)

# pairs(emmeans(model_int, ~ group| sex)) %>% confint(level = .9)
# Variance tests -------
# 
# Plot of change scores
# 
plot_delta = ggplot(df,
                    aes(x=group,y=x1rm_change)) +
  stat_dotsinterval(.width = .80)+
  labs(x = "",
       y = "Change Score (Post - Pre)") +
  ggprism::theme_prism()



# model.matrix(model)[,4]
aux_model = model.matrix(~group, contrasts = list(group = "contr.treatment"), data = df)
bp_test = breusch_pagan(model, 
              auxdesign = aux_model,
              koenker = FALSE)
## Plot of residuals
df$resid = residuals(model)

plot_resid = ggplot(df,
                    aes(x=group,y=resid)) +
  stat_dotsinterval(.width = .80)+
  ggprism::theme_prism()

# F-test for log ratio of variance
var_ratio_test = var.test(subset(df, group == "REPS")$x1rm_change,
                          subset(df, group == "LOAD")$x1rm_change)

lnvar_test = log_vr_test(
  sd1 = sd(subset(df, group == "REPS")$x1rm_change, na.rm = TRUE),
  sd2 = sd(subset(df, group == "LOAD")$x1rm_change, na.rm = TRUE),
  n1 = length((subset(df, group == "REPS")$x1rm_change)),
  n2 = length((subset(df, group == "LOAD")$x1rm_change))
)

sd_ir = diff_sdir_test(
  sd1 = sd(subset(df, group == "REPS")$x1rm_change, na.rm = TRUE),
  sd2 = sd(subset(df, group == "LOAD")$x1rm_change, na.rm =
             TRUE),
  n1 = length((subset(df, group == "REPS")$x1rm_change)),
  n2 = length((subset(df, group == "LOAD")$x1rm_change))
)
# assumed rho_xy
sd1 = sd(subset(df, group == "REPS")$x1rm_change, na.rm = TRUE)
sd2 = sd(subset(df, group == "LOAD")$x1rm_change, na.rm =TRUE)
sdir_rho_xy = (sd1^2 + sd2^2 - unname(sd_ir$estimate)^2) / (2 * sd1 * sd2)
# Heterogeneity of Treatment Effects ----------
var_d_mle <- sensitivity_analysis(control = subset(df, group == "LOAD")$x1rm_change,
                                  treatment = subset(df, group == "REPS")$x1rm_change,
                                  ATE = -1*ATE$estimate[1],
                                  ATE_SE = ATE$SE[1],
                                  method = "mle",
                                  lower.tail = FALSE)
plot_var_d_mle = var_d_mle %>%
  ggplot(aes(x=rho,y=sigma_d,
             ymin = sigma_d_lower,
             ymax = sigma_d_upper)) +
  geom_hline(yintercept = 0, color="darkred", alpha = .8) +
  geom_vline(xintercept = sdir_rho_xy, 
             linetype = "dotdash",
             color="forestgreen", alpha = .8) +
  geom_ribbon(fill = "grey",
              alpha = .2) +
  geom_line(linetype=4, #color = "red",
            linewidth = 1.1) +
  labs(x = expression(rho[LOAD*","*REPS]),
       y = expression(SD[TE])) +
  ggprism::theme_prism()

plot_pminus_mle = var_d_mle %>%
  ggplot(aes(x=rho,y=p_minus,
             ymin = p_minus_lower,
             ymax = p_minus_upper)) +
  geom_hline(yintercept = 0.5, color="darkred", alpha = .8) +
  geom_vline(xintercept = sdir_rho_xy, 
             linetype = "dotdash",
             color="forestgreen", alpha = .8) +
  geom_ribbon(fill = "grey",
              alpha = .2) +
  geom_line(linetype=4, #color = "red",
            linewidth = 1.1) +
  labs(x = expression(rho[LOAD*","*REPS]),
       y = expression(P[harmed])) +
  scale_y_continuous(labels = scales::percent,
                     breaks = seq(0.2, 1, .2))+
  ggprism::theme_prism()


library(ggdist)

indiv_level_plot = var_d_mle %>%
  mutate(ATE = ATE$estimate[[1]]) %>%
  ggplot(aes(x = rho, ydist = dist_normal(ATE, sigma_d))) +
  stat_lineribbon( .width = c(.5, .75, .9, .95, .98, .99)) +
  geom_vline(xintercept = sdir_rho_xy, 
             linetype = "dotdash",
             color="forestgreen", alpha = .8) +
  scale_fill_brewer()+
  ggprism::theme_prism()+
  labs(x = expression(rho[LOAD*","*REPS]),
       fill = "Quantile",
       y = "Difference in 1-RM back squat (kg)") 

df_long <- var_d_mle %>%
  pivot_longer(
    cols = c("sigma_d", "sigma_d_lower", "sigma_d_upper"),
    names_to = "sigma",
    values_to = "SD"
  ) %>%
  mutate(sigma_level = factor(case_when(
    sigma == "sigma_d" ~ .5,
    sigma == "sigma_d_upper" ~ .975,
    sigma == "sigma_d_lower" ~ .025
  )),ordered=TRUE)

indiv_level_plot2 = df_long %>%
  mutate(ATE = ATE$estimate[[1]]) %>%
  filter(rho %in% seq(-1,1,.5)) %>%
  filter(SD > 0) %>%
  ggplot(aes(x = rho, color = sigma_level, ydist = dist_normal(ATE, SD))) +
  geom_vline(xintercept = sdir_rho_xy, 
             linetype = "dotdash",
             color="forestgreen", alpha = .8) +
  stat_slab(fill=NA) +
 scale_color_viridis_d()+
  ggprism::theme_prism()+
  scale_x_continuous(breaks = seq(-1, 1, .5))+
  labs(x = expression(rho[LOAD*","*REPS]),
       color = "Confidence Interval Level",
       y = "Difference in 1-RM back squat (kg)") 


library(patchwork)

# p_final = plot_var_d_mle / plot_pminus_mle + plot_annotation(tag_levels = 'A')
# 
# ggsave(here("figure1_alt.png"),
#        height = 7.5,
#        width = 7.5)
# 
# p_final2 = indiv_level_plot / indiv_level_plot2 + plot_annotation(tag_levels = 'A')
# 
# ggsave(here("figure2.png"),
#        height = 7.5,
#        width = 7.5)

p_final3 = (plot_var_d_mle + indiv_level_plot) / (plot_pminus_mle + indiv_level_plot2) + plot_annotation(tag_levels = 'A')
ggsave(here("figure1.png"),
       height = 7.5,
       width = 10)
