#========================================================================
# Construct S-ABC plot
# Author: Yumi Takagi
# Creation date (dd.mm.yyyy): 16.01.2026
#========================================================================

#' @param effect_measure String indicating the label of the lower horizontal axis."difference in means", "risk difference", "log risk ratio", "log odds ratio", "log hazard ratio"
#' @param BF_type String indicating the BF hypotheses setting of plot. "simple", "composite"
#' @param estimate Numerical vector containing the estimate(s).
#' @param N1 Numerical vector containing the sample size(s). Must be equal the number of estimates.
#' @param N2 Numerical vector containing the sample size(s). Must be equal the number of estimates.
#' @param zte Numerical vector containing zero treatment effect. Must be equal null hypothesis of the null hypothesis statistical test.
#' @param ate Numerical vector containing assumed treatment effect. Must be equal treatment effect used in the sample size calculation.
#' @param mte (optional)Numerical vector containing minimum treatment effect. Do not write if the BF type is simple.
#' @param conf_level_main Numerical vector indicating the confidence level(s). Bust be between 0 and 1.
#' @param conf_level_sub (optional) Numerical vector indicating the confidence level(s). Bust be between 0 and 1.

SABCplot <- function(
  effect_measure = NULL
  , BF_type = NULL
  , estimate = NULL
  , stderr = NULL
  , N1 = NULL
  , N2 = NULL
  , zte = NULL
  , ate = NULL
  , mte = NULL
  , conf_level_main = NULL
  , conf_level_sub = NULL
  , low_axes = NULL
  , upp_axes = NULL
  , length_axes = NULL
){

  #---------------------------------------
  # Load required packages
  #---------------------------------------
  # require("ggplot2")
  # require("dplyr")

  #---------------------------------------
  # p-value calculation
  #---------------------------------------

  n_total <- N1+N2

  p_val <- function(Hyp, estimate, se){
  value <- abs((estimate - H)/se)
  p_value <- 2*(1 - pnorm(value))
  }

  p_value_zte <- p_val(zte, estimate, se)
  p_value_ate <- p_val(ate, estimate, se)
  p_value_mte <- p_val(mte, estimate, se)

  #---------------------------------------
  # ALL CI calculation
  #---------------------------------------

  CI_low <- function(estimate, se, confidence_level){
  lower <- estimate - qnorm((1 - (1 - confidence_level)/2), 0, 1)*se
  }
  CI_upp <- function(estimate, se, confidence_level){
  upper <- estimate + qnorm((1 - (1 - confidence_level)/2), 0, 1)*se
  }

  lower_main <- CI_low(estimate, se, conf_level_main)
  upper_main <- CI_upp(estimate, se, conf_level_main)
  lower_sub <- CI_low(estimate, se, conf_level_sub)
  upper_sub <- CI_upp(estimate, se, conf_level_sub)

  CC <- seq(0, (1-0.001), by = 0.0001)
  lower_list <- c()
  for (i in CC) {
   lower2 <- CI_low(estimate, se, i)
   lower_list <- c(lower_list, c(lower2))
   }
  upper_list <- c()
  for (i in CC) {
   upper2 <- CI_upp(estimate, se, i)
   upper_list <- c(upper_list, c(upper2))
   }

  CI_list <- data.frame(CC, lower_list, upper_list)

  #---------------------------------------
  # ALL CI calculation
  #---------------------------------------
  if (BF_type %in% "simple") {
    ln_LR <- function(HA, H0, estimate, se){
      likelihood_HA <- dnorm(estimate, mean = HA, sd = se)
      likelihood_H0 <- dnorm(estimate, mean = H0, sd = se)
      BF <- likelihood_H0/likelihood_HA
      ln_LR <- log(BF)
      }
    log_LR <- ln_LR(HA, H0, estimate, se)
  } else if (type %in% "composite") {


####mid####
ln_LR_mid <- function(HA, H0, estimate, mid, se, ver){
  df <- N1 + N2 - 2
  likelihood_HA <- 1-pnorm(mid, mean = estimate, sd = se)
  likelihood_H0 <- pnorm(mid, mean = estimate, sd = se)
  BF <- likelihood_H0/likelihood_HA
  ln_LR <- log(BF)
}
log_LR_mid <- ln_LR_mid(HA, H0, estimate, mid, se, ver)

#######################
#lnBF_0A category
#Caluculate the delta corresponding to lnBF_0A
#V_4(lnBF_0A = -4), V_3(lnBF_0A = -3), V_2(lnBF_0A = -2), V_1(lnBF_0A = -1),
#V0(lnBF_0A = 0), v1(lnBF_0A = 1), v2(lnBF_0A = 2), v3(lnBF_0A = 3), V4(lnBF_0A = 4)
#######################
####通常####
a <- seq(low_axes, upp_axes, by = 0.1)
b <- c()
for (i in a) {
  b2 <- ln_LR(HA, H0, i, se, ver)
  b <- c(b, c(b2))
}
lnBF_list <- data.frame(a, b)

V_num <- function(v, HA){
  lnBF_list$v_4 <- lnBF_list$b - (v)
  lnBF_list$v_4_line <- with(lnBF_list,(ifelse(v_4 < 0, 0, 1)))
  lnBF_list_v_4_minus <- subset(lnBF_list, lnBF_list$v_4_line == 0)
  lnBF_list_v_4_plus <- subset(lnBF_list, lnBF_list$v_4_line == 1)
  v_4_max <- ifelse(HA > 0,
                    max(lnBF_list_v_4_minus$v_4),
                    min(lnBF_list_v_4_plus$v_4))
  v_4_min <- ifelse(HA > 0,
                    min(lnBF_list_v_4_plus$v_4),
                    max(lnBF_list_v_4_minus$v_4))
  lnBF_list$v_4_flag_min <- with(lnBF_list, ifelse(v_4 == v_4_min, 1, 99))
  lnBF_list$v_4_flag_max <- with(lnBF_list, ifelse(v_4 == v_4_max, 1, 99))
  a_v_4_min <- lnBF_list$a[lnBF_list$v_4_flag_min == 1]
  a_v_4_max <- lnBF_list$a[lnBF_list$v_4_flag_max == 1]

  a_v_4 <- seq(a_v_4_min, a_v_4_max, by = 0.0001)
  b_v_4 <- c()

  for (i in a_v_4) {
    b2_v_4 <- ln_LR(HA, H0, i, se, ver)
    b_v_4 <- c(b_v_4, c(b2_v_4))
  }
  list_v_4 <- data.frame(a_v_4, b_v_4)

  list_v_4$c_v_4 <- list_v_4$b_v_4 - (v)
  te_v_4 <- which(abs(list_v_4$c_v_4 - 0) == min(abs(list_v_4$c_v_4 - 0)))
  V_4 <- list_v_4$a_v_4[te_v_4]
}

V_4 <- V_num(-4, HA)
V_3 <- V_num(-3, HA)
V_2 <- V_num(-2, HA)
V_1 <- V_num(-1, HA)
V0 <- V_num(0, HA)
V1 <- V_num(1, HA)
V2 <- V_num(2, HA)
V3 <- V_num(3, HA)
V4 <- V_num(4, HA)

####mid####
a <- seq(low_axes, upp_axes, by = 0.1)
b_mid <- c()
for (i in a) {
  b2_mid <- ln_LR_mid(HA, H0, i, mid, se, ver)
  b_mid <- c(b_mid, c(b2_mid))
}
lnBF_mid_list <- data.frame(a, b_mid)

V_num_mid <- function(v, HA){
  lnBF_mid_list$v_4 <- lnBF_mid_list$b - (v)
  lnBF_mid_list$v_4_line <- with(lnBF_mid_list,(ifelse(v_4 < 0, 0, 1)))
  lnBF_mid_list_v_4_minus <- subset(lnBF_mid_list, lnBF_mid_list$v_4_line == 0)
  lnBF_mid_list_v_4_plus <- subset(lnBF_mid_list, lnBF_mid_list$v_4_line == 1)
  v_4_max <- ifelse(HA > 0,
                    max(lnBF_mid_list_v_4_minus$v_4),
                    min(lnBF_mid_list_v_4_plus$v_4))
  v_4_min <- ifelse(HA > 0,
                    min(lnBF_mid_list_v_4_plus$v_4),
                    max(lnBF_mid_list_v_4_minus$v_4))
  lnBF_mid_list$v_4_flag_min <- with(lnBF_mid_list, ifelse(v_4 == v_4_min, 1, 99))
  lnBF_mid_list$v_4_flag_max <- with(lnBF_mid_list, ifelse(v_4 == v_4_max, 1, 99))
  a_v_4_min <- lnBF_mid_list$a[lnBF_mid_list$v_4_flag_min == 1]
  a_v_4_max <- lnBF_mid_list$a[lnBF_mid_list$v_4_flag_max == 1]

  a_v_4 <- seq(a_v_4_min, a_v_4_max, by = 0.0001)
  b_v_4 <- c()

  for (i in a_v_4) {
    b2_v_4 <- ln_LR_mid(HA, H0, i, mid, se, ver)
    b_v_4 <- c(b_v_4, c(b2_v_4))
  }
  list_v_4 <- data.frame(a_v_4, b_v_4)

  list_v_4$c_v_4 <- list_v_4$b_v_4 - (v)
  te_v_4 <- which(abs(list_v_4$c_v_4 - 0) == min(abs(list_v_4$c_v_4 - 0)))
  V_4 <- list_v_4$a_v_4[te_v_4]
}

V_4_mid <- V_num_mid(-4, HA)
V_3_mid <- V_num_mid(-3, HA)
V_2_mid <- V_num_mid(-2, HA)
V_1_mid <- V_num_mid(-1, HA)
V0_mid <- V_num_mid(0, HA)
V1_mid <- V_num_mid(1, HA)
V2_mid <- V_num_mid(2, HA)
V3_mid <- V_num_mid(3, HA)
V4_mid <- V_num_mid(4, HA)

#---------------------------------------
# s-value calculation
#---------------------------------------

p_value_HA <- p_val(HA, estimate, se, ver)
s_level_HA <- ifelse(log(p_value_HA)/log(0.5) < 7, log(p_value_HA)/log(0.5), 7)
p_value_mid <- p_val(mid, estimate, se, ver)
s_level_mid <- ifelse(log(p_value_mid)/log(0.5) < 7, log(p_value_mid)/log(0.5), 7)


######################
#s value
######################
s_val <- seq(0, 7, by = 0.001)
CC_cal <- c()
for (i in s_val) {
  CC_cal2 <- 1 - exp(i * log(0.5))
  CC_cal <- c(CC_cal, c(CC_cal2))
}
CC_cal_list <- data.frame(s_val, CC_cal)

CC_cal <- CC_cal_list$CC_cal
lower_list2 <- c()
for (i in CC_cal) {
  lower22 <- CI_low(HA, H0, estimate, se, i, ver)
  lower_list2 <- c(lower_list2, c(lower22))
}
upper_list2 <- c()
for (i in CC_cal) {
  upper22 <- CI_upp(HA, H0, estimate, se, i, ver)
  upper_list2 <- c(upper_list2, c(upper22))
}

CI_list2 <- data.frame(CC_cal, lower_list2, upper_list2)

CI_list3 <- merge(CC_cal_list, CI_list2, by = "CC_cal", all.x = T)

#######################
#Dataset in drawing the ABC plot
#######################
####通常####
logBF_scale <- seq(-4, 4, by = 1)
delta <- c(V_4, V_3, V_2, V_1, V0, V1, V2, V3, V4)
list_total <- data.frame(delta, logBF_scale)

####mid####
logBF_scale <- seq(-4, 4, by = 1)
delta_mid <- c(V_4_mid, V_3_mid, V_2_mid, V_1_mid, V0_mid, V1_mid, V2_mid, V3_mid, V4_mid)
list_total_mid <- data.frame(delta_mid, logBF_scale)

k <- length(estimate)
id <- k:1
xmin <- low_axes
xmax <- upp_axes
V_min <- ifelse(HA > 0, low_axes, upp_axes)
V_max <- ifelse(HA > 0, upp_axes, low_axes)
V_max2 <- V_max - length_axes

d_name <- paste("d","(", effect_measure, ")")
ci_name <- "confidence coefficient"
lnBF_name <- expression(paste("ln", BF[0][A],"(d)"))
H0HABF_nm <- c(expression(H[0],H[A]))
H0HABFmid_nm <- c("NTE","ATE","MTE")
s_name <- expression(paste(-log[2], "(1 - confidence coefficient)[S-value]"))

#s scale and CI scale
number <- 1:14
s_val_cal <- 0:7



ci_p_value_H0 <- ifelse(p_value_H0 < (1 - exp(7 * log(0.5))), 1 - p_value_H0, 1 - exp(7 * log(0.5)))
ci_p_value_HA <- ifelse(p_value_HA < (1 - exp(7 * log(0.5))), 1 - p_value_HA, 1 - exp(7 * log(0.5)))
ci_p_value_mid <- ifelse(p_value_mid  < (1 - exp(7 * log(0.5))), 1 - p_value_mid, 1 - exp(7 * log(0.5)))

ci_scale <- c(0.95, 0.99, confidence_level_main, ci_p_value_H0, ci_p_value_HA, ci_p_value_mid)

s_scale <- round(log(1 - ci_scale)/log(0.5), 2)
s_scale2 <- c(s_val_cal, s_scale)

CC_scale_cal <- 1 - exp(s_val_cal * log(0.5))
ci_scale2 <- c(CC_scale_cal, ci_scale)

bind <- data.frame(cbind(number, ci_scale2, s_scale2))
bind$ci_scale3 <- with(bind, ifelse(number ==  3 | number ==  4 | number ==  5 | number ==  6 | number ==  7 | number ==  8 | number ==  12 | number ==  13 | number ==  14, "", ci_scale2))
bind$s_scale3 <- with(bind, ifelse(number ==  9 | number ==  10 | number ==  11,  "", s_scale2))

bind <- bind[order(bind$ci_scale2),]

ci_scale_break <- bind$ci_scale2
ci_scale_label <- bind$ci_scale3
s_scale_break <- bind$s_scale2
s_scale_label <- bind$s_scale3

s_scale_max <- max(CI_list3$s_val)

s_level_H0 <- ifelse(log(p_value_H0)/log(0.5) < 7, log(p_value_H0)/log(0.5), 7)
s_level_main <- log(1 - confidence_level_main)/log(0.5)
s_level_sub <- log(1 - confidence_level_sub)/log(0.5)


Result <- data.frame(N_total = n_total, Estimate = estimate, SE = se, refH0 = H0, refHA = HA, mid = mid,
                     lower = lower_main, upper = upper_main, lower_sub = lower_sub, upper_sub = upper_sub,
                     lnBF = log_LR, lnBF_mid = log_LR_mid,
                     p_value_H0 = p_value_H0, p_value_HA = p_value_HA, s_level_H0 = s_level_H0,
                     y = id,
                     xmin = xmin, xmax = xmax, V_min = V_min, V_max = V_max, V0 = V0, length_axes = length_axes)

df1 <- data.frame(y = Inf, x = Result$refH0, group = "0")
df2 <- data.frame(y = Inf, x = Result$refHA, group = "1")
df3 <- data.frame(y = Inf, x = Result$mid, group = "2")

#######################
#Drawing the ABC plot
#######################
library(ggplot2)
library(ggpubr)
library(ggpattern)

#####信頼区間関数バージョン######
switch_CI <- ifelse(confidence_level_main == confidence_level_sub, F, T)
switch_BF <- ifelse(mid == 0, T, F)
switch_BF_mid <- ifelse(mid == 0, F, T)

ABCplot_CI <- ggplot(data = Result, aes(y = y))+
  geom_line(data = CI_list3, aes(y = s_val, x = lower_list2), linewidth = 1.1, color = "dodgerblue2")+
  geom_line(data = CI_list3, aes(y = s_val, x = upper_list2), linewidth = 1.1, color = "dodgerblue2")+
  geom_vline(aes(xintercept = estimate), color = "dodgerblue2", linetype = "dashed", linewidth = 1.1)+
  geom_vline(data = df1, aes(xintercept = x, color = group), linewidth = 1.5)+
  geom_vline(data = df2, aes(xintercept = x, color = group), linewidth = 1.5)+
  {if(switch_BF_mid)geom_vline(data = df3, aes(xintercept = x, color = group), linewidth = 1.2)}+
  {if(switch_BF)scale_color_manual(name="reference line", labels = H0HABF_nm, values=c("red", "springgreen4"))}+
  {if(switch_BF_mid)scale_color_manual(name="reference line", labels = H0HABFmid_nm, values=c("red", "springgreen4", "mediumorchid3"))}+
  geom_errorbar(mapping = aes(x = estimate, y = s_level_main, xmin = lower_main, xmax = upper_main), width = 0.2, linewidth = 1.1)+
  {if(switch_CI)geom_segment(aes(x = xmin, xend = xmax, y = s_level_sub, yend = s_level_sub), color = "grey38", linetype = "dotdash", linewidth = 1)}+
  geom_segment(aes(x = H0, xend = xmax, y = s_level_H0, yend = s_level_H0), color = "grey38", linetype = "dotted", linewidth = 1)+
  geom_segment(aes(x = HA, xend = xmax, y = s_level_HA, yend = s_level_HA), color = "grey38", linetype = "dotted", linewidth = 1)+
  geom_segment(aes(x = mid, xend = xmax, y = s_level_mid, yend = s_level_mid), color = "grey38", linetype = "dotted", linewidth = 1)+
  {if(switch_BF)scale_x_continuous(limits = c(xmin, xmax), breaks = seq(xmin, xmax, by = length_axes),
                                   sec.axis = dup_axis(breaks = delta, labels = logBF_scale, name = lnBF_name), expand = c(0, 0.001))}+
  {if(switch_BF_mid)scale_x_continuous(limits = c(xmin, xmax), breaks = seq(xmin, xmax, by = length_axes),
                                       sec.axis = dup_axis(breaks = delta_mid, labels = logBF_scale, name = lnBF_name), expand = c(0, 0.001))}+
  scale_y_reverse(breaks = s_scale_break, labels = ci_scale_label,
                  sec.axis = dup_axis(breaks = s_scale_break, labels = s_scale_label, name = s_name), expand = c(0.0005, 0))+
  labs(x = d_name, y = ci_name)+
  theme(
    panel.background = element_rect(fill = "transparent",color = NA),
    panel.grid.minor = element_line(color = NA),
    panel.grid.major = element_line(color = NA),
    plot.background = element_rect(fill = "transparent",color = NA) ,
    legend.position="bottom",
    legend.title = element_text(size=14, face="bold"),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(size = 14),
    axis.ticks = element_line(colour = "white"))+
  {if(switch_BF)annotate("rect", ymin = -Inf, ymax = Inf, xmin = V0, xmax = V1, alpha = 0.1, fill = "orangered3")}+
  {if(switch_BF)annotate("rect", ymin = -Inf, ymax = Inf, xmin = V1, xmax = V2, alpha = 0.2, fill = "orangered3")}+
  {if(switch_BF)annotate("rect", ymin = -Inf, ymax = Inf, xmin = V2, xmax = V3, alpha = 0.3, fill = "orangered3")}+
  {if(switch_BF)annotate("rect", ymin = -Inf, ymax = Inf, xmin = V3, xmax = V4, alpha = 0.4, fill = "orangered3")}+
  {if(switch_BF)annotate("rect", ymin = -Inf, ymax = Inf, xmin = V4, xmax = V_min, alpha = 0.55, fill = "orangered3")}+
  {if(switch_BF)annotate("rect", ymin = -Inf, ymax = Inf, xmin = V0, xmax = V_1, alpha = 0.1, fill = "seagreen")}+
  {if(switch_BF)annotate("rect", ymin = -Inf, ymax = Inf, xmin = V_1, xmax = V_2, alpha = 0.2, fill = "seagreen")}+
  {if(switch_BF)annotate("rect", ymin = -Inf, ymax = Inf, xmin = V_2, xmax = V_3, alpha = 0.3, fill = "seagreen")}+
  {if(switch_BF)annotate("rect", ymin = -Inf, ymax = Inf, xmin = V_3, xmax = V_4, alpha = 0.4, fill = "seagreen")}+
  {if(switch_BF)annotate("rect", ymin = -Inf, ymax = Inf, xmin = V_4, xmax = V_max, alpha = 0.55, fill = "seagreen")}+
  {if(switch_BF_mid)annotate("rect", ymin = -Inf, ymax = Inf, xmin = V0_mid, xmax = V1_mid, alpha = 0.1, fill = "orangered3")}+
  {if(switch_BF_mid)annotate("rect", ymin = -Inf, ymax = Inf, xmin = V1_mid, xmax = V2_mid, alpha = 0.2, fill = "orangered3")}+
  {if(switch_BF_mid)annotate("rect", ymin = -Inf, ymax = Inf, xmin = V2_mid, xmax = V3_mid, alpha = 0.3, fill = "orangered3")}+
  {if(switch_BF_mid)annotate("rect", ymin = -Inf, ymax = Inf, xmin = V3_mid, xmax = V4_mid, alpha = 0.4, fill = "orangered3")}+
  {if(switch_BF_mid)annotate("rect", ymin = -Inf, ymax = Inf, xmin = V4_mid, xmax = V_min, alpha = 0.55, fill = "orangered3")}+
  {if(switch_BF_mid)annotate("rect", ymin = -Inf, ymax = Inf, xmin = V0_mid, xmax = V_1_mid, alpha = 0.1, fill = "seagreen")}+
  {if(switch_BF_mid)annotate("rect", ymin = -Inf, ymax = Inf, xmin = V_1_mid, xmax = V_2_mid, alpha = 0.2, fill = "seagreen")}+
  {if(switch_BF_mid)annotate("rect", ymin = -Inf, ymax = Inf, xmin = V_2_mid, xmax = V_3_mid, alpha = 0.3, fill = "seagreen")}+
  {if(switch_BF_mid)annotate("rect", ymin = -Inf, ymax = Inf, xmin = V_3_mid, xmax = V_4_mid, alpha = 0.4, fill = "seagreen")}+
  {if(switch_BF_mid)annotate("rect", ymin = -Inf, ymax = Inf, xmin = V_4_mid, xmax = V_max, alpha = 0.55, fill = "seagreen")}

ABCplot_CI
print(ABCplot)

}

