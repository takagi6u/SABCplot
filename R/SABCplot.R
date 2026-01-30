#========================================================================
# Construct S-ABC plot
# Author: Yumi Takagi
# Creation date (dd.mm.yyyy): 30.01.2026
#========================================================================

#' @title S-ABCplot
#' @description \code{SABCplot} The function SABCplot is the function that plots Figures 3 and 4.
#' @importFrom ggplot2, dplyr
#' @param effect_measure String indicating the label of the lower horizontal axis."difference in means", "risk difference", "log risk ratio", "log odds ratio", "log hazard ratio"
#' @param BF_type String indicating the BF hypotheses setting of plot. "simple", "composite"
#' @param estimate Numerical vector containing the estimate(s).
#' @param N1 Numerical vector containing the sample size(s). Must be equal the number of estimates.
#' @param N2 Numerical vector containing the sample size(s). Must be equal the number of estimates.
#' @param zte Numerical vector containing zero treatment effect. Must be equal null hypothesis of the null hypothesis statistical test.
#' @param ate Numerical vector containing assumed treatment effect. Must be equal treatment effect used in the sample size calculation.
#' @param mte (optional)Numerical vector containing minimum treatment effect. Do not write if the BF type is simple.
#' @param conf_levels Numerical vector indicating the confidence level(s). Bust be between 0 and 1.
#' @param low_axes Numerical vector indicating lower limits of the bottom horizontal axis.
#' @param upp_axes Numerical vector indicating upper limits of the bottom horizontal axis.
#' @param length_axes Numerical vector indicating scaling of the bottom horizontal axis.

#' @examples
#'=========================================================================================
#' Graphics in the paper Figure 3, 4
#'=========================================================================================
#'
#'adjestimate <- 16.7
#'adjlower <- 4.2
#'adjupper <- 29.2
#'N1 <- 95
#'N2 <- 95
#'adjdf <- N1+N2-4
#'adjse <- (adjestimate - 4.2)/qt((1 - (1 - 0.95)/2), adjdf)
#'estimate <- adjestimate
#'se <- adjse
#'--
#'Figure3
#'--
#'SABCplot(effect_measure = "difference in means"
#'         , BF_type = "composite"
#'         , estimate = estimate
#'         , stderr = se
#'         , N1 = N1
#'         , N2 = N2
#'         , zte = 0
#'         , ate = 58
#'         , mte = 8
#'         , conf_levels = c(0.95)
#'         , low_axes = -10
#'         , upp_axes = 70
#'         , length_axes = 10)
#'ggsave("ABCplot_mte8.jpeg", dpi = 300, width = 11, height = 9)
#'--
#'Figure4
#'--
#'SABCplot(effect_measure = "difference in means"
#'         , BF_type = "composite"
#'         , estimate = estimate
#'         , stderr = se
#'         , N1 = N1
#'         , N2 = N2
#'         , zte = 0
#'         , ate = 58
#'         , mte = 20
#'         , conf_levels = c(0.95)
#'         , low_axes = -10
#'         , upp_axes = 70
#'         , length_axes = 10)
#' ggsave("ABCplot_mte20.jpeg", dpi = 300, width = 11, height = 9)
#' @export
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
  , conf_levels = NULL
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
  # p-value and S-value calculations
  #---------------------------------------

  n_total <- N1+N2

  p_val <- function(hyp, estimate, stderr){
  value <- abs((estimate - hyp)/stderr)
  p_value <- 2*(1 - pnorm(value))
  }

  p_val_zte <- p_val(zte, estimate, stderr)
  p_val_ate <- p_val(ate, estimate, stderr)

  s_val_zte <- log(p_val_zte)/log(0.5)
  s_val_ate <- log(p_val_ate)/log(0.5)

  s_val_zte_cut <- ifelse(s_val_zte < 7, s_val_zte, 7)
  s_val_ate_cut <- ifelse(s_val_ate < 7, s_val_ate, 7)

  if (!is.null(mte)) {
    p_val_mte <- p_val(mte, estimate, stderr)
    s_val_mte <- log(p_val_mte)/log(0.5)
    s_val_mte_cut <- ifelse(s_val_mte < 7, s_val_mte, 7)
  } else {
    p_val_mte <- NA
    s_val_mte <- NA
    s_val_mte_cut <- NA
  }

  #---------------------------------------
  # All CI calculation and S-values
  #---------------------------------------

  calc_CI <- function(estimate, stderr, conf_levels) {
    s_levels <- log(1 - conf_levels)/log(0.5)
    z_score <- qnorm(1 - (1 - conf_levels) / 2)
    lower <- estimate - z_score * stderr
    upper <- estimate + z_score * stderr
    return(data.frame(conf_level = conf_levels, s_levels = s_levels, estimate = estimate, lower = lower, upper = upper))
  }

  CI_results <- calc_CI(estimate, stderr, conf_levels)

  calc_all_CI <- function(estimate, stderr, s_levels) {
    cc_cals <- 1 - exp(s_levels * log(0.5))
    z_score <- qnorm(1 - (1 - cc_cals) / 2)
    lower <- estimate - z_score * stderr
    upper <- estimate + z_score * stderr
    return(data.frame(s_level = s_levels, conf_level = cc_cals, lower = lower, upper = upper))
  }

  sval <- seq(0, 7, by = 0.0001)
  CC_list <- calc_all_CI(estimate, stderr, sval)

  #---------------------------------------
  # LnBF_0A_value calculation
  # lnBF_0A category calculation. Caluculate the delta corresponding to lnBF_0A.
  # V_4(lnBF_0A = -4), V_3(lnBF_0A = -3), V_2(lnBF_0A = -2), V_1(lnBF_0A = -1),
  # V0(lnBF_0A = 0), v1(lnBF_0A = 1), v2(lnBF_0A = 2), v3(lnBF_0A = 3), V4(lnBF_0A = 4)
  #---------------------------------------

  calc_log_BF <- function(BF_type, estimate, stderr, ate = NULL, zte = NULL, mte = NULL) {
    if (BF_type == "simple") {
      log_lik_H0 <- dnorm(estimate, mean = zte, sd = stderr, log = TRUE)
      log_lik_HA <- dnorm(estimate, mean = ate, sd = stderr, log = TRUE)
      return(log_lik_H0 - log_lik_HA)
    } else if (BF_type == "composite") {
      log_lik_H0 <- pnorm(mte, mean = estimate, sd = stderr, log.p = TRUE)
      log_lik_HA <- pnorm(mte, mean = estimate, sd = stderr, lower.tail = FALSE, log.p = TRUE)
      return(log_lik_H0 - log_lik_HA)
    } else {
      stop("BF_type must be 'simple' or 'composite'")
    }
  }

  log_BF <- calc_log_BF(BF_type = BF_type, estimate = estimate, stderr = stderr, ate = ate, zte = zte, mte = mte)


  find_estimate_at_v <- function(v, type, stderr, low_axes, upp_axes, ate=NULL, zte=NULL, mte=NULL) {
    target_func <- function(x) {
      calc_log_BF(BF_type = type, estimate = x, stderr = stderr, ate = ate, zte = zte, mte = mte) - v
    }

    res <- tryCatch({
      uniroot(target_func, interval = c(low_axes, upp_axes))$root
    }, error = function(e) return(NA))
    return(res)
  }

  logBFval <- -4:4  # -4, -3, ..., 4
  est_list <- sapply(logBFval, function(v) {
    find_estimate_at_v(
      v = v,
      type = BF_type, # "simple" or "composite"
      stderr = stderr,
      low_axes = low_axes,
      upp_axes = upp_axes,
      ate = ate, zte = zte, mte = mte
    )
  })

  names(est_list) <- paste0("V", logBFval)

  BF_list <- data.frame(logBFval, est_list)

  V_4 <- with(BF_list, est_list[1])
  V_3 <- with(BF_list, est_list[2])
  V_2 <- with(BF_list, est_list[3])
  V_1 <- with(BF_list, est_list[4])
  V0 <- with(BF_list, est_list[5])
  V1 <- with(BF_list, est_list[6])
  V2 <- with(BF_list, est_list[7])
  V3 <- with(BF_list, est_list[8])
  V4 <- with(BF_list, est_list[9])

  #---------------------------------------
  # x axis Scaling
  #---------------------------------------

   xmin <- low_axes
   xmax <- upp_axes

   ref_point <- if(!is.null(mte)) mte else if(!is.null(ate)) ate else 0
   V_min <- ifelse(ref_point > 0, low_axes, upp_axes)
   V_max <- ifelse(ref_point > 0, upp_axes, low_axes)

   #---------------------------------------
   # y axis Scaling
   #---------------------------------------

   s_val_int <- 0:7
   CC_cal <- c(0, 0.5, 0.95, 0.99)

   cc_p_val_zte <- ifelse((1 - p_val_zte) < (1 - exp(7 * log(0.5))), 1 - p_val_zte, 1 - exp(7 * log(0.5)))
   cc_p_val_ate <- ifelse((1 - p_val_ate) < (1 - exp(7 * log(0.5))), 1 - p_val_ate, 1 - exp(7 * log(0.5)))
   if (!is.null(mte)) {
     cc_p_val_mte <- ifelse((1 - p_val_mte) < (1 - exp(7 * log(0.5))), 1 - p_val_mte, 1 - exp(7 * log(0.5)))
   } else {
     cc_p_val_mte <- numeric(0)
   }

   cc_defalut_scale <- c(CC_cal, cc_p_val_zte, cc_p_val_ate, cc_p_val_mte)

   cc_s_val_int <- 1 - exp(s_val_int * log(0.5))

   cc_scale <- sort(unique(c(conf_levels, cc_defalut_scale, cc_s_val_int)))

   scale_master <- data.frame(ci_break = cc_scale)
   scale_master$s_break <- round(log(1 - scale_master$ci_break)/log(0.5), 2)

   ci_to_show <- c(0, 0.5, 0.95, 0.99, conf_levels)
   s_to_show_int <- 0:7

   s_to_show_hyp <- c(s_val_zte_cut, s_val_ate_cut)
   if (!is.null(s_val_mte_cut)) {
     s_to_show_hyp <- c(s_to_show_hyp, s_val_mte_cut)
   }

   scale_master$ci_label <- ifelse(scale_master$ci_break %in% ci_to_show,
                                   as.character(round(scale_master$ci_break, 3)), "")

   scale_master$s_label  <- ifelse(round(scale_master$s_break, 2) %in% round(c(s_to_show_int, s_to_show_hyp), 2),
                                   as.character(scale_master$s_break), "")

   ci_scale_break <- scale_master$ci_break
   ci_scale_label <- scale_master$ci_label
   s_scale_break <- scale_master$s_break
   s_scale_label <- scale_master$s_label

 #---------------------------------------
 # Dataset in drawing the ABC plot
 #---------------------------------------

  k <- length(estimate)
  id <- k:1

  df1 <- data.frame(y = Inf, x = zte, group = "0")
  df2 <- data.frame(y = Inf, x = ate, group = "1")
  if (!is.null(mte)) {
    df3 <- data.frame(y = Inf, x = mte, group = "2")
  } else {
    df3 <- NULL
  }

  d_name <- paste(effect_measure)
  ci_name <- "confidence coefficient"
  lnBF_name <- expression(paste("ln", BF[0][A],"(d)"))
  zate_nm <- c("ZTE", "ATE")
  zamte_nm <- c("ZTE","ATE","MTE")
  s_name <- expression(paste(-log[2], "(1 - confidence coefficient)[S-value]"))

  switch_BF_sim <- ifelse(BF_type == "simple", T, F)
  switch_BF_com <- ifelse(BF_type == "composite", T, F)

  Result <- data.frame(n_total = n_total
                       , estimate = estimate
                       , stderr = stderr
                       , zte = zte
                       , ate = ate
                       , mte = if(is.null(mte)) NA else mte
                       , log_BF = log_BF
                       , p_val_zte = p_val_zte
                       , s_val_zte_cut = s_val_zte_cut
                       , s_val_ate_cut = s_val_ate_cut
                       , s_val_mte_cut = if(is.null(s_val_mte_cut)) NA else s_val_mte_cut
                       , y = id
                       , xmin = xmin
                       , xmax = xmax
                       , V_min = V_min
                       , V_max = V_max
                       , length_axes = length_axes)

S_ABCplot <- ggplot2::ggplot(data = Result, mapping = ggplot2::aes(y = y))+
  ggplot2::geom_line(data = CC_list, ggplot2::aes(y = s_level, x = lower), linewidth = 1.1, color = "dodgerblue2")+
  ggplot2::geom_line(data = CC_list, ggplot2::aes(y = s_level, x = upper), linewidth = 1.1, color = "dodgerblue2")+
  ggplot2::geom_vline(ggplot2::aes(xintercept = estimate), color = "dodgerblue2", linetype = "dashed", linewidth = 1.1)+
  ggplot2::geom_vline(data = df1, ggplot2::aes(xintercept = x, color = group), linewidth = 1.5)+
  ggplot2:: geom_vline(data = df2, ggplot2::aes(xintercept = x, color = group), linewidth = 1.5)+
  {if(!is.null(mte))ggplot2::geom_vline(data = df3, ggplot2::aes(xintercept = x, color = group), linewidth = 1.2)}+
  ggplot2::scale_color_manual(name = "reference line",
     labels = if(!is.null(mte)) zamte_nm else zate_nm,
     values = if(!is.null(mte)) c("red", "springgreen4", "mediumorchid3") else c("red", "springgreen4")) +
  ggplot2::geom_errorbar(data = CI_results, ggplot2::aes(x = estimate, y = s_levels, xmin = lower, xmax = upper), width = 0.2, linewidth = 1.1)+
  ggplot2::geom_segment(ggplot2::aes(x = zte, xend = xmax, y = s_val_zte_cut, yend = s_val_zte_cut), color = "grey38", linetype = "dotted", linewidth = 1)+
  ggplot2::geom_segment(ggplot2::aes(x = ate, xend = xmax, y = s_val_ate_cut, yend = s_val_ate_cut), color = "grey38", linetype = "dotted", linewidth = 1)+
  {if(!is.null(mte))ggplot2::geom_segment(ggplot2::aes(x = mte, xend = xmax, y = s_val_mte_cut, yend = s_val_mte_cut), color = "grey38", linetype = "dotted", linewidth = 1)}+
  ggplot2::scale_x_continuous(limits = c(xmin, xmax), breaks = seq(xmin, xmax, by = length_axes),
                    sec.axis = ggplot2::dup_axis(breaks = BF_list$est_list, labels = BF_list$logBFval, name = lnBF_name), expand = c(0, 0.001))+
  ggplot2::scale_y_reverse(breaks = s_scale_break, labels = ci_scale_label,
                  sec.axis = ggplot2::dup_axis(breaks = s_scale_break, labels = s_scale_label, name = s_name), expand = c(0.0005, 0))+
  ggplot2::labs(x = d_name, y = ci_name)+
  ggplot2:: theme(
    panel.background = ggplot2::element_rect(fill = "transparent",color = NA),
    panel.grid.minor = ggplot2::element_line(color = NA),
    panel.grid.major = ggplot2::element_line(color = NA),
    plot.background = ggplot2::element_rect(fill = "transparent",color = NA) ,
    legend.position="bottom",
    legend.title = ggplot2::element_text(size=14, face="bold"),
    axis.text.x = ggplot2::element_text(size = 18),
    axis.text.y = ggplot2::element_text(size = 18),
    axis.title = ggplot2::element_text(size = 14),
    axis.ticks = ggplot2::element_line(colour = "white"))+
  ggplot2::annotate("rect", ymin = -Inf, ymax = Inf, xmin = V0, xmax = V1, alpha = 0.1, fill = "orangered3")+
  ggplot2::annotate("rect", ymin = -Inf, ymax = Inf, xmin = V1, xmax = V2, alpha = 0.2, fill = "orangered3")+
  ggplot2::annotate("rect", ymin = -Inf, ymax = Inf, xmin = V2, xmax = V3, alpha = 0.3, fill = "orangered3")+
  ggplot2::annotate("rect", ymin = -Inf, ymax = Inf, xmin = V3, xmax = V4, alpha = 0.4, fill = "orangered3")+
  ggplot2::annotate("rect", ymin = -Inf, ymax = Inf, xmin = V4, xmax = V_min, alpha = 0.55, fill = "orangered3")+
  ggplot2::annotate("rect", ymin = -Inf, ymax = Inf, xmin = V0, xmax = V_1, alpha = 0.1, fill = "seagreen")+
  ggplot2::annotate("rect", ymin = -Inf, ymax = Inf, xmin = V_1, xmax = V_2, alpha = 0.2, fill = "seagreen")+
  ggplot2::annotate("rect", ymin = -Inf, ymax = Inf, xmin = V_2, xmax = V_3, alpha = 0.3, fill = "seagreen")+
  ggplot2::annotate("rect", ymin = -Inf, ymax = Inf, xmin = V_3, xmax = V_4, alpha = 0.4, fill = "seagreen")+
  ggplot2::annotate("rect", ymin = -Inf, ymax = Inf, xmin = V_4, xmax = V_max, alpha = 0.55, fill = "seagreen")

print(S_ABCplot)
}


