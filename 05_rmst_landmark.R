# 05_rmst_landmark.R
# Restricted mean survival time differences and landmark Cox analyses
# (HRR-nondef subgroup).
#
# Inputs : TALAPRO2_<endpoint>_HRRnondef_<arm>.csv (same folder)
# Outputs: result_rmst.csv
#          result_landmark.csv

suppressPackageStartupMessages({
  library(survival)
  library(survRM2)
})

load_subgroup <- function(endpoint) {
  t_f <- sprintf("TALAPRO2_%s_HRRnondef_TalaEnza.csv",    endpoint)
  c_f <- sprintf("TALAPRO2_%s_HRRnondef_PlaceboEnza.csv", endpoint)
  d   <- rbind(read.csv(t_f), read.csv(c_f))
  d$treat <- factor(d$treat, levels = c("PlaceboEnza", "TalaEnza"))
  d$arm   <- as.integer(d$treat == "TalaEnza")
  d
}

tau_grid      <- list(OS = c(12, 24, 36, 48), rPFS = c(6, 12, 18, 24, 30))
landmark_grid <- list(OS = c(24, 36),         rPFS = c(12, 18))

# --- RMST -------------------------------------------------------------------
rmst_out <- data.frame()
for (ep in names(tau_grid)) {
  if (!file.exists(sprintf("TALAPRO2_%s_HRRnondef_TalaEnza.csv", ep))) next
  d <- load_subgroup(ep)
  for (tau in tau_grid[[ep]]) {
    r <- tryCatch(
      rmst2(time = d$Survival.time, status = d$Status, arm = d$arm, tau = tau),
      error = function(e) NULL
    )
    if (is.null(r)) next
    e <- r$unadjusted.result["RMST (arm=1)-(arm=0)", ]
    rmst_out <- rbind(rmst_out, data.frame(
      endpoint  = ep,
      tau       = tau,
      RMST_diff = e["Est."],
      LCL       = e["lower .95"],
      UCL       = e["upper .95"],
      p         = e["p"]
    ))
  }
}
write.csv(rmst_out, "result_rmst.csv", row.names = FALSE)
message("wrote result_rmst.csv")

# --- Landmark Cox -----------------------------------------------------------
lm_out <- data.frame()
for (ep in names(landmark_grid)) {
  if (!file.exists(sprintf("TALAPRO2_%s_HRRnondef_TalaEnza.csv", ep))) next
  d <- load_subgroup(ep)
  for (lm_t in landmark_grid[[ep]]) {
    sub <- d[d$Survival.time > lm_t, ]
    sub$time_lm <- sub$Survival.time - lm_t
    fit <- coxph(Surv(time_lm, Status) ~ treat, data = sub)
    s   <- summary(fit)
    lm_out <- rbind(lm_out, data.frame(
      endpoint       = ep,
      landmark_month = lm_t,
      n_remaining    = nrow(sub),
      n_events       = sum(sub$Status == 1),
      HR  = s$coefficients[1, "exp(coef)"],
      LCL = s$conf.int[1, "lower .95"],
      UCL = s$conf.int[1, "upper .95"],
      p   = s$coefficients[1, "Pr(>|z|)"]
    ))
  }
}
write.csv(lm_out, "result_landmark.csv", row.names = FALSE)
message("wrote result_landmark.csv")
