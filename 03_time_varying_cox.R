# 03_time_varying_cox.R
# Cox model with a treatment-by-log(time + 1) interaction term, for each
# endpoint x subgroup.
#
# Inputs : TALAPRO2_<endpoint>_<subgroup>_<arm>.csv (same folder)
# Outputs: result_time_varying.csv

suppressPackageStartupMessages(library(survival))

load_subgroup <- function(endpoint, subgroup) {
  t_f <- sprintf("TALAPRO2_%s_%s_TalaEnza.csv",    endpoint, subgroup)
  c_f <- sprintf("TALAPRO2_%s_%s_PlaceboEnza.csv", endpoint, subgroup)
  d   <- rbind(read.csv(t_f), read.csv(c_f))
  d$treat_num <- as.integer(d$treat == "TalaEnza")
  d
}

results <- data.frame()
for (sg in c("HRRnondef", "HRRdef")) {
  for (ep in c("OS", "rPFS")) {
    if (!file.exists(sprintf("TALAPRO2_%s_%s_TalaEnza.csv", ep, sg))) next
    d   <- load_subgroup(ep, sg)
    fit <- coxph(
      Surv(Survival.time, Status) ~ treat_num + tt(treat_num),
      tt   = function(x, t, ...) x * log(t + 1),
      data = d
    )
    s <- summary(fit)$coefficients
    results <- rbind(results, data.frame(
      endpoint      = ep,
      subgroup      = sg,
      treat_coef    = s["treat_num",     "coef"],
      treat_p       = s["treat_num",     "Pr(>|z|)"],
      ttreat_coef   = s["tt(treat_num)", "coef"],
      interaction_p = s["tt(treat_num)", "Pr(>|z|)"]
    ))
  }
}

write.csv(results, "result_time_varying.csv", row.names = FALSE)
message("wrote result_time_varying.csv")
