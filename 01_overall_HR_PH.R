# 01_overall_HR_PH.R
# Overall hazard ratio and proportional-hazards diagnostics for each
# endpoint x subgroup. Reference arm: PlaceboEnza (HR < 1 favors TalaEnza).
#
# Inputs : TALAPRO2_<endpoint>_<subgroup>_<arm>.csv (same folder)
# Outputs: result_overall_HR.csv

suppressPackageStartupMessages(library(survival))

load_subgroup <- function(endpoint, subgroup) {
  t_f <- sprintf("TALAPRO2_%s_%s_TalaEnza.csv",    endpoint, subgroup)
  c_f <- sprintf("TALAPRO2_%s_%s_PlaceboEnza.csv", endpoint, subgroup)
  d   <- rbind(read.csv(t_f), read.csv(c_f))
  d$treat <- factor(d$treat, levels = c("PlaceboEnza", "TalaEnza"))
  d
}

results <- data.frame()
for (sg in c("HRRnondef", "HRRdef")) {
  for (ep in c("OS", "rPFS")) {
    if (!file.exists(sprintf("TALAPRO2_%s_%s_TalaEnza.csv", ep, sg))) next
    d   <- load_subgroup(ep, sg)
    fit <- coxph(Surv(Survival.time, Status) ~ treat, data = d)
    s   <- summary(fit)
    zph <- cox.zph(fit)
    results <- rbind(results, data.frame(
      endpoint     = ep,
      subgroup     = sg,
      n_total      = nrow(d),
      n_events     = sum(d$Status == 1),
      HR           = s$coefficients[1, "exp(coef)"],
      LCL          = s$conf.int[1, "lower .95"],
      UCL          = s$conf.int[1, "upper .95"],
      p            = s$coefficients[1, "Pr(>|z|)"],
      schoenfeld_p = zph$table["treat", "p"]
    ))
  }
}

write.csv(results, "result_overall_HR.csv", row.names = FALSE)
message("wrote result_overall_HR.csv")
