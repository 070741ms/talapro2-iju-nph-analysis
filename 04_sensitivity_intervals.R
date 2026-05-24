# 04_sensitivity_intervals.R
# Piecewise Cox regression under three alternative interval-specification
# schemes (HRR-nondef subgroup; underlies eFigure 5).
#
# Inputs : TALAPRO2_<endpoint>_HRRnondef_<arm>.csv (same folder)
# Outputs: result_piecewise_sensitivity.csv

suppressPackageStartupMessages(library(survival))

schemes <- list(
  OS = list(
    primary      = c(12, 24, 36, 48),
    sens_3mo     = seq(3, 60, by = 3),
    sens_6mo     = seq(6, 60, by = 6),
    sens_shifted = c(9, 18, 30, 42)
  ),
  rPFS = list(
    primary      = c(6, 12, 18, 24, 30),
    sens_3mo     = seq(3, 36, by = 3),
    sens_6mo     = seq(6, 36, by = 6),
    sens_shifted = c(4, 10, 16, 22)
  )
)

load_subgroup <- function(endpoint) {
  t_f <- sprintf("TALAPRO2_%s_HRRnondef_TalaEnza.csv",    endpoint)
  c_f <- sprintf("TALAPRO2_%s_HRRnondef_PlaceboEnza.csv", endpoint)
  d   <- rbind(read.csv(t_f), read.csv(c_f))
  d$treat <- factor(d$treat, levels = c("PlaceboEnza", "TalaEnza"))
  d
}

piecewise <- function(d, cut) {
  sp <- survSplit(Surv(Survival.time, Status) ~ ., data = d,
                  cut = cut, episode = "interval_id")
  lo  <- c(0, cut); hi <- c(cut, NA)
  lab <- ifelse(is.na(hi), sprintf("%g+", lo), sprintf("%g-%g", lo, hi))

  out <- data.frame()
  for (i in seq_along(lab)) {
    sub <- sp[sp$interval_id == i, ]
    ev  <- sum(sub$Status == 1)
    if (ev < 5 || length(unique(sub$treat[sub$Status == 1])) < 2) {
      out <- rbind(out, data.frame(
        interval = lab[i], n_at_risk = nrow(sub), n_events = ev,
        HR = NA_real_, LCL = NA_real_, UCL = NA_real_, p = NA_real_))
      next
    }
    f <- coxph(Surv(tstart, Survival.time, Status) ~ treat, data = sub)
    s <- summary(f)
    out <- rbind(out, data.frame(
      interval  = lab[i],
      n_at_risk = nrow(sub),
      n_events  = ev,
      HR  = s$coefficients[1, "exp(coef)"],
      LCL = s$conf.int[1, "lower .95"],
      UCL = s$conf.int[1, "upper .95"],
      p   = s$coefficients[1, "Pr(>|z|)"]
    ))
  }
  out
}

all_out <- data.frame()
for (ep in names(schemes)) {
  if (!file.exists(sprintf("TALAPRO2_%s_HRRnondef_TalaEnza.csv", ep))) next
  d <- load_subgroup(ep)
  for (sch in names(schemes[[ep]])) {
    o <- piecewise(d, schemes[[ep]][[sch]])
    o$endpoint <- ep
    o$scheme   <- sch
    all_out <- rbind(all_out, o)
  }
}

all_out <- all_out[, c("endpoint", "scheme", "interval",
                       "n_at_risk", "n_events", "HR", "LCL", "UCL", "p")]
write.csv(all_out, "result_piecewise_sensitivity.csv", row.names = FALSE)
message("wrote result_piecewise_sensitivity.csv")
