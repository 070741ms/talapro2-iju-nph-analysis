# 02_piecewise_cox.R
# Piecewise Cox regression at the primary interval scheme, for each
# endpoint x subgroup. Reference arm: PlaceboEnza.
#
# Inputs : TALAPRO2_<endpoint>_<subgroup>_<arm>.csv (same folder)
# Outputs: result_piecewise_primary.csv

suppressPackageStartupMessages(library(survival))

cuts <- list(
  OS_HRRnondef   = c(12, 24, 36, 48),
  rPFS_HRRnondef = c(6, 12, 18, 24, 30),
  OS_HRRdef      = c(12, 24, 36),
  rPFS_HRRdef    = c(6, 12, 18)
)

load_subgroup <- function(endpoint, subgroup) {
  t_f <- sprintf("TALAPRO2_%s_%s_TalaEnza.csv",    endpoint, subgroup)
  c_f <- sprintf("TALAPRO2_%s_%s_PlaceboEnza.csv", endpoint, subgroup)
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
for (key in names(cuts)) {
  parts <- strsplit(key, "_", fixed = TRUE)[[1]]
  ep <- parts[1]; sg <- parts[2]
  if (!file.exists(sprintf("TALAPRO2_%s_%s_TalaEnza.csv", ep, sg))) next
  o <- piecewise(load_subgroup(ep, sg), cuts[[key]])
  o$endpoint <- ep
  o$subgroup <- sg
  all_out <- rbind(all_out, o)
}

all_out <- all_out[, c("endpoint", "subgroup", "interval",
                       "n_at_risk", "n_events", "HR", "LCL", "UCL", "p")]
write.csv(all_out, "result_piecewise_primary.csv", row.names = FALSE)
message("wrote result_piecewise_primary.csv")
