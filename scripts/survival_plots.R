library(survival)
library(survminer)

data_usage <- read.csv('/Users/heathertealess/Downloads/gbm_tcrs.csv',
                       stringsAsFactors = FALSE)

data_usage <- subset(data_usage,
                     vital_status              != '--' &
                       days_to_last_follow_up    != '--')

data_usage$time <- with(data_usage,
                        ifelse(days_to_death != '--',
                               as.numeric(days_to_death),
                               as.numeric(days_to_last_follow_up))
)
data_usage$status <- with(data_usage,
                          ifelse(vital_status == 'alive', 0, 1)
)

threshold <- median(data_usage$TRBV4.1, na.rm = TRUE)
data_usage$gene_usage_status <- factor(
  ifelse(data_usage$TRBV4.1 > threshold, "high", "low"),
  levels = c("low","high")
)

fit <- survfit(Surv(time, status) ~ gene_usage_status, data = data_usage)
print(fit)

cox_model   <- coxph(Surv(time, status) ~ gene_usage_status, data = data_usage)
cox_summary <- summary(cox_model)

HR      <- signif(cox_summary$coefficients[1, "exp(coef)"],  2)
p_value <- signif(cox_summary$coefficients[1, "Pr(>|z|)"],   3)
CI_low  <- signif(cox_summary$conf.int[1, "lower .95"],       2)
CI_high <- signif(cox_summary$conf.int[1, "upper .95"],       2)

ggsurvplot(
  fit,
  data               = data_usage,
  title              = paste0(
    "TRBV4-1 in TCGA-GBM: HR=", HR,
    "; CI=[", CI_low, ":", CI_high, "]"
  ),
  pval               = TRUE,
  risk.table         = TRUE,
  surv.median.line   = "hv",
  palette            = c("red", "blue"),
  xlab               = "Time (days)",
  ylab               = "Survival Probability",
  legend.title       = "Expression",
  legend.labs        = c("Low", "High"),
  font.title         = c(16, "bold"),
  font.x             = 12,
  font.y             = 12,
  font.tickslab      = 11
)
