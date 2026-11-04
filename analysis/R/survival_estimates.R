write_table <- function(x, path, na = "") {
  readr::write_csv(as.data.frame(x), path, na = na)
}

# The baseline has these tables in write.csv's format, which differs from readr's.
write_base_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE)
}

# Type 2 to match the reported tables.
sample_quantile <- function(x, p) {
  as.numeric(quantile(x, p, na.rm = TRUE, names = FALSE, type = 2))
}

normal_95 <- qnorm(0.975)

number_at_risk <- function(times, months) {
  vapply(months, function(month) sum(times >= month), integer(1))
}

competing_event_label <- function(code) {
  if (code == 1L) "Local tumor progression" else "Death before local tumor progression"
}

# Normal-approximation CIs on the cmprsk variance, clipped to [0, 1].
cif_estimates <- function(data, cohort, group = "Overall", times = milestone_months) {
  incidence <- cmprsk::cuminc(data$cr_time_months, data$cr_status, cencode = 0)
  at_times <- cmprsk::timepoints(incidence, times = times)
  map_dfr(1:2, function(code) {
    key <- grep(paste0(" ", code, "$"), rownames(at_times$est), value = TRUE)
    observed <- length(key) > 0L
    estimate <- if (observed) as.numeric(at_times$est[key[1], ]) else rep(0, length(times))
    se <- if (observed) sqrt(pmax(as.numeric(at_times$var[key[1], ]), 0)) else rep(0, length(times))
    tibble(
      cohort = cohort,
      stratum = group,
      event = competing_event_label(code),
      time_months = times,
      n = nrow(data),
      events = sum(data$cr_status == code),
      estimate = estimate,
      lower = if (observed) pmax(0, estimate - normal_95 * se) else 0,
      upper = if (observed) pmin(1, estimate + normal_95 * se) else 0,
      number_at_risk = number_at_risk(data$cr_time_months, times),
      note = if (observed) {
        "Pointwise normal-approximation CIF intervals bounded at 0 and 1"
      } else {
        "No observed events of this cause; zero empirical estimate, not proof of zero risk"
      }
    )
  })
}

survival_endpoints <- list(
  OS = list(time = "os_time_months", event = "os_event"),
  Composite = list(time = "lpfs_time_months", event = "lpfs_event"),
  LocalFailure = list(time = "cr_time_months", event = "I(cr_status==1)"),
  CompetingDeath = list(time = "cr_time_months", event = "I(cr_status==2)")
)

survival_response <- function(endpoint) {
  columns <- survival_endpoints[[endpoint]]
  paste0("Surv(", columns$time, ",", columns$event, ")")
}

km_estimates <- function(data, cohort, endpoint, group = "Overall", times = milestone_months) {
  columns <- survival_endpoints[[endpoint]]
  fit <- survfit(as.formula(paste0(survival_response(endpoint), "~1")), data = data, conf.type = "log-log")
  at_times <- summary(fit, times = times, extend = TRUE)
  milestones <- tibble(
    cohort = cohort,
    stratum = group,
    endpoint = endpoint,
    time_months = times,
    n = nrow(data),
    events = sum(data[[columns$event]]),
    estimate = at_times$surv,
    lower = at_times$lower,
    upper = at_times$upper,
    number_at_risk = number_at_risk(data[[columns$time]], times)
  ) %>%
    mutate(across(c(estimate, lower, upper), ~ if_else(number_at_risk == 0, NA_real_, .x)))
  overall <- summary(fit)$table
  median_survival <- tibble(
    cohort = cohort,
    stratum = group,
    endpoint = endpoint,
    n = nrow(data),
    events = sum(data[[columns$event]]),
    median_months = unname(overall["median"]),
    lower = unname(overall["0.95LCL"]),
    upper = unname(overall["0.95UCL"])
  )
  list(milestones = milestones, median = median_survival)
}

wald_p_value <- function(coefficients, covariance) {
  statistic <- tryCatch(
    as.numeric(t(coefficients) %*% solve(covariance, coefficients)),
    error = function(e) NA_real_
  )
  pchisq(statistic, length(coefficients), lower.tail = FALSE)
}

assert_complete_covariates <- function(data, covariates, label) {
  if (any(!complete.cases(data[, covariates, drop = FALSE]))) {
    stop("Missing model covariates in ", label)
  }
}

fit_fine_gray <- function(data, label, covariates = model_covariates) {
  assert_complete_covariates(data, covariates, label)
  design <- model.matrix(reformulate(covariates), data = data)[, -1, drop = FALSE]
  fit <- cmprsk::crr(data$cr_time_months, data$cr_status, cov1 = design, failcode = 1, cencode = 0)
  se <- sqrt(diag(fit$var))
  estimates <- tibble(
    model = label,
    term = names(fit$coef),
    n = nrow(data),
    events = sum(data$cr_status == 1),
    competing_events = sum(data$cr_status == 2),
    estimate = exp(fit$coef),
    lower = exp(fit$coef - normal_95 * se),
    upper = exp(fit$coef + normal_95 * se),
    p_value = 2 * pnorm(abs(fit$coef / se), lower.tail = FALSE),
    converged = isTRUE(fit$converged),
    warnings = ifelse(isTRUE(fit$converged), "", "Fine-Gray did not converge")
  )
  list(fit = fit, table = estimates, mm = design)
}

fit_cox <- function(data, label, endpoint, covariates = model_covariates) {
  assert_complete_covariates(data, covariates, label)
  captured <- character()
  fit <- withCallingHandlers(
    coxph(
      reformulate(covariates, response = survival_response(endpoint)),
      data = data,
      ties = "efron",
      x = TRUE,
      y = TRUE
    ),
    warning = function(w) {
      captured <<- c(captured, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  coefficients <- coef(fit)
  se <- sqrt(diag(vcov(fit)))
  estimates <- tibble(
    model = label,
    term = names(coefficients),
    n = fit$n,
    events = fit$nevent,
    estimate = exp(coefficients),
    lower = exp(coefficients - normal_95 * se),
    upper = exp(coefficients + normal_95 * se),
    p_value = summary(fit)$coefficients[, "Pr(>|z|)"],
    converged = all(is.finite(coefficients)) && !length(captured),
    warnings = paste(captured, collapse = "; ")
  )
  schoenfeld <- as.data.frame(cox.zph(fit)$table) %>%
    rownames_to_column("term") %>%
    transmute(model = label, term, chisq, df, p_value = p)
  list(fit = fit, table = estimates, ph = schoenfeld)
}

# Firth because the safety outcomes are sparse.
fit_firth_safety <- function(data, outcome, include_probe = TRUE) {
  covariates <- c("tumor_size_group", "central_location", if (include_probe) "probe_count")
  fit <- logistf::logistf(reformulate(covariates, outcome), data = data, pl = TRUE)
  size_terms <- grep("^tumor_size_group", names(fit$coefficients))
  tibble(
    term = names(fit$coefficients)[-1],
    estimate = exp(fit$coefficients[-1]),
    lower = exp(fit$ci.lower[-1]),
    upper = exp(fit$ci.upper[-1]),
    p_value = fit$prob[-1],
    global_size_p = logistf::logistftest(fit, test = size_terms)$prob
  )
}
