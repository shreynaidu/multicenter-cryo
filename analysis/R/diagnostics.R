# Wald test on the two contrasts orthogonal to the linear part of the basis.
fine_gray_nonlinearity <- function(dat) {
  spline_basis <- splines::ns(dat$tumor_size_cm, df = 3)
  design <- model.matrix(
    ~ spline_basis + age_per_10y + central_location + oxygen_any,
    data = dat
  )[, -1, drop = FALSE]
  fit <- crr(dat$cr_time_months, dat$cr_status, cov1 = design, failcode = 1, cencode = 0)
  spline_terms <- grep("^spline_basis", names(fit$coef))
  linear_component <- coef(lm(dat$tumor_size_cm ~ spline_basis))[-1]
  contrasts <- t(qr.Q(qr(matrix(linear_component, ncol = 1)), complete = TRUE)[, 2:3, drop = FALSE])
  estimate <- contrasts %*% fit$coef[spline_terms]
  variance <- contrasts %*% fit$var[spline_terms, spline_terms, drop = FALSE] %*% t(contrasts)
  tibble(
    endpoint = "LTP Fine-Gray",
    n = nrow(dat),
    nonlinearity_p = pchisq(as.numeric(t(estimate) %*% solve(variance, estimate)), 2, lower.tail = FALSE),
    converged = isTRUE(fit$converged),
    spline_df = 3
  )
}

cox_nonlinearity <- function(dat, endpoint, linear_fit) {
  fit <- coxph(
    reformulate(
      c("splines::ns(tumor_size_cm,df=3)", setdiff(model_covariates, "tumor_size_cm")),
      response = survival_response(endpoint)
    ),
    data = dat,
    ties = "efron"
  )
  likelihood_ratio <- 2 * (as.numeric(logLik(fit)) - as.numeric(logLik(linear_fit)))
  tibble(
    endpoint = endpoint,
    n = nrow(dat),
    nonlinearity_p = pchisq(likelihood_ratio, 2, lower.tail = FALSE),
    converged = all(is.finite(coef(fit))),
    spline_df = 3
  )
}

run_diagnostics <- function(dat, primary, out_dir) {
  results <- bind_rows(
    fine_gray_nonlinearity(dat),
    cox_nonlinearity(dat, "Composite", primary$Composite107$fit),
    cox_nonlinearity(dat, "OS", primary$OS107$fit)
  )
  write_table(results, file.path(out_dir, "diagnostics", "continuous_size_nonlinearity.csv"), na = "NA")
}
