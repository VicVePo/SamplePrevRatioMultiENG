#' Sample Size and Power Calculation for Prevalence Ratio Studies
#' 
#' @param alpha Significance level
#' @param power Desired power (optional if n is provided)
#' @param n Sample size (optional if power is provided)
#' @param p0 Prevalence in unexposed group (required if method = "rho")
#' @param PR Prevalence Ratio (optional if p1 is provided)
#' @param p1 Prevalence in exposed group (optional if PR is provided)
#' @param CI_upper Upper confidence interval limit of PR
#' @param CI_lower Lower confidence interval limit of PR (optional if SE is provided)
#' @param SE Standard error of coefficient (optional if CI is provided)
#' @param n_previous Sample size of previous study (required for method "se")
#' @param r Ratio of group sizes (exposed/unexposed), default is 1
#' @param method Calculation method: "rho" or "se" (standard error)
#' @param rho_values Vector of rho values for sample size calculation (if method = "rho")
#' @return A data frame with results according to chosen method
#' @export
SamplePrevRatioMultiENG <- function(alpha, power = NULL, n = NULL,
                                    p0 = NULL, PR = NULL, p1 = NULL,
                                    CI_upper = NULL, CI_lower = NULL,
                                    SE = NULL, n_previous = NULL, r = 1, 
                                    method = "rho",
                                    rho_values = seq(0, 0.9, by = 0.1)) {
  
  # Verifications according to method
  if (method == "rho") {
    if (is.null(p0)) {
      stop("For 'rho' method, p0 is required")
    }
    if (is.null(PR) && is.null(p1)) {
      stop("For 'rho' method, either PR or p1 must be provided")
    }
    if (!is.null(PR) && !is.null(p1)) {
      stop("Provide either PR or p1, not both")
    }
  } else if (method == "se") {
    if (is.null(PR)) {
      stop("For 'se' method, PR is required")
    }
    if (is.null(SE) && (is.null(CI_upper) || is.null(CI_lower))) {
      stop("For 'se' method, either SE or both CI limits are required")
    }
    if (is.null(n_previous)) {
      stop("For 'se' method, n_previous is required")
    }
  } else {
    stop("Method must be either 'rho' or 'se'")
  }
  
  # Calculate SE if CI provided
  if (!is.null(CI_upper) && !is.null(CI_lower)) {
    SE <- (log(CI_upper) - log(PR)) / qnorm(0.975)
    cat("Calculated Standard Error:", SE, "
")
  }
  
  if (method == "rho") {
    # Calculate p1 if PR provided
    if (!is.null(PR)) {
      p1 <- PR * p0
    } else {
      PR <- p1/p0
    }
    
    q1 <- r / (1 + r)
    q0 <- 1 / (1 + r)
    z_alpha <- qnorm(1 - alpha/2)
    p_pooled <- (q1 * p1 + q0 * p0)
    
    if (!is.null(power)) {
      z_beta <- qnorm(power)
      calculate_n <- function(rho) {
        A <- z_alpha * sqrt(p_pooled * (1 - p_pooled) * (1/q1 + 1/q0))
        B <- z_beta * sqrt(p1 * (1 - p1) / q1 + p0 * (1 - p0) / q0)
        C <- (p1 - p0)^2
        N <- ((A + B)^2 / C) / (1 - rho^2)
        return(ceiling(N))
      }
      results <- data.frame(
        rho = rho_values,
        total_size = sapply(rho_values, calculate_n)
      )
      results$exposed_size <- ceiling(results$total_size * q1)
      results$unexposed_size <- ceiling(results$total_size * q0)
    } else {
      calculate_power <- function(rho) {
        A <- z_alpha * sqrt(p_pooled * (1 - p_pooled) * (1/q1 + 1/q0))
        C <- (p1 - p0)^2
        B <- sqrt(n * C * (1 - rho^2)) - A
        power <- pnorm(B / sqrt(p1 * (1 - p1) / q1 + p0 * (1 - p0) / q0))
        return(power)
      }
      results <- data.frame(
        rho = rho_values,
        power = sapply(rho_values, calculate_power)
      )
    }
  } else {
    # Standard error based method
    z_alpha <- qnorm(1 - alpha/2)
    if (!is.null(power)) {
      z_gamma <- qnorm(power)
      n <- ceiling((z_alpha + z_gamma)^2 * n_previous * SE^2 / (log(PR))^2)
      results <- data.frame(
        total_size = n
      )
      results$exposed_size <- ceiling(n * r/(1 + r))
      results$unexposed_size <- ceiling(n * 1/(1 + r))
    } else {
      z_gamma <- sqrt(n * (log(PR))^2 / (n_previous * SE^2)) - z_alpha
      power <- pnorm(z_gamma)
      results <- data.frame(
        power = power
      )
    }
  }
  return(results)
}

#' Cross-sectional Study Logistics
#'
#' @param final_n Final calculated sample size
#' @param rejection_rate Expected rejection rate (proportion between 0 and 1)
#' @param eligibility_rate Expected eligibility rate (proportion between 0 and 1)
#' @param subjects_per_day Number of subjects/records that can be processed per day
#' @param working_days_month Number of working days per month
#' @return A list with study logistics calculations
#' @export
cross_sectional_logistics <- function(final_n, 
                                    rejection_rate, 
                                    eligibility_rate, 
                                    subjects_per_day, 
                                    working_days_month) {
  
  # Input validation
  if (rejection_rate < 0 || rejection_rate >= 1) 
    stop("Rejection rate must be between 0 and 1")
  if (eligibility_rate <= 0 || eligibility_rate > 1) 
    stop("Eligibility rate must be between 0 and 1")
  if (subjects_per_day <= 0) 
    stop("Number of subjects per day must be positive")
  if (working_days_month <= 0) 
    stop("Number of working days per month must be positive")
  
  # Calculations
  n_evaluate <- final_n / (1 - rejection_rate)
  n_invite <- n_evaluate / eligibility_rate
  total_days <- n_invite / subjects_per_day
  duration_months <- total_days / working_days_month
  
  # Results
  results <- list(
    final_sample = final_n,
    evaluate_sample = ceiling(n_evaluate),
    invite_sample = ceiling(n_invite),
    recruitment_days = ceiling(total_days),
    recruitment_months = round(duration_months, 2)
  )
  
  # Print summary
  cat("
Study logistics summary:
")
  cat("----------------------------------------
")
  cat("Required final sample:", final_n, "
")
  cat("Subjects to evaluate:", ceiling(n_evaluate), 
      "(considering", rejection_rate*100, "% rejection rate)
")
  cat("Subjects to invite:", ceiling(n_invite), 
      "(considering", eligibility_rate*100, "% eligibility rate)
")
  cat("Days needed:", ceiling(total_days),
      "(evaluating", subjects_per_day, "subjects per day)
")
  cat("Months needed:", round(duration_months, 2),
      "(with", working_days_month, "working days per month)
")
  cat("----------------------------------------
")
  
  return(invisible(results))
}
