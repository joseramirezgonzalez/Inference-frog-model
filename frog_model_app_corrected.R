# ============================================================================
# FINITE-HORIZON FROG-MODEL INFERENCE EXPLORER
# Single-file R Shiny application
# Revised build: robust mathematical notation, a base-R score graph,
# an explicit K-T joint-growth diagnostic, and visible phase-color keys.
#
# Required packages:
#   install.packages(c("shiny", "bslib", "ggplot2", "plotly", "DT", "scales"))
#
# Run locally after saving this file as app.R inside a folder:
#   shiny::runApp("path/to/that/folder")
#
# It can also be launched directly with:
#   print(source("app.R"))
#
# The app implements the exact finite-horizon event-driven simulator,
# the origin-particle estimator, the activated-particle pooled estimator,
# and the constrained full-hazard maximum-likelihood estimator described in
# the accompanying article and Supplementary Material.
# ============================================================================

required_packages <- c("shiny", "bslib", "ggplot2", "plotly", "DT", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    paste0(
      "Missing required packages: ", paste(missing_packages, collapse = ", "),
      ". Install them with:\ninstall.packages(c(",
      paste(sprintf('"%s"', missing_packages), collapse = ", "), "))"
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(ggplot2)
  library(plotly)
  library(DT)
  library(scales)
})

options(shiny.maxRequestSize = 20 * 1024^2)

# ----------------------------------------------------------------------------
# Constants and presentation helpers
# ----------------------------------------------------------------------------

BETA_CRITICAL <- 0.5
PHASE_SURVIVAL <- "Positive-survival regime (S)"
PHASE_EXTINCTION <- "Almost-sure-extinction regime (E)"
PHASE_UNDETERMINED <- "Undetermined (U)"

APP_COLORS <- list(
  navy = "#102A43",
  blue = "#1565C0",
  teal = "#087F8C",
  green = "#16865B",
  amber = "#D97706",
  red = "#C2413B",
  purple = "#6D4CCF",
  slate = "#52606D",
  light = "#F5F8FB",
  border = "#D9E2EC"
)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

fmt_num <- function(x, digits = 4L) {
  if (length(x) == 0L) return(character(0))
  vapply(x, function(value) {
    if (is.na(value)) return("NA")
    if (is.infinite(value) && value > 0) return("∞")
    if (is.infinite(value) && value < 0) return("−∞")
    if (value == 0) return("0")
    av <- abs(value)
    if (av < 1e-3 || av >= 1e5) {
      formatC(value, format = "e", digits = max(2L, digits - 1L))
    } else {
      formatC(value, format = "fg", digits = digits, flag = "#")
    }
  }, character(1))
}

fmt_int <- function(x) {
  if (length(x) == 0L) return(character(0))
  vapply(x, function(value) {
    if (is.na(value)) return("NA")
    format(round(value), big.mark = ",", scientific = FALSE, trim = TRUE)
  }, character(1))
}

fmt_pct <- function(x, digits = 1L) {
  if (length(x) == 0L) return(character(0))
  vapply(x, function(value) {
    if (is.na(value)) return("NA")
    paste0(formatC(100 * value, format = "f", digits = digits), "%")
  }, character(1))
}

fmt_interval <- function(ci, digits = 4L) {
  if (length(ci) != 2L || any(is.na(ci))) return("Not available")
  paste0("[", fmt_num(ci[1], digits), ", ", fmt_num(ci[2], digits), "]")
}

phase_code <- function(ci, beta_c = BETA_CRITICAL) {
  if (length(ci) != 2L || any(is.na(ci))) return("NA")
  if (ci[2] < beta_c) return("S")
  if (ci[1] > beta_c) return("E")
  "U"
}

phase_label <- function(code) {
  if (length(code) == 0L) return(character(0))
  vapply(as.character(code), function(value) {
    switch(
      value,
      S = PHASE_SURVIVAL,
      E = PHASE_EXTINCTION,
      U = PHASE_UNDETERMINED,
      "Not available"
    )
  }, character(1))
}

phase_badge <- function(code) {
  class_name <- switch(
    as.character(code),
    S = "phase-badge phase-s",
    E = "phase-badge phase-e",
    U = "phase-badge phase-u",
    "phase-badge phase-na"
  )
  tags$span(class = class_name, phase_label(code))
}

phase_badge_html <- function(code) {
  class_name <- switch(
    as.character(code),
    S = "phase-badge phase-s",
    E = "phase-badge phase-e",
    U = "phase-badge phase-u",
    "phase-badge phase-na"
  )
  sprintf('<span class="%s">%s</span>', class_name, phase_label(code))
}

phase_color_key <- function() {
  tags$div(
    class = "phase-color-key",
    tags$div(
      class = "phase-key-item",
      tags$span(class = "phase-key-swatch swatch-s"),
      tags$span(tags$b("Green — Survival (S)"), tags$br(), HTML("The entire confidence region lies below <i>&beta;</i><sub>c</sub> = 0.5."))
    ),
    tags$div(
      class = "phase-key-item",
      tags$span(class = "phase-key-swatch swatch-u"),
      tags$span(tags$b("Amber / orange — Undetermined (U)"), tags$br(), HTML("The confidence region contains or crosses <i>&beta;</i><sub>c</sub> = 0.5."))
    ),
    tags$div(
      class = "phase-key-item",
      tags$span(class = "phase-key-swatch swatch-e"),
      tags$span(tags$b("Red — Extinction (E)"), tags$br(), HTML("The entire confidence region lies above <i>&beta;</i><sub>c</sub> = 0.5."))
    )
  )
}

true_phase_code <- function(beta) {
  if (is.na(beta)) return("NA")
  if (beta < BETA_CRITICAL) return("S")
  if (beta > BETA_CRITICAL) return("E")
  "U"
}

safe_softplus <- function(z) {
  pmax(z, 0) + log1p(exp(-abs(z)))
}

joint_growth_ratio <- function(K, T) {
  if (!is.finite(K) || !is.finite(T) || K <= 0 || T < 1) return(NA_real_)
  (T^4 * log(T^2 + 2)) / K
}

joint_growth_class <- function(rho) {
  if (!is.finite(rho)) return("control-caption")
  if (rho > 1) return("warning-note")
  if (rho > 0.1) return("method-note")
  "explanation"
}

# MathJax helpers.  Formulas are emitted as raw HTML delimiters and are
# typeset by the explicit MathJax 3 loader included in the UI below.
math_inline <- function(tex) {
  tags$span(class = "math-inline", HTML(paste0("\\(", tex, "\\)")))
}

math_display <- function(tex) {
  tags$div(class = "math-display", HTML(paste0("\\[", tex, "\\]")))
}

plotly_clean <- function(p) {
  plotly::config(
    p,
    displaylogo = FALSE,
    responsive = TRUE,
    modeBarButtonsToRemove = c(
      "lasso2d", "select2d", "autoScale2d", "toggleSpikelines"
    )
  )
}

# ----------------------------------------------------------------------------
# Exact simulator: Algorithm 1
# ----------------------------------------------------------------------------

simulate_one_realization <- function(beta, T, store_events = FALSE) {
  stopifnot(beta > 0, T >= 1)

  R <- integer(T)
  d <- integer(T)

  active_home <- 0L
  active_position <- 0L
  active_age <- 0L
  activated_sites <- 0L

  origin_jump <- NA_integer_

  dynamics_rows <- vector("list", T + 1L)
  dynamics_rows[[1L]] <- data.frame(
    time = 0L,
    alive = 1L,
    cumulative_activated = 1L,
    visible_transitions = 0L,
    jumps = 0L,
    deaths = 0L,
    new_activations = 0L,
    stringsAsFactors = FALSE
  )

  event_rows <- if (store_events) vector("list", max(1L, T * T)) else NULL
  event_index <- 0L
  stopped_early <- FALSE

  for (t in 0:(T - 1L)) {
    next_home <- integer(0)
    next_position <- integer(0)
    next_age <- integer(0)
    visited_positions <- integer(0)

    jumps_t <- 0L
    deaths_t <- 0L
    visible_t <- length(active_home)

    if (visible_t > 0L) {
      for (i in seq_len(visible_t)) {
        j <- active_age[i]
        idx <- j + 1L
        R[idx] <- R[idx] + 1L

        hazard <- beta / (beta + j + 1)
        dies <- stats::rbinom(1L, size = 1L, prob = hazard) == 1L

        if (active_home[i] == 0L && j == 0L && t == 0L) {
          origin_jump <- as.integer(!dies)
        }

        if (dies) {
          d[idx] <- d[idx] + 1L
          deaths_t <- deaths_t + 1L
          outcome <- "Death"
          position_after <- NA_integer_
        } else {
          direction <- sample(c(-1L, 1L), size = 1L)
          position_after <- active_position[i] + direction
          next_home <- c(next_home, active_home[i])
          next_position <- c(next_position, position_after)
          next_age <- c(next_age, j + 1L)
          visited_positions <- c(visited_positions, position_after)
          jumps_t <- jumps_t + 1L
          outcome <- if (direction < 0L) "Jump left" else "Jump right"
        }

        if (store_events) {
          event_index <- event_index + 1L
          event_rows[[event_index]] <- data.frame(
            calendar_time = t,
            home_label = active_home[i],
            age = j,
            position_before = active_position[i],
            outcome = outcome,
            position_after = position_after,
            death_hazard = hazard,
            stringsAsFactors = FALSE
          )
        }
      }
    }

    # Simultaneous activation is resolved only after every time-t transition.
    new_sites <- setdiff(unique(visited_positions), activated_sites)
    if (length(new_sites) > 0L) {
      next_home <- c(next_home, new_sites)
      next_position <- c(next_position, new_sites)
      next_age <- c(next_age, rep.int(0L, length(new_sites)))
      activated_sites <- c(activated_sites, new_sites)
    }

    dynamics_rows[[t + 2L]] <- data.frame(
      time = t + 1L,
      alive = length(next_home),
      cumulative_activated = length(activated_sites),
      visible_transitions = visible_t,
      jumps = jumps_t,
      deaths = deaths_t,
      new_activations = length(new_sites),
      stringsAsFactors = FALSE
    )

    active_home <- next_home
    active_position <- next_position
    active_age <- next_age

    if (length(active_home) == 0L && t < T - 1L) {
      for (future_time in (t + 2L):T) {
        dynamics_rows[[future_time + 1L]] <- data.frame(
          time = future_time,
          alive = 0L,
          cumulative_activated = length(activated_sites),
          visible_transitions = 0L,
          jumps = 0L,
          deaths = 0L,
          new_activations = 0L,
          stringsAsFactors = FALSE
        )
      }
      stopped_early <- TRUE
      break
    }
  }

  dynamics <- do.call(rbind, dynamics_rows[!vapply(dynamics_rows, is.null, logical(1))])
  events <- if (store_events && event_index > 0L) {
    do.call(rbind, event_rows[seq_len(event_index)])
  } else {
    data.frame(
      calendar_time = integer(0),
      home_label = integer(0),
      age = integer(0),
      position_before = integer(0),
      outcome = character(0),
      position_after = integer(0),
      death_hazard = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  extinction_candidates <- dynamics$time[dynamics$alive == 0L]
  extinction_time <- if (length(extinction_candidates) > 0L) {
    min(extinction_candidates)
  } else {
    NA_integer_
  }

  N <- R[1L]
  M <- R[1L] - d[1L]
  B <- sum(R)

  list(
    R = R,
    d = d,
    origin_jump = origin_jump,
    N = N,
    M = M,
    B = B,
    extinction_time = extinction_time,
    final_alive = tail(dynamics$alive, 1L),
    total_activated_by_T = tail(dynamics$cumulative_activated, 1L),
    dynamics = dynamics,
    events = events,
    stopped_early = stopped_early
  )
}

simulate_dataset <- function(K, T, beta, seed, detail_n = 20L, progress = NULL) {
  K <- as.integer(K)
  T <- as.integer(T)
  detail_n <- min(as.integer(detail_n), K)
  set.seed(as.integer(seed))

  R_total <- integer(T)
  d_total <- integer(T)
  origin_successes <- 0L
  summaries <- vector("list", K)
  dynamics <- vector("list", detail_n)
  events <- vector("list", detail_n)

  for (r in seq_len(K)) {
    one <- simulate_one_realization(beta, T, store_events = r <= detail_n)
    R_total <- R_total + one$R
    d_total <- d_total + one$d
    origin_successes <- origin_successes + one$origin_jump

    summaries[[r]] <- data.frame(
      realization = r,
      origin_jump = one$origin_jump,
      first_transition_rows_N = one$N,
      first_transition_jumps_M = one$M,
      visible_rows_B = one$B,
      extinction_time = one$extinction_time,
      alive_at_T = one$final_alive,
      activated_by_T = one$total_activated_by_T,
      stringsAsFactors = FALSE
    )

    if (r <= detail_n) {
      dynamics[[r]] <- one$dynamics
      events[[r]] <- one$events
    }

    if (!is.null(progress)) {
      progress(1 / K, paste("Simulating realization", r, "of", K))
    }
  }

  list(
    source = "Exact simulation",
    K = K,
    T = T,
    beta_reference = beta,
    seed = as.integer(seed),
    detail_n = detail_n,
    R = R_total,
    d = d_total,
    S_root = origin_successes,
    per_realization = do.call(rbind, summaries),
    dynamics = dynamics,
    events = events
  )
}

# ----------------------------------------------------------------------------
# Estimators and confidence sets
# ----------------------------------------------------------------------------

compute_origin_estimator <- function(S, K, alpha) {
  S <- as.integer(S)
  K <- as.integer(K)

  raw_q <- S / K
  raw_beta <- if (S == 0L) Inf else if (S == K) 0 else (K - S) / S

  smooth_q <- (S + 0.5) / (K + 1)
  smooth_beta <- (1 - smooth_q) / smooth_q

  q_lower <- if (S == 0L) 0 else stats::qbeta(alpha / 2, S, K - S + 1)
  q_upper <- if (S == K) 1 else stats::qbeta(1 - alpha / 2, S + 1, K - S)

  beta_ci <- c(
    1 / q_upper - 1,
    if (q_lower == 0) Inf else 1 / q_lower - 1
  )

  asymptotic_se_beta <- sqrt(
    smooth_beta * (smooth_beta + 1)^2 / K
  )

  list(
    S = S,
    K = K,
    raw_q = raw_q,
    raw_beta = raw_beta,
    smooth_q = smooth_q,
    estimate = smooth_beta,
    q_ci = c(q_lower, q_upper),
    beta_ci = beta_ci,
    phase = phase_code(beta_ci),
    asymptotic_se_beta = asymptotic_se_beta,
    lower_boundary = S == 0L,
    upper_boundary = S == K
  )
}

compute_pooled_estimator <- function(M, N, alpha) {
  M <- as.integer(M)
  N <- as.integer(N)

  raw_q <- M / N
  raw_beta <- if (M == 0L) Inf else if (M == N) 0 else (N - M) / M

  smooth_q <- (M + 0.5) / (N + 1)
  theta <- log((1 - smooth_q) / smooth_q)
  beta <- exp(theta)
  se_theta <- 1 / sqrt(N * smooth_q * (1 - smooth_q))
  z <- stats::qnorm(1 - alpha / 2)
  beta_ci <- exp(theta + c(-1, 1) * z * se_theta)

  list(
    M = M,
    N = N,
    raw_q = raw_q,
    raw_beta = raw_beta,
    smooth_q = smooth_q,
    theta = theta,
    estimate = beta,
    se_theta = se_theta,
    se_beta = beta * se_theta,
    beta_ci = beta_ci,
    phase = phase_code(beta_ci),
    lower_boundary = M == 0L,
    upper_boundary = M == N
  )
}

compute_hazard_mle <- function(R, d, alpha, beta_bounds) {
  R <- as.numeric(R)
  d <- as.numeric(d)
  T <- length(R)
  j <- 0:(T - 1L)

  if (length(d) != T) stop("R and d must have the same length.")
  if (any(!is.finite(R)) || any(!is.finite(d))) {
    stop("Risk and death counts must be finite.")
  }
  if (any(R < 0) || any(d < 0) || any(d > R)) {
    stop("Every age must satisfy 0 <= d_j <= R_j.")
  }
  if (length(beta_bounds) != 2L || any(!is.finite(beta_bounds)) ||
      beta_bounds[1] <= 0 || beta_bounds[2] <= beta_bounds[1]) {
    stop("The beta bounds must satisfy 0 < lower < upper.")
  }

  theta_bounds <- log(beta_bounds)
  B_plus <- sum(R)
  D_plus <- sum(d)
  if (B_plus <= 0) stop("The risk table contains no visible transitions.")

  loglik <- function(theta) {
    z <- theta - log(j + 1)
    sum(d * z - R * safe_softplus(z))
  }
  score <- function(theta) {
    D_plus - sum(R * stats::plogis(theta - log(j + 1)))
  }
  information <- function(theta) {
    h <- stats::plogis(theta - log(j + 1))
    sum(R * h * (1 - h))
  }

  score_lower <- score(theta_bounds[1])
  score_upper <- score(theta_bounds[2])

  if (score_lower <= 0) {
    theta_hat <- theta_bounds[1]
    boundary_flag <- "Lower constrained boundary"
  } else if (score_upper >= 0) {
    theta_hat <- theta_bounds[2]
    boundary_flag <- "Upper constrained boundary"
  } else {
    theta_hat <- stats::uniroot(
      score,
      interval = theta_bounds,
      tol = 1e-11,
      maxiter = 1000L
    )$root
    boundary_flag <- "Interior"
  }

  beta_hat <- exp(theta_hat)
  info_hat <- information(theta_hat)
  se_theta <- if (info_hat > 0) 1 / sqrt(info_hat) else NA_real_
  se_beta <- beta_hat * se_theta
  zcrit <- stats::qnorm(1 - alpha / 2)

  wald_ci <- if (boundary_flag == "Interior" && is.finite(se_theta)) {
    exp(theta_hat + c(-1, 1) * zcrit * se_theta)
  } else {
    c(NA_real_, NA_real_)
  }

  ll_hat <- loglik(theta_hat)
  chi_crit <- stats::qchisq(1 - alpha, df = 1)
  lr_equation <- function(theta) {
    2 * (ll_hat - loglik(theta)) - chi_crit
  }

  lower_lr <- if (theta_hat <= theta_bounds[1] + 1e-10) {
    theta_bounds[1]
  } else if (lr_equation(theta_bounds[1]) <= 0) {
    theta_bounds[1]
  } else {
    stats::uniroot(
      lr_equation,
      interval = c(theta_bounds[1], theta_hat),
      tol = 1e-11,
      maxiter = 1000L
    )$root
  }

  upper_lr <- if (theta_hat >= theta_bounds[2] - 1e-10) {
    theta_bounds[2]
  } else if (lr_equation(theta_bounds[2]) <= 0) {
    theta_bounds[2]
  } else {
    stats::uniroot(
      lr_equation,
      interval = c(theta_hat, theta_bounds[2]),
      tol = 1e-11,
      maxiter = 1000L
    )$root
  }

  lr_ci <- exp(c(lower_lr, upper_lr))

  direct_optimum <- stats::optimize(
    function(theta) -loglik(theta),
    interval = theta_bounds,
    tol = 1e-11
  )$minimum

  list(
    R = R,
    d = d,
    j = j,
    B_plus = B_plus,
    D_plus = D_plus,
    theta_bounds = theta_bounds,
    beta_bounds = beta_bounds,
    theta_hat = theta_hat,
    estimate = beta_hat,
    boundary_flag = boundary_flag,
    information = info_hat,
    se_theta = se_theta,
    se_beta = se_beta,
    wald_ci = wald_ci,
    lr_ci = lr_ci,
    wald_phase = phase_code(wald_ci),
    lr_phase = phase_code(lr_ci),
    loglik = loglik,
    score = score,
    information_function = information,
    loglik_hat = ll_hat,
    score_at_hat = score(theta_hat),
    score_lower = score_lower,
    score_upper = score_upper,
    chi_crit = chi_crit,
    direct_optimum = direct_optimum,
    optimization_difference = abs(theta_hat - direct_optimum)
  )
}

compute_all_estimators <- function(data, confidence_level, beta_bounds) {
  alpha <- 1 - confidence_level
  R <- as.numeric(data$R)
  d <- as.numeric(data$d)

  origin <- compute_origin_estimator(data$S_root, data$K, alpha)
  pooled <- compute_pooled_estimator(R[1] - d[1], R[1], alpha)
  hazard <- compute_hazard_mle(R, d, alpha, beta_bounds)

  list(
    origin = origin,
    pooled = pooled,
    hazard = hazard,
    confidence_level = confidence_level,
    alpha = alpha
  )
}

# ----------------------------------------------------------------------------
# Imported risk-table support and deterministic validation
# ----------------------------------------------------------------------------

read_risk_table <- function(path, T) {
  tab <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(tab) == 0L) stop("The uploaded CSV is empty.")

  lower_names <- tolower(trimws(names(tab)))
  age_col <- match("age", lower_names)
  r_col <- match("r", lower_names)
  d_col <- match("d", lower_names)

  if (any(is.na(c(age_col, r_col, d_col)))) {
    stop("The CSV must contain columns named age, R, and d.")
  }

  tab <- data.frame(
    age = as.integer(tab[[age_col]]),
    R = as.numeric(tab[[r_col]]),
    d = as.numeric(tab[[d_col]]),
    stringsAsFactors = FALSE
  )

  if (any(is.na(tab$age)) || any(is.na(tab$R)) || any(is.na(tab$d))) {
    stop("The age, R, and d columns must contain numeric values only.")
  }
  if (any(tab$age < 0L | tab$age >= T)) {
    stop("Every uploaded age must lie between 0 and T - 1.")
  }
  if (anyDuplicated(tab$age)) {
    stop("Each age may appear at most once in the uploaded table.")
  }
  if (any(tab$R != round(tab$R)) || any(tab$d != round(tab$d))) {
    stop("Risk and death counts must be nonnegative integers.")
  }

  full <- data.frame(age = 0:(T - 1L), R = 0, d = 0)
  full$R[match(tab$age, full$age)] <- tab$R
  full$d[match(tab$age, full$age)] <- tab$d

  if (any(full$R < 0) || any(full$d < 0) || any(full$d > full$R)) {
    stop("Every age must satisfy 0 <= d_j <= R_j.")
  }
  if (full$R[1] <= 0) stop("Age zero must have at least one risk row.")

  full
}

validate_dataset <- function(data, estimates) {
  R <- as.numeric(data$R)
  d <- as.numeric(data$d)
  K <- data$K
  T <- data$T
  N <- R[1]
  M <- R[1] - d[1]
  B <- sum(R)
  hazard <- estimates$hazard

  checks <- list(
    list(
      check = "Risk counts are nonnegative",
      status = all(R >= 0),
      detail = "<i>R</i><sub>j</sub> &ge; 0 for every visible age."
    ),
    list(
      check = "Deaths lie within the risk sets",
      status = all(d >= 0 & d <= R),
      detail = "0 &le; <i>d</i><sub>j</sub> &le; <i>R</i><sub>j</sub> for every visible age."
    ),
    list(
      check = "Age compression identity",
      status = identical(as.numeric(sum(R)), as.numeric(B)),
      detail = paste0("&sum;<sub>j=0</sub><sup>T&minus;1</sup> <i>R</i><sub>K,T,j</sub> = <i>B</i><sub>K,T</sub> = ", fmt_int(B), ".")
    ),
    list(
      check = "Age-zero identity",
      status = identical(as.numeric(R[1]), as.numeric(N)),
      detail = paste0("<i>R</i><sub>K,T,0</sub> = <i>N</i><sub>K,T</sub> = ", fmt_int(N), ".")
    ),
    list(
      check = "First-jump identity",
      status = identical(as.numeric(R[1] - d[1]), as.numeric(M)),
      detail = paste0("<i>R</i><sub>K,T,0</sub> &minus; <i>d</i><sub>K,T,0</sub> = <i>M</i><sub>K,T</sub> = ", fmt_int(M), ".")
    ),
    list(
      check = "Origin outcomes are compatible with age-zero totals",
      status = data$S_root <= M && (K - data$S_root) <= d[1],
      detail = "Origin jumps are a subset of age-zero jumps, and origin deaths are a subset of age-zero deaths."
    ),
    list(
      check = "Nested age-risk sets",
      status = if (T <= 1L) TRUE else all(R[-1L] <= R[-T] - d[-T]),
      detail = "For <i>j</i> &lt; <i>T</i>&minus;1, <i>R</i><sub>j+1</sub> &le; <i>R</i><sub>j</sub> &minus; <i>d</i><sub>j</sub>; later-age rows require an earlier observed jump."
    ),
    list(
      check = "Global finite-speed risk bound",
      status = B >= K && B <= K * T^2,
      detail = paste0("<i>K</i> &le; <i>B</i><sub>K,T</sub> &le; <i>K</i><i>T</i><sup>2</sup>; observed <i>B</i><sub>K,T</sub> = ", fmt_int(B), ".")
    ),
    list(
      check = "Global age-zero bound",
      status = N >= K && N <= K * (2 * T - 1),
      detail = paste0("<i>K</i> &le; <i>N</i><sub>K,T</sub> &le; <i>K</i>(2<i>T</i>&minus;1); observed <i>N</i><sub>K,T</sub> = ", fmt_int(N), ".")
    ),
    list(
      check = "Positive observed information",
      status = is.finite(hazard$information) && hazard$information > 0,
      detail = paste0("&#119973;<sub>K,T</sub>(&theta;&#770;) = ", fmt_num(hazard$information), ".")
    ),
    list(
      check = "Root-solver / direct-maximizer agreement",
      status = hazard$optimization_difference < 1e-6,
      detail = paste0(
        "|&theta;&#770;<sub>score</sub> &minus; &theta;&#770;<sub>direct</sub>| = ",
        formatC(hazard$optimization_difference, format = "e", digits = 2), "."
      )
    )
  )

  if (!is.null(data$per_realization)) {
    pr <- data$per_realization
    checks <- c(
      checks,
      list(
        list(
          check = "Per-realization age-zero finite-speed bounds",
          status = all(pr$first_transition_rows_N >= 1 &
                         pr$first_transition_rows_N <= 2 * T - 1),
          detail = "1 &le; <i>N</i><sub>r,T</sub> &le; 2<i>T</i> &minus; 1 for every simulated realization."
        ),
        list(
          check = "Per-realization visible-row bounds",
          status = all(pr$visible_rows_B >= 1 & pr$visible_rows_B <= T^2),
          detail = "1 &le; <i>B</i><sub>r,T</sub> &le; <i>T</i><sup>2</sup> for every simulated realization."
        )
      )
    )
  }

  out <- do.call(
    rbind,
    lapply(checks, function(x) {
      data.frame(
        Check = x$check,
        Status = if (isTRUE(x$status)) "PASS" else "FAIL",
        Detail = x$detail,
        stringsAsFactors = FALSE
      )
    })
  )

  out$Status <- ifelse(
    out$Status == "PASS",
    '<span class="validation-pass">PASS</span>',
    '<span class="validation-fail">FAIL</span>'
  )
  out
}

# ----------------------------------------------------------------------------
# Repeated-sampling explorer
# ----------------------------------------------------------------------------

run_repeated_sampling <- function(B, K, T, beta, confidence_level,
                                  beta_bounds, seed, progress = NULL) {
  B <- as.integer(B)
  alpha <- 1 - confidence_level
  rows <- vector("list", B)

  set.seed(as.integer(seed))
  dataset_seeds <- sample.int(.Machine$integer.max, B)

  for (b in seq_len(B)) {
    dat <- simulate_dataset(
      K = K,
      T = T,
      beta = beta,
      seed = dataset_seeds[b],
      detail_n = 0L,
      progress = NULL
    )
    est <- compute_all_estimators(dat, confidence_level, beta_bounds)

    origin_ci <- est$origin$beta_ci
    pooled_ci <- est$pooled$beta_ci
    wald_ci <- est$hazard$wald_ci
    lr_ci <- est$hazard$lr_ci

    rows[[b]] <- data.frame(
      replication = b,
      root_estimate = est$origin$estimate,
      pooled_estimate = est$pooled$estimate,
      hazard_estimate = est$hazard$estimate,
      root_cover = beta >= origin_ci[1] && beta <= origin_ci[2],
      pooled_cover = beta >= pooled_ci[1] && beta <= pooled_ci[2],
      hazard_wald_cover = if (any(is.na(wald_ci))) NA else beta >= wald_ci[1] && beta <= wald_ci[2],
      hazard_lr_cover = beta >= lr_ci[1] && beta <= lr_ci[2],
      root_length = origin_ci[2] - origin_ci[1],
      pooled_length = diff(pooled_ci),
      hazard_wald_length = if (any(is.na(wald_ci))) NA else diff(wald_ci),
      hazard_lr_length = diff(lr_ci),
      root_phase = est$origin$phase,
      pooled_phase = est$pooled$phase,
      hazard_wald_phase = est$hazard$wald_phase,
      hazard_lr_phase = est$hazard$lr_phase,
      hazard_boundary = est$hazard$boundary_flag != "Interior",
      N = est$pooled$N,
      B_rows = est$hazard$B_plus,
      stringsAsFactors = FALSE
    )

    if (!is.null(progress)) {
      progress(1 / B, paste("Repeated-sampling data set", b, "of", B))
    }
  }

  result <- do.call(rbind, rows)
  true_phase <- true_phase_code(beta)

  summarize_method <- function(name, estimate, cover, interval_length, phase) {
    finite_est <- is.finite(estimate)
    finite_len <- is.finite(interval_length)
    correct_phase <- if (true_phase %in% c("S", "E")) {
      mean(phase == true_phase, na.rm = TRUE)
    } else {
      NA_real_
    }
    decisive <- mean(phase %in% c("S", "E"), na.rm = TRUE)

    data.frame(
      Method = name,
      Bias = mean(estimate[finite_est] - beta, na.rm = TRUE),
      RMSE = sqrt(mean((estimate[finite_est] - beta)^2, na.rm = TRUE)),
      Median = stats::median(estimate[finite_est], na.rm = TRUE),
      Coverage = mean(cover, na.rm = TRUE),
      Mean_finite_length = mean(interval_length[finite_len], na.rm = TRUE),
      Unbounded_or_unavailable = mean(!finite_len),
      Correct_phase = correct_phase,
      Decisive_rate = decisive,
      stringsAsFactors = FALSE
    )
  }

  summary <- rbind(
    summarize_method(
      "Origin: exact Clopper–Pearson",
      result$root_estimate,
      result$root_cover,
      result$root_length,
      result$root_phase
    ),
    summarize_method(
      "Activated particles: Wald",
      result$pooled_estimate,
      result$pooled_cover,
      result$pooled_length,
      result$pooled_phase
    ),
    summarize_method(
      "Full hazard MLE: Wald",
      result$hazard_estimate,
      result$hazard_wald_cover,
      result$hazard_wald_length,
      result$hazard_wald_phase
    ),
    summarize_method(
      "Full hazard MLE: likelihood ratio",
      result$hazard_estimate,
      result$hazard_lr_cover,
      result$hazard_lr_length,
      result$hazard_lr_phase
    )
  )

  list(
    draws = result,
    summary = summary,
    beta = beta,
    confidence_level = confidence_level,
    B = B,
    boundary_rate = mean(result$hazard_boundary),
    mean_N_per_realization = mean(result$N / K),
    mean_rows_per_realization = mean(result$B_rows / K)
  )
}

# ----------------------------------------------------------------------------
# User interface
# ----------------------------------------------------------------------------

app_css <- paste0(
  "\n",
  ":root { --frog-navy: ", APP_COLORS$navy, "; --frog-teal: ", APP_COLORS$teal,
  "; --frog-border: ", APP_COLORS$border, "; }\n",
  "body { background: #F4F7FA; }\n",
  ".navbar { box-shadow: 0 4px 18px rgba(16,42,67,.10); }\n",
  ".navbar-brand { font-weight: 800; letter-spacing: -.02em; }\n",
  ".app-brand-subtitle { display:block; font-size:.70rem; opacity:.78; font-weight:500; }\n",
  ".hero-panel { background: linear-gradient(135deg,#102A43 0%,#087F8C 100%); color:white;",
  " border-radius:18px; padding:28px 30px; margin-bottom:18px; box-shadow:0 12px 30px rgba(16,42,67,.18); }\n",
  ".hero-panel h2 { font-weight:800; letter-spacing:-.03em; margin-top:0; }\n",
  ".hero-panel p { max-width:1000px; font-size:1.03rem; opacity:.94; margin-bottom:0; }\n",
  ".section-kicker { text-transform:uppercase; letter-spacing:.10em; font-size:.74rem; font-weight:800; color:#087F8C; }\n",
  ".explanation { background:#F8FBFD; border-left:4px solid #087F8C; border-radius:10px; padding:14px 16px; color:#334E68; }\n",
  ".warning-note { background:#FFF8E8; border-left:4px solid #D97706; border-radius:10px; padding:14px 16px; color:#7C4A03; }\n",
  ".method-note { background:#EEF5FF; border-left:4px solid #1565C0; border-radius:10px; padding:14px 16px; color:#243B53; }\n",
  ".phase-badge { display:inline-block; border-radius:999px; padding:5px 10px; font-weight:750; font-size:.80rem; white-space:nowrap; }\n",
  ".phase-s { background:#DCFCE7; color:#166534; border:1px solid #86EFAC; }\n",
  ".phase-e { background:#FEE2E2; color:#991B1B; border:1px solid #FCA5A5; }\n",
  ".phase-u { background:#FEF3C7; color:#92400E; border:1px solid #FCD34D; }\n",
  ".phase-na { background:#E5E7EB; color:#374151; border:1px solid #CBD5E1; }\n",
  ".phase-color-key { display:grid; grid-template-columns:repeat(3,minmax(0,1fr)); gap:10px; padding:12px 14px; margin:4px 12px 2px; background:#F8FBFD; border:1px solid #D9E2EC; border-radius:12px; }\n",
  ".phase-key-item { display:flex; align-items:flex-start; gap:9px; color:#334E68; font-size:.82rem; line-height:1.25; }\n",
  ".phase-key-swatch { width:18px; height:18px; min-width:18px; margin-top:2px; border-radius:4px; border:1px solid rgba(0,0,0,.14); }\n",
  ".swatch-s { background:#16865B; } .swatch-u { background:#D97706; } .swatch-e { background:#C2413B; }\n",
  ".joint-ratio-box { border:2px solid #D97706; background:#FFF8E8; border-radius:12px; padding:12px 14px; margin:8px 0 12px; color:#6B3F00; }\n",
  ".joint-ratio-value { display:block; font-size:1.35rem; font-weight:900; color:#9A4E00; margin:.15rem 0 .3rem; overflow-wrap:anywhere; }\n",
  ".validation-pass { color:#166534; background:#DCFCE7; border-radius:999px; padding:3px 9px; font-weight:800; }\n",
  ".validation-fail { color:#991B1B; background:#FEE2E2; border-radius:999px; padding:3px 9px; font-weight:800; }\n",
  ".metric-glyph { font-size:2.15rem; font-weight:800; line-height:1; opacity:.90; }\n",
  ".formula-box { background:white; border:1px solid #D9E2EC; border-radius:14px; padding:16px 18px; height:100%; }\n",
  ".formula-box h5 { font-weight:800; color:#102A43; }\n",
  ".math-display { overflow-x:auto; overflow-y:hidden; padding:.20rem 0 .35rem; }\n",
  ".math-inline { white-space:nowrap; }\n",
  "mjx-container[jax='CHTML'][display='true'] { margin:.65rem 0 !important; max-width:100%; overflow-x:auto; overflow-y:hidden; }\n",
  ".card { border:1px solid #E3EAF1; box-shadow:0 6px 18px rgba(16,42,67,.06); border-radius:14px; }\n",
  ".card-header { background:white; font-weight:800; color:#102A43; border-bottom:1px solid #E8EEF4; }\n",
  ".control-caption { color:#627D98; font-size:.80rem; margin-top:-4px; margin-bottom:10px; }\n",
  ".sidebar .form-label, .sidebar label { font-weight:700; color:#243B53; }\n",
  ".btn-primary { font-weight:800; box-shadow:0 5px 12px rgba(8,127,140,.20); }\n",
  ".dataTables_wrapper { font-size:.91rem; }\n",
  ".small-caps { font-size:.72rem; font-weight:800; text-transform:uppercase; letter-spacing:.08em; color:#627D98; }\n",
  ".result-callout { border:1px solid #D9E2EC; border-radius:12px; padding:14px; background:#FFFFFF; }\n",
  ".result-callout .value { font-size:1.35rem; font-weight:850; color:#102A43; }\n",
  ".nav-link { font-weight:650; }\n",
  ".footer-note { color:#627D98; font-size:.82rem; padding:12px 0 24px; text-align:center; }\n",
  "@media (max-width: 768px) { .hero-panel { padding:20px; } .hero-panel h2 { font-size:1.55rem; } .phase-color-key { grid-template-columns:1fr; } }\n"
)

shared_sidebar <- sidebar(
  width = 340,
  open = "desktop",
  title = "Experiment controls",
  radioButtons(
    "data_mode",
    "Data source",
    choices = c(
      "Exact finite-horizon simulation" = "simulate",
      "Analyze an age-risk CSV" = "upload"
    ),
    selected = "simulate"
  ),
  conditionalPanel(
    condition = "input.data_mode == 'simulate'",
    numericInput(
      "beta_value",
      HTML("True parameter &beta;"),
      value = 0.35,
      min = 0.001,
      step = 0.01
    ),
    numericInput(
      "seed",
      "Random-number seed",
      value = 20260621,
      min = 1,
      step = 1
    ),
    numericInput(
      "detail_n",
      "Realizations retained for detailed dynamics",
      value = 12,
      min = 1,
      max = 30,
      step = 1
    )
  ),
  conditionalPanel(
    condition = "input.data_mode == 'upload'",
    fileInput(
      "risk_file",
      "Age-risk table (.csv)",
      accept = c(".csv", "text/csv")
    ),
    numericInput(
      "S_root_upload",
      "Origin-particle jumps among K realizations",
      value = 50,
      min = 0,
      step = 1
    ),
    numericInput(
      "beta_reference_upload",
      HTML("Reference &beta; for overlays (optional)"),
      value = 0.35,
      min = 0.001,
      step = 0.01
    ),
    downloadButton("download_template", "Download CSV template", class = "btn-outline-secondary btn-sm")
  ),
  tags$hr(),
  numericInput("K", "Independent realizations K", value = 100, min = 1, max = 5000, step = 1),
  numericInput("T", "Calendar horizon T", value = 12, min = 1, max = 150, step = 1),
  sliderInput(
    "confidence_level",
    "Confidence level",
    min = 0.80,
    max = 0.999,
    value = 0.95,
    step = 0.005,
    ticks = FALSE
  ),
  tags$div(class = "control-caption", "The origin interval is exact; pooled and likelihood intervals use fixed-horizon asymptotic calibration."),
  numericInput(
    "beta_lower",
    HTML("Lower constrained bound for &beta;"),
    value = 0.005,
    min = 1e-6,
    step = 0.005
  ),
  numericInput(
    "beta_upper",
    HTML("Upper constrained bound for &beta;"),
    value = 20,
    min = 0.01,
    step = 0.5
  ),
  uiOutput("workload_text"),
  uiOutput("joint_growth_text"),
  actionButton(
    "run_analysis",
    "Run experiment and inference",
    class = "btn-primary w-100",
    icon = icon("play")
  ),
  tags$hr(),
  tags$div(
    class = "small text-muted",
    HTML("The phase threshold is fixed at <b><i>&beta;</i><sub>c</sub> = 1/2</b>.")
  )
)

ui <- tagList(
  tags$head(
    tags$script(HTML(
      "window.MathJax = {tex: {inlineMath: [['\\\\(', '\\\\)']], displayMath: [['\\\\[', '\\\\]']]}, options: {skipHtmlTags: ['script','noscript','style','textarea','pre','code']}, chtml: {scale: 1.02}};"
    )),
    tags$script(
      id = "MathJax-script",
      async = NA,
      src = "https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js"
    ),
    tags$script(HTML(
      "(function(){
         function typesetMath(){
           if(window.MathJax && window.MathJax.typesetPromise){
             window.MathJax.typesetPromise().catch(function(err){console.warn('MathJax typesetting warning:', err);});
           }
         }
         document.addEventListener('DOMContentLoaded', function(){
           setTimeout(typesetMath, 150);
           document.addEventListener('shiny:connected', function(){setTimeout(typesetMath, 80);});
           document.addEventListener('shiny:value', function(){setTimeout(typesetMath, 80);});
           document.body.addEventListener('shown.bs.tab', function(){setTimeout(typesetMath, 80);});
         });
         window.frogTypesetMath = typesetMath;
       })();"
    ))
  ),
  page_navbar(
    title = tags$div(
      "Frog-Model Inference Explorer",
      tags$span(class = "app-brand-subtitle", "Finite-horizon likelihood and phase classification")
    ),
    window_title = "Frog-Model Inference Explorer",
    id = "main_navigation",
    fillable = FALSE,
    sidebar = shared_sidebar,
    theme = bs_theme(
      version = 5,
      bootswatch = "flatly",
      primary = APP_COLORS$teal,
      secondary = APP_COLORS$slate,
      success = APP_COLORS$green,
      danger = APP_COLORS$red,
      warning = APP_COLORS$amber,
      info = APP_COLORS$blue
    ),
    header = tagList(tags$style(HTML(app_css))),

    nav_panel(
      "Overview",
      icon = icon("house"),
      div(
        class = "hero-panel",
        tags$div(class = "section-kicker", style = "color:#B9F2EE;", "Interactive companion"),
        tags$h2("Inference from an endogenously generated observation set"),
        tags$p(
          "Simulate or import a finite-horizon labeled experiment, compute all three estimators, compare confidence regions, inspect the exact likelihood, and translate uncertainty about β into an explicit three-way phase decision."
        )
      ),
      layout_columns(
        col_widths = c(3, 3, 3, 3),
        uiOutput("value_beta"),
        uiOutput("value_design"),
        uiOutput("value_rows"),
        uiOutput("value_mle")
      ),
      uiOutput("visible_rows_definition"),
      card(
        card_header("Joint-growth asymptotic diagnostic"),
        card_body(uiOutput("joint_growth_overview"))
      ),
      layout_columns(
        col_widths = c(8, 4),
        card(
          full_screen = TRUE,
          card_header("Confidence regions and the phase threshold"),
          phase_color_key(),
          plotlyOutput("overview_interval_plot", height = "430px"),
          card_footer(
            tags$div(
              class = "small text-muted",
              HTML("Intervals are shown on the transformed display scale <b>log(1 + <i>&beta;</i>)</b>, so that zero, moderate estimates, and long upper tails remain visible together. Color code: <b style='color:#16865B;'>green</b> = survival report (S), <b style='color:#D97706;'>amber/orange</b> = undetermined (U), and <b style='color:#C2413B;'>red</b> = extinction report (E).")
            )
          )
        ),
        card(
          card_header("Interpretation of the current experiment"),
          card_body(uiOutput("overview_narrative"))
        )
      ),
      layout_columns(
        col_widths = c(7, 5),
        card(
          full_screen = TRUE,
          card_header("First-order variance ordering"),
          plotlyOutput("variance_order_plot", height = "380px"),
          card_footer(
            tags$div(
              class = "small text-muted",
              "The displayed coefficients use the chosen reference β and the observed risk table. The sample-level plug-in ordering mirrors the nested information structure of the article."
            )
          )
        ),
        card(
          card_header("What each procedure retains"),
          card_body(
            tags$div(class = "explanation", tags$b("Origin benchmark."), " One guaranteed age-zero transition per realization; exact binomial confidence coverage."),
            tags$br(),
            tags$div(class = "method-note", tags$b("Activated-particle pooling."), " Every visible age-zero transition; the random sample size is endogenous and is not treated as conditionally binomial."),
            tags$br(),
            tags$div(class = "warning-note", tags$b("Full hazard likelihood."), " Every visible age-specific death or jump plus terminal right censoring; regular Wald calibration is not asserted at a constrained boundary fit.")
          )
        )
      )
    ),

    nav_panel(
      "Estimator comparison",
      icon = icon("scale-balanced"),
      layout_columns(
        col_widths = c(7, 5),
        card(
          full_screen = TRUE,
          card_header("Estimator and confidence-set comparison"),
          phase_color_key(),
          plotlyOutput("interval_plot", height = "470px"),
          card_footer(
            tags$div(
              class = "small text-muted",
              HTML("Color code relative to <i>&beta;</i><sub>c</sub> = 0.5: <b style='color:#16865B;'>green</b> = the confidence region lies entirely below <i>&beta;</i><sub>c</sub>, so the app reports the positive-survival regime (S); <b style='color:#D97706;'>amber/orange</b> = the region contains or crosses <i>&beta;</i><sub>c</sub>, so the report is undetermined (U); <b style='color:#C2413B;'>red</b> = the region lies entirely above <i>&beta;</i><sub>c</sub>, so the app reports the almost-sure-extinction regime (E).")
            )
          )
        ),
        card(
          card_header("Observed counts and reductions"),
          card_body(uiOutput("count_summary"))
        )
      ),
      card(
        full_screen = TRUE,
        card_header("Complete inferential summary"),
        DTOutput("estimator_table")
      ),
      layout_columns(
        col_widths = c(4, 4, 4),
        card(
          card_header("Origin-particle estimator"),
          card_body(uiOutput("origin_explanation"))
        ),
        card(
          card_header("Activated-particle estimator"),
          card_body(uiOutput("pooled_explanation"))
        ),
        card(
          card_header("Full-hazard MLE"),
          card_body(uiOutput("hazard_explanation"))
        )
      )
    ),

    nav_panel(
      "Likelihood & hazards",
      icon = icon("chart-line"),
      layout_columns(
        col_widths = c(7, 5),
        card(
          full_screen = TRUE,
          card_header("Likelihood-ratio profile"),
          plotlyOutput("likelihood_plot", height = "430px"),
          card_footer(uiOutput("likelihood_plot_explanation"))
        ),
        card(
          full_screen = TRUE,
          card_header("Monotone score and numerical maximizer"),
          plotlyOutput("score_plot", height = "430px"),
          card_footer(uiOutput("score_plot_explanation"))
        )
      ),
      layout_columns(
        col_widths = c(7, 5),
        card(
          full_screen = TRUE,
          card_header("Observed age-specific hazards"),
          plotlyOutput("hazard_plot", height = "430px"),
          card_footer(
            tags$div(
              class = "small text-muted",
              tagList(
                "Empirical death fractions are descriptive, and point size reflects the number at risk. The fitted curve is ",
                math_inline("h_\\beta(j)=\\beta/(\\beta+j+1)"),
                "."
              )
            )
          )
        ),
        card(
          full_screen = TRUE,
          card_header("Observed information by age"),
          plotlyOutput("information_plot", height = "430px"),
          card_footer(
            tags$div(
              class = "small text-muted",
              tagList(
                "Each bar equals ",
                math_inline("R_j h(j)\\{1-h(j)\\}"),
                " at the fitted parameter. Positive-age rows quantify the information discarded by the first two reductions."
              )
            )
          )
        )
      ),
      card(
        full_screen = TRUE,
        card_header("Age-compressed sufficient table"),
        DTOutput("risk_table")
      )
    ),

    nav_panel(
      "Process dynamics",
      icon = icon("route"),
      uiOutput("dynamics_availability"),
      conditionalPanel(
        condition = "input.data_mode == 'simulate'",
        layout_columns(
          col_widths = c(3, 9),
          card(
            card_header("Detailed realization"),
            card_body(
              uiOutput("realization_selector"),
              uiOutput("realization_metrics")
            )
          ),
          card(
            full_screen = TRUE,
            card_header("Alive particles and cumulative activation"),
            plotlyOutput("dynamics_plot", height = "420px"),
            card_footer(
              tags$div(
                class = "small text-muted",
                "Newly reached home sites are activated simultaneously after all transitions in a calendar slice have been processed; they first face a transition in the following interval."
              )
            )
          )
        ),
        layout_columns(
          col_widths = c(6, 6),
          card(
            full_screen = TRUE,
            card_header("Calendar-time event flow"),
            plotlyOutput("event_flow_plot", height = "390px")
          ),
          card(
            full_screen = TRUE,
            card_header("Visible transition log"),
            DTOutput("event_table")
          )
        )
      )
    ),

    nav_panel(
      "Repeated-sampling explorer",
      icon = icon("repeat"),
      div(
        class = "warning-note mb-3",
        tags$b("Purpose."),
        " This optional module is a numerical operating-characteristic explorer. It is not needed for the proofs and does not replace the exact or asymptotic guarantees stated in the article."
      ),
      layout_columns(
        col_widths = c(3, 9),
        card(
          card_header("Study controls"),
          card_body(
            numericInput("mc_reps", "Number of repeated data sets", value = 100, min = 20, max = 500, step = 10),
            numericInput("mc_seed", "Repeated-sampling seed", value = 87321, min = 1, step = 1),
            uiOutput("mc_workload"),
            actionButton("run_mc", "Run repeated-sampling study", class = "btn-primary w-100", icon = icon("flask"))
          )
        ),
        card(
          full_screen = TRUE,
          card_header("Sampling distributions of the three point estimators"),
          plotlyOutput("mc_distribution_plot", height = "440px")
        )
      ),
      layout_columns(
        col_widths = c(7, 5),
        card(
          full_screen = TRUE,
          card_header("Empirical coverage against the requested level"),
          plotlyOutput("mc_coverage_plot", height = "390px")
        ),
        card(
          card_header("Study diagnostics"),
          card_body(uiOutput("mc_diagnostics"))
        )
      ),
      card(
        full_screen = TRUE,
        card_header("Repeated-sampling numerical summary"),
        DTOutput("mc_summary_table")
      )
    ),

    nav_panel(
      "Data & validation",
      icon = icon("check-double"),
      layout_columns(
        col_widths = c(7, 5),
        card(
          full_screen = TRUE,
          card_header("Automated internal invariants"),
          DTOutput("validation_table")
        ),
        card(
          card_header("Import format and reproducibility"),
          card_body(
            tags$div(
              class = "explanation",
              tags$b("CSV format."),
              " Supply three columns named age, R, and d. Ages may be omitted when both counts are zero; the app fills the complete range 0,…,T−1."
            ),
            tags$br(),
            tags$div(
              class = "method-note",
              tags$b("Origin information."),
              " The age table does not identify which age-zero rows came from the K origin particles, so uploaded analyses must also provide the number of origin jumps."
            ),
            tags$br(),
            tags$div(
              class = "warning-note",
              tags$b("Reproducibility."),
              " Exact simulations use the displayed seed. Re-running with the same inputs and software environment reproduces tables, estimates, and figures."
            ),
            tags$hr(),
            downloadButton("download_risk", "Download current risk table", class = "btn-outline-primary w-100"),
            tags$br(), tags$br(),
            downloadButton("download_summary", "Download inferential summary", class = "btn-outline-primary w-100")
          )
        )
      ),
      conditionalPanel(
        condition = "input.data_mode == 'simulate'",
        card(
          full_screen = TRUE,
          card_header("Per-realization finite-horizon summary"),
          DTOutput("per_realization_table")
        )
      )
    ),

    nav_panel(
      "Methods guide",
      icon = icon("book-open"),
      div(
        class = "hero-panel",
        tags$div(class = "section-kicker", style = "color:#B9F2EE;", "Model and computation"),
        tags$h2("A compact guide to every quantity shown in the app"),
        tags$p("The app is designed as a self-contained computational companion: every graph is linked to an observable count, an estimator, or a confidence-set statement from the article.")
      ),
      layout_columns(
        col_widths = c(4, 4, 4),
        div(
          class = "formula-box",
          tags$h5("Lifetime law and hazard"),
          math_display("P_\\beta(\\Xi\\ge j)=\\frac{\\Gamma(\\beta+1)\\Gamma(j+1)}{\\Gamma(\\beta+j+1)},\\qquad h_\\beta(j)=\\frac{\\beta}{\\beta+j+1}."),
          tags$p("The simulator draws the next death-or-jump outcome directly from this age-dependent hazard; latent survival marks and complete lifetimes need not be generated.")
        ),
        div(
          class = "formula-box",
          tags$h5("Origin benchmark"),
          math_display("\\widehat\\beta^{\\mathrm{root}}=\\frac{K-S^{\\mathrm{root}}}{S^{\\mathrm{root}}},\\qquad \\widetilde q^{\\mathrm{root}}=\\frac{S^{\\mathrm{root}}+1/2}{K+1}."),
          tags$p(tagList(
            "The confidence set transforms an equal-tailed Clopper–Pearson interval through the decreasing map ",
            math_inline("\\beta=q^{-1}-1"),
            ", reversing endpoints."
          ))
        ),
        div(
          class = "formula-box",
          tags$h5("Activated-particle pooling"),
          math_display("\\widehat\\beta^{\\mathrm{act}}=\\frac{N_{K,T}}{M_{K,T}}-1,\\qquad \\widehat{\\mathrm{se}}(\\widetilde\\theta)=\\{N_{K,T}\\widetilde q(1-\\widetilde q)\\}^{-1/2}."),
          tags$p(tagList(
            "The apparent Bernoulli form is justified by predictable truncation across endogenous activations—not by assuming ",
            math_inline("M_{K,T}\\mid N_{K,T}"),
            " is binomial."
          ))
        )
      ),
      tags$br(),
      layout_columns(
        col_widths = c(6, 6),
        div(
          class = "formula-box",
          tags$h5("Age-compressed likelihood"),
          math_display("\\ell^*(\\theta)=\\sum_{j=0}^{T-1}\\left[d_j z_j(\\theta)-R_j\\,\\operatorname{softplus}\\{z_j(\\theta)\\}\\right],\\quad z_j(\\theta)=\\theta-\\log(j+1)."),
          math_display("U(\\theta)=D^+-\\sum_jR_jh_\\theta(j),\\qquad \\mathcal J(\\theta)=\\sum_jR_jh_\\theta(j)\\{1-h_\\theta(j)\\}."),
          tags$p(tagList(
            "The app evaluates ", math_inline("h_\\theta(j)"),
            " through plogis() and softplus through a stable max-plus-log1p identity. The score is continuous and strictly decreasing."
          ))
        ),
        div(
          class = "formula-box",
          tags$h5("Confidence-set phase rule"),
          math_display("\\delta(C)=\\begin{cases}\\mathsf S,&\\sup C<1/2,\\\\\\mathsf E,&\\inf C>1/2,\\\\\\mathsf U,&1/2\\in C.\\end{cases}"),
          tags$p(tagList(
            "The origin rule has finite-sample phase-error control. Pooled and likelihood rules are asymptotically calibrated for fixed ",
            math_inline("T"), " as ", math_inline("K\\to\\infty"),
            ". The dynamical theorem intentionally makes no assertion at ", math_inline("\\beta=1/2"), "."
          ))
        )
      ),
      tags$br(),
      card(
        card_header("Growing-horizon asymptotic condition"),
        card_body(
          tags$div(
            class = "method-note",
            tags$b("Why the app reports a K-T diagnostic."),
            " For the triangular-array theorem in which both the number of realizations and the time horizon increase, the sufficient condition is that the deterministic ratio below converges to zero. Larger T is therefore not automatically better when K is kept fixed or grows too slowly."
          ),
          tags$br(),
          math_display("\\rho(K,T)=\\frac{T^4\\log(T^2+2)}{K}."),
          tags$p(tagList(
            "The fixed-", math_inline("T"), " theory does not require this condition. The displayed number is instead a growing-horizon diagnostic: smaller values indicate that the current pair ", math_inline("(K,T)"), " is more compatible with the joint-growth regime used in the proof."
          ))
        )
      ),
      tags$br(),
      card(
        card_header("Reading the figures"),
        card_body(
          tags$ol(
            tags$li(tags$b("Interval comparison:"), " determines whether each procedure makes a survival, extinction, or undetermined report."),
            tags$li(tags$b("Likelihood profile:"), " shows all β values whose log-likelihood drop lies within the χ² threshold."),
            tags$li(tags$b("Hazard plot:"), " compares age-specific empirical death fractions with the fitted and reference hazards."),
            tags$li(tags$b("Information plot:"), " attributes observed curvature to individual ages and reveals the contribution of right-censored positive-age rows."),
            tags$li(tags$b("Dynamics plots:"), " make the timing convention and endogenous activation mechanism visible in a selected realization."),
            tags$li(tags$b("Repeated-sampling plots:"), " illustrate numerical coverage, boundary rates, and estimator dispersion under chosen finite K and T.")
          )
        )
      )
    ),

    footer = tags$div(
      class = "footer-note",
      "Finite-Horizon Frog-Model Inference Explorer · single-file Shiny implementation"
    )
  )
)

# ----------------------------------------------------------------------------
# Server
# ----------------------------------------------------------------------------

server <- function(input, output, session) {

  output$workload_text <- renderUI({
    K <- max(1, as.integer(input$K %||% 1))
    T <- max(1, as.integer(input$T %||% 1))
    worst <- K * T^2
    class_name <- if (worst > 2e6) "warning-note" else "control-caption"
    tags$div(
      class = class_name,
      HTML(paste0(
        "Deterministic work bound: at most <b>", fmt_int(worst),
        "</b> visible-row updates (K T<sup>2</sup>)."
      ))
    )
  })

  output$joint_growth_text <- renderUI({
    K <- max(1, as.integer(input$K %||% 1))
    T <- max(1, as.integer(input$T %||% 1))
    rho <- joint_growth_ratio(K, T)
    rho_next <- joint_growth_ratio(K, T + 1)
    tags$div(
      class = "joint-ratio-box",
      tags$div(class = "small-caps", "K-T joint-growth diagnostic"),
      tags$span(class = "joint-ratio-value", paste0("rho(K,T) = ", fmt_num(rho, 6))),
      tags$div(HTML(paste0(
        "Current inputs: <b>K = ", fmt_int(K), "</b>, <b>T = ", fmt_int(T), "</b>.<br>",
        "The joint-growth theorem requires <b>rho(K<sub>n</sub>,T<sub>n</sub>) &rarr; 0</b>, where ",
        "rho(K,T) = T<sup>4</sup> log(T<sup>2</sup> + 2) / K.<br>",
        "At the same K, changing the horizon to T + 1 would give rho = <b>", fmt_num(rho_next, 6), "</b>. ",
        "Thus, increasing T alone can move the design farther from the proved joint-growth regime."
      )))
    )
  })

  output$joint_growth_overview <- renderUI({
    K <- max(1, as.integer(input$K %||% 1))
    T <- max(1, as.integer(input$T %||% 1))
    rho <- joint_growth_ratio(K, T)
    tagList(
      layout_columns(
        col_widths = c(4, 8),
        tags$div(
          class = "result-callout",
          tags$div(class = "small-caps", "Joint-growth ratio (not required for fixed T)"),
          tags$div(class = "value", fmt_num(rho, 7)),
          tags$div(class = "small text-muted", HTML("<i>&rho;</i>(<i>K</i>,<i>T</i>) = <i>T</i><sup>4</sup> log(<i>T</i><sup>2</sup> + 2) / <i>K</i>"))
        ),
        tags$div(
          class = "warning-note",
          tags$b("Do not interpret a longer horizon as automatically better."),
          tags$p(
            "For fixed T, consistency follows by increasing the number K of independent realizations, with no K-T rate restriction. When T also increases, the theorem instead assumes a sequence for which this displayed ratio converges to zero."
          ),
          tags$p(
            "The number is a conservative theoretical diagnostic rather than a finite-sample pass/fail cutoff. At fixed K, increasing T raises the numerator approximately at fourth-order speed, apart from the logarithmic factor."
          )
        )
      )
    )
  })

  output$mc_workload <- renderUI({
    B <- max(1, as.integer(input$mc_reps %||% 1))
    K <- max(1, as.integer(input$K %||% 1))
    T <- max(1, as.integer(input$T %||% 1))
    total <- B * K * T^2
    tags$div(
      class = if (total > 2e7) "warning-note" else "control-caption",
      HTML(paste0(
        "Worst-case repeated-sampling work: <b>", fmt_int(total),
        "</b> visible-row updates."
      ))
    )
  })

  output$download_template <- downloadHandler(
    filename = function() paste0("frog_model_risk_table_T", input$T, ".csv"),
    content = function(file) {
      T <- max(1L, as.integer(input$T))
      utils::write.csv(
        data.frame(age = 0:(T - 1L), R = c(input$K, rep(0, T - 1L)), d = 0),
        file,
        row.names = FALSE
      )
    }
  )

  analysis_data <- eventReactive(
    input$run_analysis,
    {
      tryCatch({
        K <- as.integer(input$K)
        T <- as.integer(input$T)
        if (!is.finite(K) || K < 1L) stop("K must be a positive integer.")
        if (!is.finite(T) || T < 1L) stop("T must be a positive integer.")
        if (input$beta_lower <= 0 || input$beta_upper <= input$beta_lower) {
          stop("The constrained beta bounds must satisfy 0 < lower < upper.")
        }

        if (identical(input$data_mode, "simulate")) {
          beta <- as.numeric(input$beta_value)
          if (!is.finite(beta) || beta <= 0) stop("The simulation parameter beta must be positive.")
          detail_n <- min(K, max(1L, as.integer(input$detail_n)))

          withProgress(message = "Exact finite-horizon simulation", value = 0, {
            simulate_dataset(
              K = K,
              T = T,
              beta = beta,
              seed = input$seed,
              detail_n = detail_n,
              progress = function(amount, detail) {
                incProgress(amount, detail = detail)
              }
            )
          })
        } else {
          req(input$risk_file)
          tab <- read_risk_table(input$risk_file$datapath, T)
          S_root <- as.integer(input$S_root_upload)
          if (!is.finite(S_root) || S_root < 0L || S_root > K) {
            stop("Origin-particle jumps must be an integer between 0 and K.")
          }
          if (tab$R[1] < K) {
            stop("A valid K-realization labeled experiment must have R_0 >= K.")
          }
          M_zero <- tab$R[1] - tab$d[1]
          if (S_root > M_zero || (K - S_root) > tab$d[1]) {
            stop("The supplied origin outcomes are incompatible with the age-zero jump/death totals.")
          }
          if (T > 1L && any(tab$R[-1L] > tab$R[-T] - tab$d[-T])) {
            stop("The uploaded table violates R_{j+1} <= R_j - d_j for at least one age.")
          }
          list(
            source = "Uploaded age-risk table",
            K = K,
            T = T,
            beta_reference = as.numeric(input$beta_reference_upload),
            seed = NA_integer_,
            detail_n = 0L,
            R = tab$R,
            d = tab$d,
            S_root = S_root,
            per_realization = NULL,
            dynamics = NULL,
            events = NULL
          )
        }
      }, error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = 10)
        NULL
      })
    },
    ignoreNULL = FALSE
  )

  estimates <- reactive({
    dat <- analysis_data()
    req(dat)
    tryCatch(
      compute_all_estimators(
        dat,
        confidence_level = input$confidence_level,
        beta_bounds = c(input$beta_lower, input$beta_upper)
      ),
      error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = 10)
        NULL
      }
    )
  })

  reference_beta <- reactive({
    dat <- analysis_data()
    est <- estimates()
    req(dat, est)
    beta <- dat$beta_reference
    if (!is.finite(beta) || beta <= 0) beta <- est$hazard$estimate
    beta
  })

  interval_data <- reactive({
    est <- estimates()
    req(est)

    rows <- list(
      data.frame(
        Method = "Origin · exact CP",
        Estimate = est$origin$estimate,
        Lower = est$origin$beta_ci[1],
        Upper = est$origin$beta_ci[2],
        Phase = est$origin$phase,
        Calibration = "Exact finite-sample",
        stringsAsFactors = FALSE
      ),
      data.frame(
        Method = "Activated · Wald",
        Estimate = est$pooled$estimate,
        Lower = est$pooled$beta_ci[1],
        Upper = est$pooled$beta_ci[2],
        Phase = est$pooled$phase,
        Calibration = "Asymptotic fixed-T",
        stringsAsFactors = FALSE
      ),
      data.frame(
        Method = "Full hazard · LR",
        Estimate = est$hazard$estimate,
        Lower = est$hazard$lr_ci[1],
        Upper = est$hazard$lr_ci[2],
        Phase = est$hazard$lr_phase,
        Calibration = "Asymptotic fixed-T",
        stringsAsFactors = FALSE
      )
    )

    if (!any(is.na(est$hazard$wald_ci))) {
      rows <- append(
        rows,
        list(data.frame(
          Method = "Full hazard · Wald",
          Estimate = est$hazard$estimate,
          Lower = est$hazard$wald_ci[1],
          Upper = est$hazard$wald_ci[2],
          Phase = est$hazard$wald_phase,
          Calibration = "Asymptotic fixed-T",
          stringsAsFactors = FALSE
        )),
        after = 2L
      )
    }

    do.call(rbind, rows)
  })

  output$value_beta <- renderUI({
    beta <- reference_beta()
    value_box(
      title = if (analysis_data()$source == "Exact simulation") "True parameter β" else "Reference parameter β",
      value = fmt_num(beta, 5),
      showcase = tags$div(class = "metric-glyph", "β"),
      theme = "bg-gradient-blue-teal"
    )
  })

  output$value_design <- renderUI({
    dat <- analysis_data()
    req(dat)
    value_box(
      title = "Experimental design",
      value = paste0("K = ", fmt_int(dat$K), " · T = ", fmt_int(dat$T)),
      showcase = tags$div(class = "metric-glyph", "K×T"),
      theme = "bg-gradient-indigo-blue"
    )
  })

  output$value_rows <- renderUI({
    dat <- analysis_data()
    est <- estimates()
    req(dat, est)
    value_box(
      title = HTML("Total visible particle–age rows: <i>B</i><sub>K,T</sub>"),
      value = fmt_int(est$hazard$B_plus),
      tags$div(
        style = "font-size:.78rem; line-height:1.25; margin-top:.25rem;",
        HTML(paste0(
          "<i>B</i><sub>K,T</sub> = &sum;<sub>j=0</sub><sup>T&minus;1</sup> <i>R</i><sub>K,T,j</sub>",
          "<br><i>K</i> = ", fmt_int(dat$K),
          " &le; <i>B</i><sub>K,T</sub> &le; <i>K</i><i>T</i><sup>2</sup> = ",
          fmt_int(dat$K * dat$T^2)
        ))
      ),
      showcase = tags$div(class = "metric-glyph", HTML("<i>B</i><sub>K,T</sub>")),
      theme = "bg-gradient-green-teal"
    )
  })

  output$visible_rows_definition <- renderUI({
    dat <- analysis_data()
    est <- estimates()
    req(dat, est)

    tags$div(
      class = "explanation mb-3",
      HTML(paste0(
        "<b>Total visible particle–age rows:</b> ",
        "<i>B</i><sub>K,T</sub> = &sum;<sub>j=0</sub><sup>T&minus;1</sup>",
        " <i>R</i><sub>K,T,j</sub> = <b>", fmt_int(est$hazard$B_plus), "</b>.",
        " &nbsp; The deterministic finite-speed bounds are ",
        "<b><i>K</i> &le; <i>B</i><sub>K,T</sub> &le; <i>K</i><i>T</i><sup>2</sup></b>, ",
        "which here become <b>", fmt_int(dat$K), " &le; ",
        fmt_int(est$hazard$B_plus), " &le; ", fmt_int(dat$K * dat$T^2), "</b>.",
        " The symbol <i>R</i><sub>K,T,j</sub> denotes the age-<i>j</i> risk count; ",
        "their sum is the scalar <i>B</i><sub>K,T</sub> displayed by the app."
      ))
    )
  })

  output$value_mle <- renderUI({
    est <- estimates()
    req(est)
    value_box(
      title = "Full-hazard MLE",
      value = fmt_num(est$hazard$estimate, 5),
      tags$span(style = "font-size:.82rem;", est$hazard$boundary_flag),
      showcase = tags$div(class = "metric-glyph", "ℓ"),
      theme = if (est$hazard$boundary_flag == "Interior") "bg-gradient-purple-blue" else "warning"
    )
  })

  make_interval_plot <- function() {
    df <- interval_data()
    beta_ref <- reference_beta()

    finite_estimates <- df$Estimate[is.finite(df$Estimate)]
    cap_beta <- max(c(2, beta_ref * 1.35, finite_estimates * 1.25, BETA_CRITICAL * 2), na.rm = TRUE)
    cap_beta <- min(max(cap_beta, 1), max(100, cap_beta))

    df$Lower_plot <- pmin(df$Lower, cap_beta)
    df$Upper_plot <- ifelse(is.infinite(df$Upper), cap_beta, pmin(df$Upper, cap_beta))
    df$Estimate_plot <- ifelse(is.infinite(df$Estimate), cap_beta, pmin(df$Estimate, cap_beta))
    df$Lower_t <- log1p(df$Lower_plot)
    df$Upper_t <- log1p(df$Upper_plot)
    df$Estimate_t <- log1p(df$Estimate_plot)
    df$Truncated <- is.infinite(df$Upper) | df$Upper > cap_beta
    df$Interval_text <- mapply(
      function(lower, upper) fmt_interval(c(lower, upper)),
      df$Lower,
      df$Upper,
      USE.NAMES = FALSE
    )
    df$Hover <- paste0(
      "<b>", df$Method, "</b><br>",
      "Point estimate: ", fmt_num(df$Estimate), "<br>",
      "Confidence region: ", df$Interval_text, "<br>",
      "Phase report: ", phase_label(df$Phase), "<br>",
      "Calibration: ", df$Calibration
    )
    method_order <- rev(df$Method)
    df$Method <- factor(df$Method, levels = method_order)

    phase_colors <- c(
      S = APP_COLORS$green,
      E = APP_COLORS$red,
      U = APP_COLORS$amber,
      "NA" = APP_COLORS$slate
    )

    p <- plot_ly()
    for (i in seq_len(nrow(df))) {
      phase_key <- as.character(df$Phase[i])
      if (is.na(phase_key) || !phase_key %in% names(phase_colors)) phase_key <- "NA"
      line_color <- unname(phase_colors[[phase_key]])
      method_i <- as.character(df$Method[i])

      p <- p %>%
        add_segments(
          x = df$Lower_t[i], xend = df$Upper_t[i],
          y = method_i, yend = method_i,
          line = list(color = line_color, width = 10),
          hoverinfo = "text",
          text = df$Hover[i],
          showlegend = FALSE,
          name = method_i
        ) %>%
        add_markers(
          x = df$Estimate_t[i], y = method_i,
          marker = list(
            size = 12,
            color = "white",
            line = list(color = line_color, width = 3)
          ),
          hoverinfo = "text",
          text = df$Hover[i],
          showlegend = FALSE,
          name = method_i
        )

      if (isTRUE(df$Truncated[i])) {
        p <- p %>% add_markers(
          x = df$Upper_t[i], y = method_i,
          marker = list(
            size = 11,
            symbol = "triangle-up",
            color = "white",
            line = list(color = line_color, width = 2)
          ),
          hoverinfo = "text",
          text = paste0(df$Hover[i], "<br><i>Upper endpoint truncated for display.</i>"),
          showlegend = FALSE,
          name = method_i
        )
      }
    }

    tick_beta <- unique(sort(c(0, 0.1, 0.25, 0.5, 1, 2, 5, 10, cap_beta)))
    plotly_clean(
      p %>% layout(
        xaxis = list(
          title = "Parameter <i>&beta;</i> on the log(1 + <i>&beta;</i>) display scale",
          range = c(0, log1p(cap_beta * 1.03)),
          tickvals = log1p(tick_beta),
          ticktext = fmt_num(tick_beta, 3),
          zeroline = FALSE
        ),
        yaxis = list(
          title = "",
          categoryorder = "array",
          categoryarray = method_order,
          automargin = TRUE
        ),
        shapes = list(
          list(
            type = "line", xref = "x", yref = "paper",
            x0 = log1p(BETA_CRITICAL), x1 = log1p(BETA_CRITICAL), y0 = 0, y1 = 1,
            line = list(color = APP_COLORS$red, width = 2, dash = "dash")
          ),
          list(
            type = "line", xref = "x", yref = "paper",
            x0 = log1p(beta_ref), x1 = log1p(beta_ref), y0 = 0, y1 = 1,
            line = list(color = APP_COLORS$navy, width = 2, dash = "dot")
          )
        ),
        annotations = list(
          list(
            x = log1p(BETA_CRITICAL), y = 1.04, xref = "x", yref = "paper",
            text = "Critical threshold: <i>&beta;</i><sub>c</sub> = 0.5",
            showarrow = FALSE,
            font = list(color = APP_COLORS$red, size = 12)
          )
        ),
        hoverlabel = list(align = "left"),
        showlegend = FALSE,
        margin = list(l = 150, r = 25, t = 45, b = 70)
      )
    )
  }

  output$overview_interval_plot <- renderPlotly(make_interval_plot())
  output$interval_plot <- renderPlotly(make_interval_plot())

  output$overview_narrative <- renderUI({
    dat <- analysis_data()
    est <- estimates()
    beta_ref <- reference_beta()
    req(dat, est)

    lr_phase <- est$hazard$lr_phase
    source_text <- if (dat$source == "Exact simulation") {
      paste0(
        "The data were generated exactly from β = ", fmt_num(beta_ref),
        " with ", fmt_int(dat$K), " independent realizations observed through T = ",
        fmt_int(dat$T), "."
      )
    } else {
      paste0(
        "The analysis uses an uploaded sufficient age table for K = ", fmt_int(dat$K),
        " and T = ", fmt_int(dat$T), ". The reference line is not used to fit the model."
      )
    }

    tagList(
      tags$p(source_text),
      tags$p(
        HTML(paste0(
          "The experiment contains <b>", fmt_int(est$pooled$N),
          "</b> visible age-zero rows and <b><i>B</i><sub>K,T</sub> = ", fmt_int(est$hazard$B_plus),
          "</b> total visible particle–age rows."
        ))
      ),
      tags$p(
        HTML(paste0(
          "The full-hazard estimate is <b>", fmt_num(est$hazard$estimate),
          "</b>, with likelihood-ratio region <b>", fmt_interval(est$hazard$lr_ci),
          "</b>."
        ))
      ),
      tags$div(class = "result-callout", tags$div(class = "small-caps", "Likelihood-based phase report"), tags$br(), phase_badge(lr_phase)),
      tags$br(),
      if (est$hazard$boundary_flag != "Interior") {
        tags$div(
          class = "warning-note",
          "The constrained maximum lies on a computational boundary. The app reports the likelihood-ratio region, while the regular interior Wald interval is withheld."
        )
      } else {
        tags$div(
          class = "explanation",
          "The score crosses zero inside the selected parameter range, so the MLE is an interior, unique maximizer of the strictly concave log-likelihood."
        )
      }
    )
  })

  output$variance_order_plot <- renderPlotly({
    dat <- analysis_data()
    est <- estimates()
    beta_eval <- reference_beta()
    req(dat, est)

    mu_hat <- est$pooled$N / dat$K
    info_eval <- est$hazard$information_function(log(beta_eval)) / dat$K

    coeff <- data.frame(
      Method = factor(
        c("Full hazard MLE", "Activated particles", "Origin particle"),
        levels = c("Full hazard MLE", "Activated particles", "Origin particle")
      ),
      Coefficient = c(
        beta_eval^2 / info_eval,
        beta_eval * (beta_eval + 1)^2 / mu_hat,
        beta_eval * (beta_eval + 1)^2
      ),
      Detail = c(
        paste0("Observed information per realization: ", fmt_num(info_eval)),
        paste0("Observed age-zero rows per realization: ", fmt_num(mu_hat)),
        "One guaranteed origin row per realization"
      )
    )

    p <- plot_ly(
      coeff,
      x = ~Method,
      y = ~Coefficient,
      type = "bar",
      text = ~paste0(
        "<b>", Method, "</b><br>",
        "Variance coefficient: ", fmt_num(Coefficient), "<br>", Detail
      ),
      hoverinfo = "text",
      marker = list(color = c(APP_COLORS$purple, APP_COLORS$teal, APP_COLORS$slate))
    ) %>%
      layout(
        xaxis = list(title = ""),
        yaxis = list(title = "Plug-in asymptotic variance coefficient", rangemode = "tozero"),
        margin = list(l = 70, r = 20, t = 20, b = 70)
      )

    plotly_clean(p)
  })

  output$count_summary <- renderUI({
    dat <- analysis_data()
    est <- estimates()
    req(dat, est)
    N <- est$pooled$N
    M <- est$pooled$M
    B <- est$hazard$B_plus
    D <- est$hazard$D_plus

    tagList(
      tags$div(class = "result-callout", tags$div(class = "small-caps", "Origin sample"), tags$div(class = "value", paste0(est$origin$S, " jumps / ", dat$K, " origins"))),
      tags$br(),
      tags$div(class = "result-callout", tags$div(class = "small-caps", "Activated-particle reduction"), tags$div(class = "value", paste0(fmt_int(M), " jumps / ", fmt_int(N), " first transitions"))),
      tags$br(),
      tags$div(
        class = "result-callout",
        tags$div(class = "small-caps", HTML("Full hazard reduction: <i>B</i><sub>K,T</sub>")),
        tags$div(class = "value", paste0(fmt_int(D), " deaths / ", fmt_int(B), " visible particle–age rows")),
        tags$div(class = "small text-muted", HTML(
          "<i>B</i><sub>K,T</sub> = &sum;<sub>j=0</sub><sup>T&minus;1</sup> <i>R</i><sub>K,T,j</sub>, with <i>K</i> &le; <i>B</i><sub>K,T</sub> &le; <i>K</i><i>T</i><sup>2</sup>."
        ))
      ),
      tags$br(),
      tags$div(
        class = "explanation",
        HTML(paste0(
          "Observed enlargement from endogenous activation: <b>N/K = ",
          fmt_num(N / dat$K), "</b> age-zero rows per independent realization."
        ))
      )
    )
  })

  estimator_summary_frame <- reactive({
    est <- estimates()
    req(est)

    data.frame(
      Procedure = c(
        "Origin particle",
        "Activated particles",
        "Full hazard MLE",
        "Full hazard MLE"
      ),
      Point_estimate = c(
        fmt_num(est$origin$estimate),
        fmt_num(est$pooled$estimate),
        fmt_num(est$hazard$estimate),
        fmt_num(est$hazard$estimate)
      ),
      Raw_estimate = c(
        fmt_num(est$origin$raw_beta),
        fmt_num(est$pooled$raw_beta),
        fmt_num(est$hazard$estimate),
        fmt_num(est$hazard$estimate)
      ),
      Confidence_region = c(
        fmt_interval(est$origin$beta_ci),
        fmt_interval(est$pooled$beta_ci),
        fmt_interval(est$hazard$wald_ci),
        fmt_interval(est$hazard$lr_ci)
      ),
      Region_type = c(
        "Exact Clopper–Pearson transform",
        "Log-scale Wald",
        "Observed-information Wald",
        "Likelihood ratio"
      ),
      Standard_error = c(
        paste0("Asymptotic β-scale: ", fmt_num(est$origin$asymptotic_se_beta)),
        paste0("θ-scale: ", fmt_num(est$pooled$se_theta)),
        if (est$hazard$boundary_flag == "Interior") paste0("θ-scale: ", fmt_num(est$hazard$se_theta)) else "Withheld at boundary",
        "Profile-based region"
      ),
      Phase_decision = c(
        phase_badge_html(est$origin$phase),
        phase_badge_html(est$pooled$phase),
        phase_badge_html(est$hazard$wald_phase),
        phase_badge_html(est$hazard$lr_phase)
      ),
      Calibration = c(
        "Finite-sample exact coverage ≥ nominal",
        "Fixed-T asymptotic",
        "Fixed-T asymptotic; interior fits",
        "Fixed-T asymptotic"
      ),
      stringsAsFactors = FALSE
    )
  })

  output$estimator_table <- renderDT({
    datatable(
      estimator_summary_frame(),
      escape = FALSE,
      rownames = FALSE,
      filter = "none",
      options = list(
        pageLength = 8,
        dom = "tip",
        scrollX = TRUE,
        columnDefs = list(list(className = "dt-left", targets = "_all"))
      ),
      class = "stripe hover compact"
    )
  })

  output$origin_explanation <- renderUI({
    est <- estimates()$origin
    tagList(
      tags$p(HTML(paste0(
        "The app observed <b>S = ", est$S, "</b> origin jumps among <b>K = ", est$K,
        "</b> independent realizations."
      ))),
      math_display("\\widetilde q=(S+1/2)/(K+1),\\qquad \\widetilde\\beta=\\widetilde q^{-1}-1."),
      tags$p(HTML(paste0(
        "Smoothed estimate: <b>", fmt_num(est$estimate),
        "</b>. Exact confidence region for beta: <b>", fmt_interval(est$beta_ci), "</b>."
      ))),
      phase_badge(est$phase),
      tags$hr(),
      tags$p(class = "small text-muted", "Exactness comes from one Bernoulli origin transition per independent realization. The add-one-half estimate is used only to avoid infinite numerical summaries at finite-sample boundaries.")
    )
  })

  output$pooled_explanation <- renderUI({
    est <- estimates()$pooled
    tagList(
      tags$p(HTML(paste0(
        "The pooled reduction contains <b>N = ", fmt_int(est$N),
        "</b> visible first transitions, of which <b>M = ", fmt_int(est$M), "</b> are jumps."
      ))),
      math_display("\\widetilde q=(M+1/2)/(N+1),\\qquad \\widetilde\\theta=\\log\\{(1-\\widetilde q)/\\widetilde q\\}."),
      tags$p(HTML(paste0(
        "Regularized estimate: <b>", fmt_num(est$estimate),
        "</b>. Wald confidence region for beta: <b>", fmt_interval(est$beta_ci), "</b>."
      ))),
      phase_badge(est$phase),
      tags$hr(),
      tags$p(class = "small text-muted", "The interval is justified through the predictable-truncation martingale identity and independent-realization asymptotics. M conditional on N is generally not binomial.")
    )
  })

  output$hazard_explanation <- renderUI({
    est <- estimates()$hazard
    tagList(
      tags$p(HTML(paste0(
        "The sufficient table contains <b><i>B</i><sub>K,T</sub> = ", fmt_int(est$B_plus),
        "</b> visible particle–age rows and <b><i>D</i><sup>+</sup><sub>K,T</sub> = ", fmt_int(est$D_plus), "</b> observed deaths."
      ))),
      tags$p(HTML(paste0(
        "MLE: <b>", fmt_num(est$estimate), "</b>; observed information: <b>",
        fmt_num(est$information), "</b>."
      ))),
      tags$p(HTML(paste0(
        "Likelihood-ratio confidence region for beta: <b>", fmt_interval(est$lr_ci), "</b>."
      ))),
      phase_badge(est$lr_phase),
      tags$hr(),
      tags$p(class = "small text-muted", paste0(
        "Numerical status: ", est$boundary_flag,
        ". The stable root-solver and direct one-dimensional maximizer differ by ",
        formatC(est$optimization_difference, format = "e", digits = 2), "."
      ))
    )
  })

  output$likelihood_plot <- renderPlotly({
    est <- estimates()$hazard
    beta_ref <- reference_beta()
    req(est)

    theta_grid <- seq(est$theta_bounds[1], est$theta_bounds[2], length.out = 500L)
    beta_grid <- exp(theta_grid)
    deviance <- 2 * (est$loglik_hat - vapply(theta_grid, est$loglik, numeric(1)))
    y_max <- max(est$chi_crit * 3.5, min(max(deviance, na.rm = TRUE), est$chi_crit * 12))

    df <- data.frame(beta = beta_grid, deviance = deviance)
    p <- plot_ly(
      df,
      x = ~beta,
      y = ~deviance,
      type = "scatter",
      mode = "lines",
      line = list(color = APP_COLORS$purple, width = 3),
      text = ~paste0("<i>&beta;</i> = ", fmt_num(beta), "<br>2 log-likelihood drop = ", fmt_num(deviance)),
      hoverinfo = "text",
      name = "Likelihood-ratio statistic"
    ) %>%
      add_segments(
        x = est$lr_ci[1], xend = est$lr_ci[2],
        y = est$chi_crit, yend = est$chi_crit,
        line = list(color = APP_COLORS$amber, width = 9),
        opacity = 0.35,
        hoverinfo = "text",
        text = paste0("Accepted LR region: ", fmt_interval(est$lr_ci)),
        name = "LR confidence region"
      ) %>%
      add_lines(
        x = c(min(beta_grid), max(beta_grid)),
        y = c(est$chi_crit, est$chi_crit),
        line = list(color = APP_COLORS$amber, dash = "dash", width = 2),
        hoverinfo = "text",
        text = paste0("χ² threshold = ", fmt_num(est$chi_crit)),
        name = "χ² threshold"
      ) %>%
      add_segments(
        x = est$estimate, xend = est$estimate,
        y = 0, yend = y_max,
        line = list(color = APP_COLORS$navy, dash = "dot", width = 2),
        hoverinfo = "text",
        text = paste0("MLE <i>&beta;</i> = ", fmt_num(est$estimate)),
        name = "MLE"
      ) %>%
      add_segments(
        x = BETA_CRITICAL, xend = BETA_CRITICAL,
        y = 0, yend = y_max,
        line = list(color = APP_COLORS$red, dash = "dash", width = 2),
        hoverinfo = "text",
        text = "Critical threshold <i>&beta;</i><sub>c</sub> = 0.5",
        name = "Phase threshold"
      ) %>%
      add_segments(
        x = beta_ref, xend = beta_ref,
        y = 0, yend = y_max,
        line = list(color = APP_COLORS$teal, dash = "dashdot", width = 2),
        hoverinfo = "text",
        text = paste0("Reference <i>&beta;</i> = ", fmt_num(beta_ref)),
        name = "Reference β"
      ) %>%
      layout(
        xaxis = list(title = "Parameter <i>&beta;</i> (logarithmic axis)", type = "log"),
        yaxis = list(title = "Likelihood-ratio statistic  2{&ell;(&beta;&#770;<sub>K,T</sub>) &minus; &ell;(&beta;)}", range = c(0, y_max)),
        legend = list(orientation = "h", y = -0.22),
        margin = list(l = 75, r = 20, t = 20, b = 85)
      )

    plotly_clean(p)
  })

  output$likelihood_plot_explanation <- renderUI({
    est <- estimates()$hazard
    req(est)
    tags$div(
      class = "small text-muted",
      HTML(paste0(
        "The vertical coordinate is <b>2{&ell;(<i>&beta;</i>&#770;<sub>K,T</sub>) &minus; &ell;(<i>&beta;</i>)}</b>. ",
        "Values at or below <b>&chi;<sup>2</sup><sub>1,1&minus;&alpha;</sub></b> form the likelihood-ratio confidence region. ",
        "Current region: <b>", fmt_interval(est$lr_ci), "</b>; phase report: ",
        phase_badge_html(est$lr_phase), "."
      ))
    )
  })

  output$score_plot <- renderPlotly({
    est <- estimates()$hazard
    req(est)

    theta_grid <- seq(est$theta_bounds[1], est$theta_bounds[2], length.out = 700L)
    beta_grid <- exp(theta_grid)
    score_grid <- vapply(theta_grid, est$score, numeric(1))
    keep <- is.finite(beta_grid) & is.finite(score_grid) & beta_grid > 0

    validate(need(sum(keep) >= 2L, "The score curve could not be evaluated on the selected beta range."))

    beta_plot <- beta_grid[keep]
    score_plot <- score_grid[keep]
    score_hat <- as.numeric(est$score_at_hat)
    if (!is.finite(score_hat)) score_hat <- 0

    y_range <- range(c(score_plot, 0, score_hat), finite = TRUE)
    if (!all(is.finite(y_range)) || diff(y_range) <= 0) y_range <- c(-1, 1)
    y_pad <- 0.08 * diff(y_range)
    y_range <- y_range + c(-y_pad, y_pad)

    score_df <- data.frame(
      beta = beta_plot,
      score = score_plot,
      hover = paste0(
        "<i>&beta;</i> = ", fmt_num(beta_plot, 5),
        "<br><i>U</i>(log <i>&beta;</i>) = ", fmt_num(score_plot, 5)
      ),
      stringsAsFactors = FALSE
    )

    p <- plot_ly(
      data = score_df,
      x = ~beta,
      y = ~score,
      type = "scatter",
      mode = "lines",
      line = list(color = APP_COLORS$blue, width = 3),
      text = ~hover,
      hoverinfo = "text",
      name = "Score U(log beta)"
    ) %>%
      add_segments(
        x = min(beta_plot), xend = max(beta_plot),
        y = 0, yend = 0,
        line = list(color = APP_COLORS$slate, width = 2, dash = "dash"),
        hoverinfo = "skip",
        name = "Zero-score reference"
      ) %>%
      add_segments(
        x = est$estimate, xend = est$estimate,
        y = y_range[1], yend = y_range[2],
        line = list(color = APP_COLORS$red, width = 3, dash = "dot"),
        text = paste0("Constrained MLE <i>&beta;</i> = ", fmt_num(est$estimate, 5)),
        hoverinfo = "text",
        name = "Numerical maximizer"
      ) %>%
      add_markers(
        x = est$estimate,
        y = score_hat,
        marker = list(
          size = 11,
          color = "white",
          line = list(color = APP_COLORS$red, width = 3)
        ),
        text = paste0(
          "<b>Numerical maximizer</b><br>",
          "<i>&beta;</i>&#770;<sub>K,T</sub> = ", fmt_num(est$estimate, 5),
          "<br><i>U</i>(log <i>&beta;</i>&#770;<sub>K,T</sub>) = ", fmt_num(score_hat, 5),
          "<br>Status: ", est$boundary_flag
        ),
        hoverinfo = "text",
        showlegend = FALSE
      ) %>%
      layout(
        xaxis = list(
          title = "Parameter β (logarithmic axis)",
          type = "log",
          zeroline = FALSE
        ),
        yaxis = list(
          title = "Score U(log β)",
          range = y_range,
          zeroline = FALSE
        ),
        legend = list(orientation = "h", y = -0.22),
        hoverlabel = list(align = "left"),
        margin = list(l = 80, r = 25, t = 20, b = 90)
      )

    plotly_clean(p)
  })

  output$score_plot_explanation <- renderUI({
    est <- estimates()$hazard
    req(est)
    status_text <- if (est$boundary_flag == "Interior") {
      "The score changes sign inside the selected parameter interval, so the red line marks the unique unconstrained maximum."
    } else {
      paste0(
        "The score does not change sign inside the selected parameter interval; therefore the constrained maximizer is the ",
        tolower(est$boundary_flag), "."
      )
    }

    tags$div(
      class = "small text-muted",
      HTML(paste0(
        status_text,
        " The plotted function is <b><i>U</i><sub>K,T</sub>(&theta;) = <i>D</i><sup>+</sup><sub>K,T</sub> &minus; &sum;<sub>j=0</sub><sup>T&minus;1</sup> <i>R</i><sub>K,T,j</sub><i>h</i><sub>&theta;</sub>(j)</b>, ",
        "with <b>&theta; = log(&beta;)</b>. Its derivative is ",
        "<b><i>U</i>&prime;<sub>K,T</sub>(&theta;) = &minus;&#119973;<sub>K,T</sub>(&theta;) &lt; 0</b> whenever observed information is positive."
      ))
    )
  })


  risk_frame <- reactive({
    dat <- analysis_data()
    est <- estimates()
    beta_ref <- reference_beta()
    req(dat, est)

    j <- 0:(dat$T - 1L)
    h_fit <- est$hazard$estimate / (est$hazard$estimate + j + 1)
    h_ref <- beta_ref / (beta_ref + j + 1)
    empirical <- ifelse(dat$R > 0, dat$d / dat$R, NA_real_)
    info <- dat$R * h_fit * (1 - h_fit)

    data.frame(
      Age = j,
      Risk_R = as.numeric(dat$R),
      Deaths_d = as.numeric(dat$d),
      Jumps = as.numeric(dat$R - dat$d),
      Empirical_death_fraction = empirical,
      Fitted_hazard = h_fit,
      Reference_hazard = h_ref,
      Information_contribution = info,
      stringsAsFactors = FALSE
    )
  })

  output$hazard_plot <- renderPlotly({
    rf <- risk_frame()
    observed <- subset(rf, Risk_R > 0)
    size_values <- 7 + 18 * sqrt(observed$Risk_R / max(observed$Risk_R))

    p <- plot_ly() %>%
      add_markers(
        data = observed,
        x = ~Age,
        y = ~Empirical_death_fraction,
        marker = list(color = APP_COLORS$navy, size = size_values, opacity = 0.75, line = list(color = "white", width = 1)),
        text = ~paste0(
          "<b>Age ", Age, "</b><br>",
          "At risk R_j: ", fmt_int(Risk_R), "<br>",
          "Deaths d_j: ", fmt_int(Deaths_d), "<br>",
          "Empirical death fraction: ", fmt_num(Empirical_death_fraction)
        ),
        hoverinfo = "text",
        name = "Observed d_j / R_j"
      ) %>%
      add_lines(
        data = rf,
        x = ~Age,
        y = ~Fitted_hazard,
        line = list(color = APP_COLORS$purple, width = 3),
        text = ~paste0("Age ", Age, "<br>Fitted hazard: ", fmt_num(Fitted_hazard)),
        hoverinfo = "text",
        name = "Fitted hazard"
      ) %>%
      add_lines(
        data = rf,
        x = ~Age,
        y = ~Reference_hazard,
        line = list(color = APP_COLORS$teal, width = 2, dash = "dash"),
        text = ~paste0("Age ", Age, "<br>Reference hazard: ", fmt_num(Reference_hazard)),
        hoverinfo = "text",
        name = "Reference hazard"
      ) %>%
      layout(
        xaxis = list(title = "Individual age j", dtick = 1),
        yaxis = list(title = "Death probability", range = c(0, min(1, max(c(observed$Empirical_death_fraction, rf$Fitted_hazard, rf$Reference_hazard), na.rm = TRUE) * 1.15))),
        legend = list(orientation = "h", y = -0.23),
        margin = list(l = 70, r = 20, t = 20, b = 85)
      )

    plotly_clean(p)
  })

  output$information_plot <- renderPlotly({
    rf <- risk_frame()
    info_total <- sum(rf$Information_contribution)
    rf$Share <- if (info_total > 0) rf$Information_contribution / info_total else 0

    p <- plot_ly(
      rf,
      x = ~Age,
      y = ~Information_contribution,
      type = "bar",
      marker = list(color = APP_COLORS$teal),
      text = ~paste0(
        "<b>Age ", Age, "</b><br>",
        "Risk rows: ", fmt_int(Risk_R), "<br>",
        "Information: ", fmt_num(Information_contribution), "<br>",
        "Share of total: ", fmt_pct(Share)
      ),
      hoverinfo = "text",
      name = "Observed information"
    ) %>%
      layout(
        xaxis = list(title = "Individual age j", dtick = 1),
        yaxis = list(title = "Observed-information contribution", rangemode = "tozero"),
        margin = list(l = 75, r = 20, t = 20, b = 60)
      )

    plotly_clean(p)
  })

  output$risk_table <- renderDT({
    rf <- risk_frame()
    display <- transform(
      rf,
      Empirical_death_fraction = ifelse(is.na(Empirical_death_fraction), "—", fmt_num(Empirical_death_fraction)),
      Fitted_hazard = fmt_num(Fitted_hazard),
      Reference_hazard = fmt_num(Reference_hazard),
      Information_contribution = fmt_num(Information_contribution)
    )
    datatable(
      display,
      rownames = FALSE,
      options = list(pageLength = 12, scrollX = TRUE, dom = "tip"),
      class = "stripe hover compact"
    )
  })

  output$dynamics_availability <- renderUI({
    dat <- analysis_data()
    req(dat)
    if (dat$source != "Exact simulation") {
      tags$div(
        class = "warning-note mb-3",
        "Process-dynamics graphics require the exact simulator. An imported age-risk table is sufficient for inference but does not retain calendar-time trajectories or labeled transition logs."
      )
    } else {
      tags$div(
        class = "explanation mb-3",
        paste0(
          "Detailed dynamics were retained for the first ", dat$detail_n,
          " realization(s). Aggregate inference still uses all ", dat$K, " realizations."
        )
      )
    }
  })

  output$realization_selector <- renderUI({
    dat <- analysis_data()
    req(dat, dat$source == "Exact simulation")
    sliderInput(
      "selected_realization",
      "Realization to inspect",
      min = 1,
      max = max(1, dat$detail_n),
      value = 1,
      step = 1,
      ticks = FALSE
    )
  })

  selected_realization <- reactive({
    dat <- analysis_data()
    req(dat, dat$source == "Exact simulation", input$selected_realization)
    r <- min(max(1L, as.integer(input$selected_realization)), dat$detail_n)
    list(
      r = r,
      summary = dat$per_realization[dat$per_realization$realization == r, , drop = FALSE],
      dynamics = dat$dynamics[[r]],
      events = dat$events[[r]]
    )
  })

  output$realization_metrics <- renderUI({
    sr <- selected_realization()
    s <- sr$summary
    extinction_text <- if (is.na(s$extinction_time)) "Not extinct by T" else paste0("Time ", s$extinction_time)
    tagList(
      tags$hr(),
      tags$div(class = "small-caps", "Visible first transitions"),
      tags$div(class = "value", paste0("N = ", s$first_transition_rows_N, ", M = ", s$first_transition_jumps_M)),
      tags$hr(),
      tags$div(class = "small-caps", "All visible rows"),
      tags$div(class = "value", paste0("B = ", s$visible_rows_B)),
      tags$hr(),
      tags$div(class = "small-caps", "Extinction status"),
      tags$div(class = "value", extinction_text),
      tags$hr(),
      tags$div(class = "small-caps", "Activated by terminal time"),
      tags$div(class = "value", fmt_int(s$activated_by_T))
    )
  })

  output$dynamics_plot <- renderPlotly({
    sr <- selected_realization()
    dyn <- sr$dynamics

    long <- rbind(
      data.frame(time = dyn$time, value = dyn$alive, Series = "Activated and alive"),
      data.frame(time = dyn$time, value = dyn$cumulative_activated, Series = "Cumulative activated")
    )

    p <- plot_ly(
      long,
      x = ~time,
      y = ~value,
      color = ~Series,
      colors = c("Activated and alive" = APP_COLORS$purple, "Cumulative activated" = APP_COLORS$teal),
      type = "scatter",
      mode = "lines+markers",
      line = list(width = 3),
      marker = list(size = 6),
      text = ~paste0("Time: ", time, "<br>", Series, ": ", fmt_int(value)),
      hoverinfo = "text"
    ) %>%
      layout(
        xaxis = list(title = "Calendar time", dtick = 1),
        yaxis = list(title = "Number of particles", rangemode = "tozero"),
        legend = list(orientation = "h", y = -0.22),
        margin = list(l = 70, r = 20, t = 20, b = 80)
      )

    plotly_clean(p)
  })

  output$event_flow_plot <- renderPlotly({
    sr <- selected_realization()
    dyn <- sr$dynamics
    flow <- subset(dyn, time > 0)
    flow_long <- rbind(
      data.frame(time = flow$time, Count = flow$jumps, Event = "Jumps"),
      data.frame(time = flow$time, Count = flow$deaths, Event = "Deaths"),
      data.frame(time = flow$time, Count = flow$new_activations, Event = "New activations")
    )

    p <- plot_ly(
      flow_long,
      x = ~time,
      y = ~Count,
      color = ~Event,
      colors = c("Jumps" = APP_COLORS$blue, "Deaths" = APP_COLORS$red, "New activations" = APP_COLORS$green),
      type = "bar",
      text = ~paste0("Calendar time: ", time, "<br>", Event, ": ", Count),
      hoverinfo = "text"
    ) %>%
      layout(
        barmode = "group",
        xaxis = list(title = "Arrival / state-record time", dtick = 1),
        yaxis = list(title = "Event count", rangemode = "tozero"),
        legend = list(orientation = "h", y = -0.25),
        margin = list(l = 65, r = 20, t = 20, b = 85)
      )

    plotly_clean(p)
  })

  output$event_table <- renderDT({
    sr <- selected_realization()
    events <- sr$events
    if (nrow(events) == 0L) {
      events <- data.frame(Message = "No visible transitions were retained.")
    } else {
      events$death_hazard <- vapply(events$death_hazard, fmt_num, character(1))
    }
    datatable(
      events,
      rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE, dom = "tip"),
      class = "stripe hover compact"
    )
  })

  mc_results <- eventReactive(input$run_mc, {
    tryCatch({
      if (input$data_mode != "simulate") {
        stop("Repeated-sampling exploration is available only in exact-simulation mode.")
      }
      B <- as.integer(input$mc_reps)
      K <- as.integer(input$K)
      T <- as.integer(input$T)
      beta <- as.numeric(input$beta_value)
      if (B < 1L || K < 1L || T < 1L || beta <= 0) stop("Invalid repeated-sampling controls.")

      withProgress(message = "Repeated-sampling exploration", value = 0, {
        run_repeated_sampling(
          B = B,
          K = K,
          T = T,
          beta = beta,
          confidence_level = input$confidence_level,
          beta_bounds = c(input$beta_lower, input$beta_upper),
          seed = input$mc_seed,
          progress = function(amount, detail) {
            incProgress(amount, detail = detail)
          }
        )
      })
    }, error = function(e) {
      showNotification(conditionMessage(e), type = "error", duration = 10)
      NULL
    })
  }, ignoreNULL = TRUE)

  output$mc_distribution_plot <- renderPlotly({
    mc <- mc_results()
    validate(need(!is.null(mc), "Run the repeated-sampling study to populate this figure."))

    long <- rbind(
      data.frame(Method = "Origin", Estimate = mc$draws$root_estimate),
      data.frame(Method = "Activated", Estimate = mc$draws$pooled_estimate),
      data.frame(Method = "Full hazard", Estimate = mc$draws$hazard_estimate)
    )
    long <- subset(long, is.finite(Estimate))

    p <- ggplot(long, aes(x = Method, y = Estimate, fill = Method)) +
      geom_violin(alpha = 0.35, trim = TRUE, color = NA) +
      geom_boxplot(width = 0.17, outlier.alpha = 0.25, color = APP_COLORS$navy) +
      geom_hline(yintercept = mc$beta, linetype = "dotdash", linewidth = 1.1, color = APP_COLORS$teal) +
      geom_hline(yintercept = BETA_CRITICAL, linetype = "dashed", linewidth = 1, color = APP_COLORS$red) +
      scale_fill_manual(values = c("Origin" = APP_COLORS$slate, "Activated" = APP_COLORS$teal, "Full hazard" = APP_COLORS$purple), guide = "none") +
      labs(x = NULL, y = "Point estimate β") +
      theme_minimal(base_size = 13) +
      theme(panel.grid.major.x = element_blank(), axis.text.x = element_text(face = "bold"))

    suppressWarnings(ggplotly(p, tooltip = c("x", "y"))) %>% plotly_clean()
  })

  output$mc_coverage_plot <- renderPlotly({
    mc <- mc_results()
    validate(need(!is.null(mc), "Run the repeated-sampling study to populate this figure."))

    s <- mc$summary
    cov <- data.frame(Method = s$Method, Coverage = s$Coverage)

    p <- plot_ly(
      cov,
      x = ~Coverage,
      y = ~reorder(Method, Coverage),
      type = "bar",
      orientation = "h",
      marker = list(color = c(APP_COLORS$slate, APP_COLORS$teal, APP_COLORS$blue, APP_COLORS$purple)),
      text = ~paste0("<b>", Method, "</b><br>Empirical coverage: ", fmt_pct(Coverage)),
      hoverinfo = "text"
    ) %>%
      add_segments(
        x = mc$confidence_level, xend = mc$confidence_level,
        y = -0.5, yend = nrow(cov) - 0.5,
        line = list(color = APP_COLORS$red, dash = "dash", width = 2),
        text = paste0("Requested level: ", fmt_pct(mc$confidence_level)),
        hoverinfo = "text",
        name = "Requested level"
      ) %>%
      layout(
        xaxis = list(title = "Empirical coverage", range = c(0, 1), tickformat = ".0%"),
        yaxis = list(title = ""),
        margin = list(l = 210, r = 20, t = 20, b = 60),
        showlegend = FALSE
      )

    plotly_clean(p)
  })

  output$mc_diagnostics <- renderUI({
    mc <- mc_results()
    validate(need(!is.null(mc), "Run the repeated-sampling study to populate diagnostics."))

    tagList(
      tags$div(class = "result-callout", tags$div(class = "small-caps", "Data sets"), tags$div(class = "value", fmt_int(mc$B))),
      tags$br(),
      tags$div(class = "result-callout", tags$div(class = "small-caps", "Hazard boundary-fit rate"), tags$div(class = "value", fmt_pct(mc$boundary_rate))),
      tags$br(),
      tags$div(class = "result-callout", tags$div(class = "small-caps", "Mean age-zero rows per realization"), tags$div(class = "value", fmt_num(mc$mean_N_per_realization))),
      tags$br(),
      tags$div(class = "result-callout", tags$div(class = "small-caps", "Mean visible rows per realization"), tags$div(class = "value", fmt_num(mc$mean_rows_per_realization))),
      tags$br(),
      tags$div(
        class = "small text-muted",
        if (mc$beta == BETA_CRITICAL) {
          "At β = 1/2, the key phase diagnostic is the decisive-output rate rather than a correct-phase rate."
        } else {
          "Away from β = 1/2, the correct-phase column estimates the probability that an interval lies wholly on the true side of the threshold."
        }
      )
    )
  })

  output$mc_summary_table <- renderDT({
    mc <- mc_results()
    validate(need(!is.null(mc), "Run the repeated-sampling study to populate this table."))
    display <- mc$summary
    display$Bias <- vapply(display$Bias, fmt_num, character(1))
    display$RMSE <- vapply(display$RMSE, fmt_num, character(1))
    display$Median <- vapply(display$Median, fmt_num, character(1))
    display$Coverage <- vapply(display$Coverage, fmt_pct, character(1))
    display$Mean_finite_length <- vapply(display$Mean_finite_length, fmt_num, character(1))
    display$Unbounded_or_unavailable <- vapply(display$Unbounded_or_unavailable, fmt_pct, character(1))
    display$Correct_phase <- vapply(display$Correct_phase, fmt_pct, character(1))
    display$Decisive_rate <- vapply(display$Decisive_rate, fmt_pct, character(1))

    names(display) <- c(
      "Method", "Bias", "RMSE", "Median", "Coverage",
      "Mean finite interval length", "Unbounded / unavailable",
      "Correct phase", "Decisive rate"
    )

    datatable(
      display,
      rownames = FALSE,
      options = list(pageLength = 8, dom = "tip", scrollX = TRUE),
      class = "stripe hover compact"
    )
  })

  output$validation_table <- renderDT({
    dat <- analysis_data()
    est <- estimates()
    req(dat, est)
    checks <- validate_dataset(dat, est)
    datatable(
      checks,
      escape = FALSE,
      rownames = FALSE,
      options = list(pageLength = 20, dom = "tip", scrollX = TRUE),
      class = "stripe hover compact"
    )
  })

  output$per_realization_table <- renderDT({
    dat <- analysis_data()
    req(dat, dat$source == "Exact simulation")
    pr <- dat$per_realization
    pr$extinction_time <- ifelse(is.na(pr$extinction_time), "Not extinct by T", pr$extinction_time)
    datatable(
      pr,
      rownames = FALSE,
      filter = "top",
      options = list(pageLength = 12, scrollX = TRUE),
      class = "stripe hover compact"
    )
  })

  output$download_risk <- downloadHandler(
    filename = function() {
      dat <- analysis_data()
      paste0("frog_model_risk_table_K", dat$K, "_T", dat$T, ".csv")
    },
    content = function(file) {
      utils::write.csv(
        data.frame(age = 0:(analysis_data()$T - 1L), R = analysis_data()$R, d = analysis_data()$d),
        file,
        row.names = FALSE
      )
    }
  )

  output$download_summary <- downloadHandler(
    filename = function() {
      dat <- analysis_data()
      paste0("frog_model_inference_summary_K", dat$K, "_T", dat$T, ".csv")
    },
    content = function(file) {
      out <- estimator_summary_frame()
      out$Phase_decision <- gsub("<[^>]+>", "", out$Phase_decision)
      utils::write.csv(out, file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
