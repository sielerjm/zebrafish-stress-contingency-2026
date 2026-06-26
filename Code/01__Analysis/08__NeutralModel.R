# 08__NeutralModel.R
# Created by: Michael Sieler
# Date last updated: 2026-04-24
#
# Description: Sloan neutral community model (Sloan et al. 2006; Burns et al. 2016) — global
#   metacommunity fit, partitions, bootstrap CIs, Culicoidibacter focal figures, and regime-wise
#   partition summaries (Sloan fit + regime partitions are inlined below). Writes to
#   Results/08__NeutralModel/{Figures,Tables,Stats}.
#
# Manuscript display: Code/02__Results/08__NeutralModel.Rmd loads Stats/neutral_model__bundle.rds.
#
# Expected input: Run from Sieler2026 project root; Data/r_objects/ps-list__*.rds from
#   04__DataPreProcess.R (ps.list loaded via 00__InitializeEnvironment.R).
# Expected output: Results/08__NeutralModel/Figures, Tables, Stats; bundle RDS for Rmd.

init_rel <- file.path("Code", "00__Setup", "00__InitializeEnvironment.R")
if (!file.exists(init_rel)) {
  stop(
    "Set working directory to the Sieler2026 repository root, or run:\n",
    "  Rscript Code/01__Analysis/08__NeutralModel.R\n",
    "from the project root so `", init_rel, "` resolves."
  )
}
source(init_rel)

path_res <- file.path(path.results, "08__NeutralModel")
path_fig <- file.path(path_res, "Figures")
path_tbl <- file.path(path_res, "Tables")
path_rds <- file.path(path_res, "Stats")
out_dir <- path_res
dir.create(path_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(path_tbl, recursive = TRUE, showWarnings = FALSE)
dir.create(path_rds, recursive = TRUE, showWarnings = FALSE)

# Manuscript Fig 6.1 ggplot (Genus focal four); set when the Genus focal panel is built.
p_neutral_focal_four_genus <- NULL

# Archive previous outputs (copy → verify → clear) so each run is clean but reproducible.
sieler2026_archive_module_outputs(
  path_res_module = path_res,
  module_name = "08__NeutralModel",
  subdirs = c("Figures", "Tables", "Stats"),
  clear_original = TRUE
)

# --- Analysis constants (match prior testing pipeline; day 60 only for production) ---------------
target_time <- 60L
neutral_model_all_times <- FALSE
tax_rank_runs <- list(ASV = NULL, Genus = "Genus")
highlight_genus_pattern <- "Culicoidibacter"
# Additional focal genera for reviewer appendix outputs (plots/tables in inlined sections below)
focal_genera <- c("Shewanella", "Culicoidibacter", "Flavobacterium", "Cetobacterium")
.focal_four_taxon_outline_colors <- c(
  Culicoidibacter = "#E41A1C",
  Shewanella = "#377EB8",
  Flavobacterium = "#4DAF4A",
  Cetobacterium = "#984EA3"
)
partition_colors <- partition_colors_neutral_model
n_bootstrap <- 1000L
equalize_partition_sizes <- TRUE


# =============================================================================
# Inlined: 08__NeutralModel__sloan_fit.R (formerly sourced; same globals as parent)
# =============================================================================

# spp: numeric matrix, samples = rows, taxa = columns; rows must share rarefaction depth.
# light = TRUE: nls + Rsqr only (for bootstrap resamples; skips slow confint/mle).
# Returns list with fit_stats (1-row tibble), predictions (per taxon), or error message.
sncm_fit_burns <- function(spp, pool = NULL, light = FALSE) {
  if (!is.matrix(spp) || !is.numeric(spp)) {
    return(list(success = FALSE, error = "spp must be a numeric matrix", fit_stats = NULL, predictions = NULL))
  }
  if (nrow(spp) < 2L) {
    return(list(success = FALSE, error = "need >= 2 samples", fit_stats = NULL, predictions = NULL))
  }

  # Mean individuals per local community (Burns et al.)
  N <- mean(rowSums(spp))
  if (!is.finite(N) || N <= 0) {
    return(list(success = FALSE, error = "invalid N from row sums", fit_stats = NULL, predictions = NULL))
  }

  # Metacommunity relative abundance p (taxon names preserved)
  if (is.null(pool)) {
    p_m <- colMeans(spp)
  } else {
    p_m <- colMeans(pool)
  }
  p_m <- p_m[p_m != 0]

  # Occurrence frequency (proportion of samples where taxon present)
  spp_bi <- (spp > 0) * 1
  freq <- colMeans(spp_bi)
  freq <- freq[freq != 0]

  # Inner join taxa present in both p and freq with no zeros (Burns et al.)
  df_p <- tibble::tibble(taxon_id = names(p_m), p_mean = as.numeric(p_m))
  df_f <- tibble::tibble(taxon_id = names(freq), freq = as.numeric(freq))
  C <- dplyr::inner_join(df_p, df_f, by = "taxon_id") %>%
    dplyr::filter(p_mean > 0, freq > 0) %>%
    dplyr::arrange(p_mean) %>%
    dplyr::mutate(p = p_mean / N)

  if (nrow(C) < 3L) {
    return(list(success = FALSE, error = "need >= 3 taxa after filtering", fit_stats = NULL, predictions = NULL))
  }

  p <- C$p
  names(p) <- C$taxon_id
  freq <- C$freq
  names(freq) <- C$taxon_id
  n_samples <- nrow(spp)
  d <- 1 / N

  # Nonlinear least squares for migration parameter m
  m_fit <- tryCatch(
    suppressWarnings(
      minpack.lm::nlsLM(
        freq ~ stats::pbeta(d, N * m * p, N * m * (1 - p), lower.tail = FALSE),
        start = list(m = 0.1)
      )
    ),
    error = function(e) NULL
  )
  if (is.null(m_fit)) {
    return(list(success = FALSE, error = "nlsLM failed", fit_stats = NULL, predictions = NULL))
  }

  m_est <- stats::coef(m_fit)[["m"]]

  # Fast path for bootstrap: avoid confint profiling and Gaussian MLE blocks (very slow at scale).
  if (isTRUE(light)) {
    freq_pred <- stats::pbeta(d, N * m_est * p, N * m_est * (1 - p), lower.tail = FALSE)
    ss_res <- sum((freq - freq_pred)^2)
    ss_tot <- sum((freq - mean(freq))^2)
    rsqr <- if (ss_tot > 0) 1 - ss_res / ss_tot else NA_real_
    fit_stats <- tibble::tibble(
      m_nls = m_est,
      m_minus_ci_lower = NA_real_,
      m_mle = NA_real_,
      loglik_sncm = NA_real_,
      loglik_bino = NA_real_,
      loglik_pois = NA_real_,
      rsqr = rsqr,
      rsqr_bino = NA_real_,
      rsqr_pois = NA_real_,
      rmse = if (length(freq) > 1L) sqrt(ss_res / (length(freq) - 1L)) else NA_real_,
      rmse_bino = NA_real_,
      rmse_pois = NA_real_,
      aic_sncm = NA_real_,
      bic_sncm = NA_real_,
      aic_bino = NA_real_,
      bic_bino = NA_real_,
      aic_pois = NA_real_,
      bic_pois = NA_real_,
      N = N,
      n_samples = n_samples,
      n_taxa_fit = length(p),
      detect_d = d
    )
    return(list(success = TRUE, error = NULL, fit_stats = fit_stats, predictions = NULL))
  }

  m_ci_mat <- tryCatch(
    suppressWarnings(stats::confint(m_fit, "m", level = 0.95)),
    error = function(e) NULL
  )
  m_ci_lower <- if (is.matrix(m_ci_mat) && nrow(m_ci_mat) >= 1L && ncol(m_ci_mat) >= 1L) {
    v <- suppressWarnings(as.numeric(m_ci_mat[1L, 1L]))
    if (is.finite(v)) v else NA_real_
  } else {
    NA_real_
  }

  # MLE for AIC (Burns et al. likelihood on Gaussian residuals — legacy recipe)
  sncm_ll <- function(m, sigma) {
    R <- freq - stats::pbeta(d, N * m * p, N * m * (1 - p), lower.tail = FALSE)
    R <- stats::dnorm(R, 0, sigma)
    -sum(log(R))
  }
  m_mle <- tryCatch(
    stats4::mle(sncm_ll, start = list(m = 0.1, sigma = 0.1), nobs = length(p)),
    error = function(e) NULL
  )

  aic_fit <- if (!is.null(m_mle)) stats::AIC(m_mle, k = 2) else NA_real_
  bic_fit <- if (!is.null(m_mle)) stats::BIC(m_mle) else NA_real_

  freq_pred <- stats::pbeta(d, N * m_est * p, N * m_est * (1 - p), lower.tail = FALSE)
  ss_res <- sum((freq - freq_pred)^2)
  ss_tot <- sum((freq - mean(freq))^2)
  rsqr <- if (ss_tot > 0) 1 - ss_res / ss_tot else NA_real_
  rmse <- if (length(freq) > 1L) sqrt(ss_res / (length(freq) - 1L)) else NA_real_

  pred_ci <- Hmisc::binconf(freq_pred * n_samples, n_samples, alpha = 0.05, method = "wilson", return.df = TRUE)

  bino_ll <- function(mu, sigma) {
    R <- freq - stats::pbinom(d, N, p, lower.tail = FALSE)
    R <- stats::dnorm(R, mu, sigma)
    -sum(log(R))
  }
  bino_mle <- tryCatch(
    stats4::mle(bino_ll, start = list(mu = 0, sigma = 0.1), nobs = length(p)),
    error = function(e) NULL
  )
  aic_bino <- if (!is.null(bino_mle)) stats::AIC(bino_mle, k = 2) else NA_real_
  bic_bino <- if (!is.null(bino_mle)) stats::BIC(bino_mle) else NA_real_

  bino_pred <- stats::pbinom(d, N, p, lower.tail = FALSE)
  rsqr_bino <- if (ss_tot > 0) 1 - sum((freq - bino_pred)^2) / ss_tot else NA_real_
  rmse_bino <- if (length(freq) > 1L) sqrt(sum((freq - bino_pred)^2) / (length(freq) - 1L)) else NA_real_
  bino_pred_ci <- Hmisc::binconf(bino_pred * n_samples, n_samples, alpha = 0.05, method = "wilson", return.df = TRUE)

  pois_ll <- function(mu, sigma) {
    R <- freq - stats::ppois(d, N * p, lower.tail = FALSE)
    R <- stats::dnorm(R, mu, sigma)
    -sum(log(R))
  }
  pois_mle <- tryCatch(
    stats4::mle(pois_ll, start = list(mu = 0, sigma = 0.1), nobs = length(p)),
    error = function(e) NULL
  )
  aic_pois <- if (!is.null(pois_mle)) stats::AIC(pois_mle, k = 2) else NA_real_
  bic_pois <- if (!is.null(pois_mle)) stats::BIC(pois_mle) else NA_real_

  pois_pred <- stats::ppois(d, N * p, lower.tail = FALSE)
  rsqr_pois <- if (ss_tot > 0) 1 - sum((freq - pois_pred)^2) / ss_tot else NA_real_
  rmse_pois <- if (length(freq) > 1L) sqrt(sum((freq - pois_pred)^2) / (length(freq) - 1L)) else NA_real_

  negloglik <- function(mle_obj) {
    if (is.null(mle_obj)) {
      return(NA_real_)
    }
    tryCatch(as.numeric(stats4::logLik(mle_obj)), error = function(e) NA_real_)
  }

  fit_stats <- tibble::tibble(
    m_nls = m_est,
    m_minus_ci_lower = m_est - m_ci_lower,
    m_mle = if (!is.null(m_mle)) as.numeric(stats4::coef(m_mle)[["m"]]) else NA_real_,
    loglik_sncm = negloglik(m_mle),
    loglik_bino = negloglik(bino_mle),
    loglik_pois = negloglik(pois_mle),
    rsqr = rsqr,
    rsqr_bino = rsqr_bino,
    rsqr_pois = rsqr_pois,
    rmse = rmse,
    rmse_bino = rmse_bino,
    rmse_pois = rmse_pois,
    aic_sncm = aic_fit,
    bic_sncm = bic_fit,
    aic_bino = aic_bino,
    bic_bino = bic_bino,
    aic_pois = aic_pois,
    bic_pois = bic_pois,
    N = N,
    n_samples = n_samples,
    n_taxa_fit = length(p),
    detect_d = d
  )

  predictions <- tibble::tibble(
    taxon_id = names(p),
    p = as.numeric(p),
    freq = as.numeric(freq),
    freq_pred = as.numeric(freq_pred),
    pred_lwr = pred_ci[, "Lower"],
    pred_upr = pred_ci[, "Upper"],
    bino_pred = as.numeric(bino_pred),
    bino_lwr = bino_pred_ci[, "Lower"],
    bino_upr = bino_pred_ci[, "Upper"],
    partition = dplyr::case_when(
      freq > pred_ci[, "Upper"] ~ "above",
      freq < pred_ci[, "Lower"] ~ "below",
      TRUE ~ "neutral"
    )
  )

  list(success = TRUE, error = NULL, fit_stats = fit_stats, predictions = predictions)
}

# --- Bootstrap: resample samples (rows) with replacement per stratum -------------------------
bootstrap_sncm_stats <- function(spp, n_boot = 1000L) {
  n <- nrow(spp)
  if (n < 2L) {
    return(tibble::tibble())
  }
  out <- vector("list", n_boot)
  for (b in seq_len(n_boot)) {
    set.seed(42L + b)
    idx <- sample.int(n, size = n, replace = TRUE)
    res <- sncm_fit_burns(spp[idx, , drop = FALSE], light = TRUE)
    if (!isTRUE(res$success) || is.null(res$fit_stats)) {
      out[[b]] <- tibble::tibble(m_nls = NA_real_, rsqr = NA_real_)
    } else {
      out[[b]] <- res$fit_stats %>%
        dplyr::select(m_nls, rsqr)
    }
  }
  dplyr::bind_rows(out)
}

# --- Equalize partition lists to smallest set size (Burns et al. 2016) ------------------------
equalize_partition_taxa <- function(partition_ids, seed = 42L) {
  sizes <- lengths(partition_ids)
  n_min <- min(sizes)
  if (!is.finite(n_min) || n_min < 1L) {
    return(partition_ids)
  }
  set.seed(seed)
  purrr::map(partition_ids, function(ids) {
    if (length(ids) <= n_min) {
      ids
    } else {
      sample(ids, n_min)
    }
  })
}

# --- Annotate predictions for highlight taxon (ASV via tax_table Genus; Genus run via taxon_id) --
annotate_predictions_highlight <- function(ps_rare, predictions, tax_rank, pattern) {
  if (is.null(predictions) || nrow(predictions) < 1L) {
    return(predictions)
  }
  rx <- stringr::regex(pattern, ignore_case = TRUE)
  if (!is.null(tax_rank)) {
    return(predictions %>%
      dplyr::mutate(
        genus_label = as.character(.data$taxon_id),
        is_culici = stringr::str_detect(.data$taxon_id, rx)
      ))
  }
  tt <- as(phyloseq::tax_table(ps_rare), "matrix")
  taxa <- phyloseq::taxa_names(ps_rare)
  if (!"Genus" %in% colnames(tt)) {
    return(predictions %>%
      dplyr::mutate(genus_label = NA_character_, is_culici = stringr::str_detect(.data$taxon_id, rx)))
  }
  gen <- as.character(tt[, "Genus"])
  names(gen) <- taxa
  predictions %>%
    dplyr::mutate(
      genus_label = unname(gen[as.character(.data$taxon_id)]),
      genus_label = dplyr::if_else(is.na(.data$genus_label), "", .data$genus_label),
      is_culici = stringr::str_detect(.data$genus_label, rx) | stringr::str_detect(.data$taxon_id, rx)
    )
}

# --- Manual scale for partition colors (ASV + Genus) -----------------------------------
# partition_colour_guide: optional guide_legend(); use when a second colour scale follows
#   (e.g. ggnewscale::new_scale_color) so legend styling binds to this scale, not the last one.
apply_partition_color_scale <- function(
    p,
    partition_colors,
    partition_colour_guide = NULL) {
  if (is.null(partition_colors)) {
    return(p)
  }
  part_lab <- c(above = "Above", below = "Below", neutral = "Neutral")
  brk <- names(partition_colors)
  lbls <- unname(part_lab[brk])
  na_l <- is.na(lbls)
  lbls[na_l] <- brk[na_l]
  sc_args <- list(
    values = partition_colors,
    breaks = brk,
    labels = lbls,
    drop = FALSE
  )
  if (!is.null(partition_colour_guide)) {
    sc_args$guide <- partition_colour_guide
  }
  p + do.call(ggplot2::scale_color_manual, sc_args)
}

# Fills for partition *colour*-scale legend keys (shape 21 + black rim): Pastel1 red/blue for
# above/below; neutral uses grey85 (not Pastel1[[3]], which reads green).
partition_legend_soft_fills <- function(partition_breaks) {
  if (length(partition_breaks) < 1L) {
    return(NULL)
  }
  pal <- RColorBrewer::brewer.pal(9L, "Pastel1")
  key <- c(
    above = pal[[1L]],
    below = pal[[2L]],
    neutral = "grey85"
  )
  out <- unname(key[as.character(partition_breaks)])
  if (any(is.na(out))) {
    return(NULL)
  }
  out
}

# --- Bottom, horizontal legends (neutral-model figures + bar stacks) -------------------------
neutral_plot_legend_theme <- function() {
  ggplot2::theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.title = ggplot2::element_text(face = "bold")
  )
}

# --- Facet strips aligned with manuscript figures (white fill, black frame, bold labels) -----
manuscript_facet_strip_theme <- function() {
  ggplot2::theme(
    strip.background = ggplot2::element_rect(fill = "white", colour = "black", linewidth = 0.4),
    strip.text = ggplot2::element_text(face = "bold", colour = "black"),
    panel.spacing.x = grid::unit(1.1, "lines")
  )
}

# --- Neutral-model plot + optional Culicoidibacter highlight layer -----------------------------
plot_neutral_fitted <- function(
  pred,
  curve_df,
  tv,
  m_est,
  tag,
  highlight = FALSE,
  highlight_label = "Culicoidibacter",
  partition_colors = NULL,
  rsqr = NULL
) {
  if (isTRUE(highlight) && "is_culici" %in% names(pred)) {
    pal_pastel <- RColorBrewer::brewer.pal(9L, "Pastel1")
    fill_part <- c(
      above = pal_pastel[[1L]],
      below = pal_pastel[[2L]],
      neutral = pal_pastel[[3L]]
    )
    shp_part <- c(above = 24L, below = 21L, neutral = 22L)

    pred_bg <- pred %>% dplyr::filter(!.data$is_culici)
    pred_hi <- pred %>% dplyr::filter(.data$is_culici)

    p <- ggplot2::ggplot(pred, ggplot2::aes(x = .data$p, y = .data$freq))
    if (all(c("freq_lwr", "freq_upr") %in% names(curve_df))) {
      p <- p +
        ggplot2::geom_line(
          data = curve_df,
          ggplot2::aes(x = .data$p, y = .data$freq_lwr),
          inherit.aes = FALSE,
          linetype = "dashed",
          color = "grey25",
          linewidth = 0.78,
          show.legend = FALSE
        ) +
        ggplot2::geom_line(
          data = curve_df,
          ggplot2::aes(x = .data$p, y = .data$freq_upr),
          inherit.aes = FALSE,
          linetype = "dashed",
          color = "grey25",
          linewidth = 0.78,
          show.legend = FALSE
        )
    }
    p <- p +
      ggplot2::geom_line(
        data = curve_df,
        ggplot2::aes(x = .data$p, y = .data$freq_pred),
        inherit.aes = FALSE,
        linewidth = 1.05
      ) +
      ggplot2::geom_point(
        data = pred_bg,
        ggplot2::aes(color = .data$partition),
        alpha = 0.25,
        size = 0.9
      ) +
      ggplot2::scale_fill_manual(
        values = fill_part,
        limits = names(fill_part),
        drop = FALSE,
        guide = "none"
      ) +
      ggplot2::scale_shape_manual(values = shp_part, guide = "none") +
      ggplot2::scale_x_log10() +
      ggplot2::labs(
        title = "Sloan Neutral Model (Day 60)",
        subtitle = NULL,
        x = "log(Mean Relative Abundance)",
        y = "Occurrence frequency",
        color = "Model Predictions"
      ) +
      theme_sieler2026_publication(base_size = 14, legend_position = "bottom")
    p <- apply_partition_color_scale(p, partition_colors)
    model_pred_leg_fill <- if (!is.null(partition_colors)) {
      partition_legend_soft_fills(names(partition_colors))
    } else {
      NULL
    }
    model_pred_leg_override <- list(
      alpha = 1,
      size = 3.2,
      shape = 21L,
      colour = "black",
      stroke = 0.35
    )
    if (!is.null(model_pred_leg_fill)) {
      model_pred_leg_override$fill <- model_pred_leg_fill
    }
    p <- p +
      neutral_plot_legend_theme() +
      ggplot2::guides(
        color = ggplot2::guide_legend(
          order = 1L,
          nrow = 1L,
          byrow = TRUE,
          override.aes = model_pred_leg_override
        )
      ) +
      ggplot2::theme(
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.spacing.y = grid::unit(0.35, "cm"),
        legend.box.spacing = grid::unit(0.45, "cm"),
        plot.margin = ggplot2::margin(t = 8, r = 10, b = 10, l = 10, unit = "pt")
      ) +
      ggplot2::coord_cartesian(clip = "off")

    n_lab <- sum(pred$is_culici, na.rm = TRUE)
    if (n_lab > 0L && n_lab <= 40L) {
      pred_hi_lab <- pred_hi %>%
        dplyr::mutate(
          .lab_y = .data$freq - 0.09,
          .lab_x = 1e-2,
          .lab_x_right = .data$.lab_x * 10^(
            0.045 * nchar(as.character(.data$taxon_id))
          ),
          .seg_x0 = 10^(
            (log10(.data$.lab_x) + log10(.data$.lab_x_right)) / 2
          ),
          .seg_y0 = .data$.lab_y + 0.028,
          .k_edge = 0.38,
          .seg_x1 = 10^(
            log10(.data$p) -
              .data$.k_edge * (log10(.data$p) - log10(.data$.seg_x0))
          ),
          .seg_y1 = .data$freq -
            .data$.k_edge * (.data$freq - .data$.seg_y0)
        )
      p <- p +
        ggplot2::geom_segment(
          data = pred_hi_lab,
          ggplot2::aes(
            x = .data$.seg_x0,
            y = .data$.seg_y0,
            xend = .data$.seg_x1,
            yend = .data$.seg_y1
          ),
          inherit.aes = FALSE,
          colour = "black",
          linewidth = 1,
          alpha = 1,
          lineend = "round",
          show.legend = FALSE
        ) +
        ggplot2::geom_text(
          data = pred_hi_lab,
          ggplot2::aes(
            x = .data$.lab_x,
            y = .data$.lab_y,
            label = as.character(.data$taxon_id)
          ),
          inherit.aes = FALSE,
          hjust = 0,
          size = 5.2,
          fontface = "bold",
          colour = "steelblue",
          show.legend = FALSE
        )
    }
    if (n_lab > 0L) {
      p <- p +
        ggplot2::geom_point(
          data = pred_hi,
          ggplot2::aes(shape = .data$partition, fill = .data$partition),
          color = "black",
          stroke = 0.35,
          size = 5
        )
    }

    if (!is.null(rsqr) && is.finite(rsqr)) {
      r2_tbl <- tibble::tibble(
        x = 1e-6,
        y = 0.99,
        lab = paste0("R² = ", signif(rsqr, 3))
      )
      p <- p + ggplot2::geom_text(
        data = r2_tbl,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$lab),
        inherit.aes = FALSE,
        parse = FALSE,
        fontface = "bold",
        hjust = 0,
        vjust = 1,
        size = 4.2,
        color = "grey20",
        show.legend = FALSE
      )
    }
    return(p)
  }


  p <- ggplot2::ggplot(pred, ggplot2::aes(x = .data$p, y = .data$freq)) +
    ggplot2::geom_linerange(
      ggplot2::aes(ymin = .data$pred_lwr, ymax = .data$pred_upr),
      linewidth = 0.2,
      alpha = 0.35,
      color = "grey50"
    ) +
    ggplot2::geom_line(
      data = curve_df,
      ggplot2::aes(x = .data$p, y = .data$freq_pred),
      inherit.aes = FALSE,
      linewidth = 0.8
    ) +
    ggplot2::geom_point(ggplot2::aes(color = .data$partition), alpha = 0.7, size = 1.8) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(
      title = paste0("Sloan neutral model (Time = ", tv, " d) — ", tag),
      subtitle = paste0("m = ", signif(m_est, 4), "; taxa = ", nrow(pred)),
      x = "Mean relative abundance (metacommunity)",
      y = "Occurrence frequency",
      color = "Partition"
    ) +
    theme_sieler2026_publication(base_size = 14, legend_position = "bottom")

  p <- apply_partition_color_scale(p, partition_colors)
  p + neutral_plot_legend_theme()
}

# --- Single panel: community neutral fit + multiple focal genera (large points + labels) -----
plot_neutral_fitted_multi_focal <- function(
    pred,
    curve_df,
    tv,
    m_est,
    run_tag,
    partition_colors = NULL,
    rsqr = NULL) {
  if (!"is_focal" %in% names(pred)) {
    stop("plot_neutral_fitted_multi_focal: pred must contain is_focal (logical).")
  }
  pred_bg <- pred %>% dplyr::filter(!.data$is_focal)
  pred_hi <- pred %>%
    dplyr::filter(.data$is_focal) %>%
    dplyr::mutate(
      label_color = {
        lc <- unname(.focal_four_taxon_outline_colors[as.character(.data$taxon_id)])
        dplyr::if_else(is.na(lc), "grey35", lc)
      }
    )

  pal_pastel <- RColorBrewer::brewer.pal(9L, "Pastel1")
  fill_part <- c(
    above = pal_pastel[[1L]],
    below = pal_pastel[[2L]],
    neutral = pal_pastel[[3L]]
  )
  shp_part <- c(above = 24L, below = 21L, neutral = 22L)

  p <- ggplot2::ggplot(pred, ggplot2::aes(x = .data$p, y = .data$freq))
  if (all(c("freq_lwr", "freq_upr") %in% names(curve_df))) {
    p <- p +
      ggplot2::geom_line(
        data = curve_df,
        ggplot2::aes(x = .data$p, y = .data$freq_lwr),
        inherit.aes = FALSE,
        linetype = "dashed",
        color = "grey25",
        linewidth = 0.78,
        show.legend = FALSE
      ) +
      ggplot2::geom_line(
        data = curve_df,
        ggplot2::aes(x = .data$p, y = .data$freq_upr),
        inherit.aes = FALSE,
        linetype = "dashed",
        color = "grey25",
        linewidth = 0.78,
        show.legend = FALSE
      )
  }
  p <- p +
    ggplot2::geom_line(
      data = curve_df,
      ggplot2::aes(x = .data$p, y = .data$freq_pred),
      inherit.aes = FALSE,
      linewidth = 1.05
    ) +
    ggplot2::geom_point(
      data = pred_bg,
      ggplot2::aes(color = .data$partition),
      alpha = 0.25,
      size = 0.9
    ) +
    ggplot2::scale_fill_manual(
      values = fill_part,
      limits = names(fill_part),
      drop = FALSE,
      guide = "none"
    ) +
    ggplot2::scale_shape_manual(values = shp_part, guide = "none") +
    ggplot2::scale_x_log10(expand = ggplot2::expansion(mult = c(0.08, 0.32))) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      expand = ggplot2::expansion(mult = c(0.03, 0.18))
    ) +
    ggplot2::labs(
      title = "Sloan Neutral Model (Day 60)",
      subtitle = NULL,
      x = "log(Mean Relative Abundance)",
      y = "Occurrence frequency",
      color = "Model Predictions"
    ) +
    theme_sieler2026_publication(base_size = 14, legend_position = "bottom")
  # One colour scale only for Model Predictions (partition). Focal legend uses new_scale_fill.
  # Legend keys: Pastel1 above/below, grey85 neutral, shape 21 + thin black rim (Focal Taxa row style).
  part_brk <- if (!is.null(partition_colors)) names(partition_colors) else character()
  part_model_guide_override <- list(
    alpha = 1,
    size = 3.2,
    shape = 21L,
    colour = "black",
    stroke = 0.35
  )
  model_pred_leg_fill <- partition_legend_soft_fills(part_brk)
  if (!is.null(model_pred_leg_fill)) {
    part_model_guide_override$fill <- model_pred_leg_fill
  }
  p <- apply_partition_color_scale(
    p,
    partition_colors,
    partition_colour_guide = ggplot2::guide_legend(
      order = 1L,
      nrow = 1L,
      byrow = TRUE,
      override.aes = part_model_guide_override
    )
  )
  focal_taxon_legend_order <- c("Culicoidibacter", "Shewanella", "Flavobacterium", "Cetobacterium")
  faux_focal_legend <- tibble::tibble(
    taxon_id = factor(focal_taxon_legend_order, levels = focal_taxon_legend_order),
    p = stats::median(pred$p, na.rm = TRUE),
    freq = stats::median(pred$freq, na.rm = TRUE)
  )
  focal_fill_vals <- .focal_four_taxon_outline_colors[as.character(focal_taxon_legend_order)]

  p <- p +
    neutral_plot_legend_theme() +
    ggplot2::theme(
      legend.direction = "horizontal",
      legend.box = "vertical",
      legend.spacing.y = grid::unit(0.35, "cm"),
      legend.box.spacing = grid::unit(0.45, "cm"),
      axis.title.x = ggplot2::element_text(face = "bold"),
      axis.title.y = ggplot2::element_text(face = "bold"),
      panel.grid.major = ggplot2::element_line(colour = "grey88", linewidth = 0.35),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 14, r = 22, b = 12, l = 12, unit = "pt")
    )

  if (nrow(pred_hi) > 0L) {
    p <- p +
      ggplot2::geom_point(
        data = pred_hi,
        ggplot2::aes(shape = .data$partition, fill = .data$partition),
        color = "black",
        stroke = 0.35,
        size = 5
      )
    set.seed(42)
    p <- p +
      ggrepel::geom_label_repel(
        data = pred_hi,
        ggplot2::aes(
          x = .data$p,
          y = .data$freq,
          label = as.character(.data$taxon_id)
        ),
        inherit.aes = FALSE,
        colour = pred_hi$label_color,
        fill = ggplot2::alpha("white", 0.94),
        label.size = 0.28,
        size = 3.5,
        fontface = "bold",
        segment.color = ggplot2::alpha("black", 0.55),
        segment.size = 0.45,
        min.segment.length = 0,
        box.padding = ggplot2::unit(1.05, "lines"),
        point.padding = ggplot2::unit(0.95, "lines"),
        max.overlaps = Inf,
        force = 8,
        force_pull = 2.2,
        max.iter = 40000L,
        max.time = 4,
        seed = 42L,
        show.legend = FALSE,
        lineend = "round"
      )
  }

  p <- p +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_point(
      data = faux_focal_legend,
      mapping = ggplot2::aes(x = .data$p, y = .data$freq, fill = .data$taxon_id),
      inherit.aes = FALSE,
      colour = "black",
      stroke = 0.35,
      shape = 21L,
      alpha = 0,
      size = 0.25,
      show.legend = TRUE
    ) +
    ggplot2::scale_fill_manual(
      name = "Focal Taxa",
      values = focal_fill_vals,
      breaks = focal_taxon_legend_order,
      limits = focal_taxon_legend_order,
      drop = FALSE,
      guide = ggplot2::guide_legend(
        order = 2L,
        nrow = 1L,
        byrow = TRUE,
        override.aes = list(
          alpha = 1,
          size = 3.2,
          colour = "black",
          stroke = 0.35,
          shape = 21L
        )
      )
    ) +
    ggplot2::coord_cartesian(clip = "off")

  if (!is.null(rsqr) && is.finite(rsqr)) {
    r2_tbl <- tibble::tibble(
      x = 1e-6,
      y = 0.99,
      lab = paste0("R² = ", signif(rsqr, 3))
    )
    p <- p +
      ggplot2::geom_text(
        data = r2_tbl,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$lab),
        inherit.aes = FALSE,
        parse = FALSE,
        fontface = "bold",
        hjust = 0,
        vjust = 1,
        size = 4.2,
        color = "grey20",
        show.legend = FALSE
      )
  }
  p
}

# --- Neutral curve + only the two focal Culicoidibacter ASVs (above / below) -------------------
plot_neutral_focal_two_asvs <- function(
  pred,
  curve_df,
  tv,
  m_est,
  taxon_above,
  taxon_below,
  partition_colors = NULL
) {
  pal_pastel <- RColorBrewer::brewer.pal(9L, "Pastel1")
  pred <- pred %>%
    dplyr::mutate(
      focal = dplyr::case_when(
        .data$taxon_id == taxon_above ~ paste0(.data$taxon_id, " (above)"),
        .data$taxon_id == taxon_below ~ paste0(.data$taxon_id, " (below)"),
        TRUE ~ NA_character_
      )
    )
  pred2 <- pred %>% dplyr::filter(!is.na(.data$focal))

  lab_above <- paste0(taxon_above, " (above)")
  lab_below <- paste0(taxon_below, " (below)")
  # Pastel1: light red / light blue for above- vs below-neutral focal ASVs (mock-up alignment)
  focal_fill <- stats::setNames(
    c(pal_pastel[[1L]], pal_pastel[[2L]]),
    c(lab_above, lab_below)
  )
  # Triangle (above) + filled circle (below); both use fill = Pastel1 with black outline
  focal_shapes <- stats::setNames(c(24L, 21L), c(lab_above, lab_below))

  # Wilson 95% band on the neutral curve (same binomial logic as per-taxon predictions; Burns-style envelope)
  p <- ggplot2::ggplot(pred, ggplot2::aes(x = .data$p, y = .data$freq))
  if (all(c("freq_lwr", "freq_upr") %in% names(curve_df))) {
    p <- p +
      ggplot2::geom_line(
        data = curve_df,
        ggplot2::aes(x = .data$p, y = .data$freq_lwr),
        inherit.aes = FALSE,
        linetype = "dashed",
        color = "grey25",
        linewidth = 0.55,
        show.legend = FALSE
      ) +
      ggplot2::geom_line(
        data = curve_df,
        ggplot2::aes(x = .data$p, y = .data$freq_upr),
        inherit.aes = FALSE,
        linetype = "dashed",
        color = "grey25",
        linewidth = 0.55,
        show.legend = FALSE
      )
  }
  p <- p +
    ggplot2::geom_line(
      data = curve_df,
      ggplot2::aes(x = .data$p, y = .data$freq_pred),
      inherit.aes = FALSE,
      linewidth = 0.8
    ) +
    ggplot2::geom_point(
      data = pred %>% dplyr::filter(is.na(.data$focal)),
      ggplot2::aes(color = .data$partition),
      alpha = 0.25,
      size = 0.9
    ) +
    ggplot2::geom_point(
      data = pred2,
      ggplot2::aes(shape = .data$focal, fill = .data$focal),
      color = "black",
      stroke = 0.35,
      size = 5
    ) +
    ggplot2::scale_fill_manual(
      name = "Culicoidibacter-associated ASVs",
      values = focal_fill,
      breaks = c(lab_above, lab_below),
      drop = FALSE
    ) +
    ggplot2::scale_shape_manual(values = focal_shapes, guide = "none") +
    ggplot2::scale_x_log10() +
    ggplot2::labs(
      title = paste0("Culicoidibacter ASVs (neutral partitions) — Time ", tv, " d"),
      subtitle = paste0(
        taxon_above, " = above-neutral; ", taxon_below, " = below-neutral | m = ", signif(m_est, 4)
      ),
      x = "Mean relative abundance (metacommunity)",
      y = "Occurrence frequency",
      color = "All ASVs"
    ) +
    theme_sieler2026_publication(base_size = 14, legend_position = "bottom")

  p <- apply_partition_color_scale(p, partition_colors)

  # Partition legend keys: Pastel1 above/below, grey85 neutral, shape 21 + black rim.
  leg_pt <- 3.2
  all_asv_leg_override <- list(
    alpha = 1,
    size = leg_pt,
    shape = 21L,
    colour = "black",
    stroke = 0.35
  )
  all_asv_leg_fill <- if (!is.null(partition_colors)) {
    partition_legend_soft_fills(names(partition_colors))
  } else {
    NULL
  }
  if (!is.null(all_asv_leg_fill)) {
    all_asv_leg_override$fill <- all_asv_leg_fill
  }
  p +
    neutral_plot_legend_theme() +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        order = 1L,
        nrow = 1L,
        byrow = TRUE,
        override.aes = all_asv_leg_override
      ),
      fill = ggplot2::guide_legend(
        order = 2L,
        nrow = 1L,
        byrow = TRUE,
        override.aes = list(
          alpha = 1,
          size = leg_pt,
          shape = c(24L, 21L),
          colour = "black",
          stroke = 0.35
        )
      )
    ) +
    ggplot2::theme(
      legend.direction = "horizontal",
      legend.box = "vertical",
      legend.spacing.y = grid::unit(0.35, "cm"),
      legend.box.spacing = grid::unit(0.45, "cm"),
      axis.title = ggplot2::element_text(size = 12),
      axis.text = ggplot2::element_text(size = 10.5),
      plot.margin = ggplot2::margin(t = 5.5, r = 8, b = 10, l = 8, unit = "pt")
    )
}

# --- Load phyloseq and prepare counts --------------------------------------------------------
if (!exists("ps.list", inherits = TRUE) || is.null(ps.list) || !"All" %in% names(ps.list)) {
  stop("ps.list not found or missing 'All'. Run Code/00__Setup/04__DataPreProcess.R.")
}

ps_all <- ps.list[["All"]]

if (!"Time" %in% colnames(as(phyloseq::sample_data(ps_all), "data.frame"))) {
  stop("Sample metadata must include column 'Time'.")
}

# Metacommunity time strata (one neutral fit per sampling day)
sd_time <- as(phyloseq::sample_data(ps_all), "data.frame")
time_levels_neutral <- sort(unique(as.numeric(sd_time$Time)))
time_levels_neutral <- time_levels_neutral[!is.na(time_levels_neutral)]
if (!isTRUE(neutral_model_all_times)) {
  time_levels_neutral <- target_time
}

message(
  "Neutral model time strata: ",
  paste(time_levels_neutral, collapse = ", "),
  if (isTRUE(neutral_model_all_times)) " (all)" else paste0(" (day ", target_time, " only)")
)

# Retain ASV-level rarefied matrix + metadata for Culicoidibacter bar plot (tv == target_time, ASV only)
ps_rare_asv <- NULL
asv_otu_mat <- NULL
asv_meta_df <- NULL
asv_fit_for_plot <- NULL

for (tv in time_levels_neutral) {
  ps_timed <- ps_all %>%
    microViz::ps_filter(Time == tv)

  if (phyloseq::nsamples(ps_timed) < 2L) {
    message("Skip Time ", tv, ": need >= 2 samples; got ", phyloseq::nsamples(ps_timed))
    next
  }

  message("Neutral model: Time == ", tv, " d; n samples = ", phyloseq::nsamples(ps_timed))

  for (run_tag in names(tax_rank_runs)) {
  tax_rank <- tax_rank_runs[[run_tag]]

  ps_work <- ps_timed
  if (!is.null(tax_rank)) {
    ps_work <- microViz::tax_agg(ps_work, rank = tax_rank)
  }

  set.seed(42L)
  ps_rare <- phyloseq::rarefy_even_depth(ps_work, rngseed = 42L, replace = FALSE, trimOTUs = TRUE, verbose = FALSE)

  meta_df <- microViz::samdat_tbl(ps_rare) %>%
    dplyr::mutate(.sample_name = as.character(.sample_name))

  otu_mat <- as(phyloseq::otu_table(ps_rare), "matrix")
  if (phyloseq::taxa_are_rows(ps_rare)) {
    otu_mat <- t(otu_mat)
  }
  colnames(otu_mat) <- phyloseq::taxa_names(ps_rare)
  sample_ids <- phyloseq::sample_names(ps_rare)
  otu_mat <- otu_mat[sample_ids, , drop = FALSE]
  meta_df <- meta_df %>% dplyr::filter(.sample_name %in% sample_ids)
  meta_df <- meta_df[match(sample_ids, meta_df$.sample_name), ]

  rs <- rowSums(otu_mat)
  if (stats::sd(rs) > 1e-6) {
    message("Note: row sums not perfectly equal (sd = ", stats::sd(rs), ") — rarefaction check.")
  }

  spp <- otu_mat
  fits_by_time <- list()
  boot_by_time <- list()
  predictions_all <- list()
  partitions_long <- list()
  partition_equalized <- list()

  if (nrow(spp) < 2L) {
    message("Skipping ", run_tag, ": fewer than 2 samples.")
    next
  }

  fit <- sncm_fit_burns(spp)
  fits_by_time[[as.character(tv)]] <- fit

  if (isTRUE(fit$success)) {
    pred <- fit$predictions %>%
      dplyr::mutate(Time = tv) %>%
      annotate_predictions_highlight(ps_rare, ., tax_rank, highlight_genus_pattern)
    predictions_all[[1L]] <- pred

    if (run_tag == "ASV" && tv == target_time) {
      ps_rare_asv <- ps_rare
      asv_otu_mat <- otu_mat
      asv_meta_df <- meta_df
      asv_fit_for_plot <- fit
    }

    part_tab <- pred %>%
      dplyr::count(partition, name = "n_taxa")
    partitions_long[[1L]] <- part_tab %>% dplyr::mutate(Time = tv)

    if (isTRUE(equalize_partition_sizes)) {
      by_part <- split(pred$taxon_id, pred$partition)
      partition_equalized[[as.character(tv)]] <- equalize_partition_taxa(by_part, seed = 42L)
    }

    focal_list <- if (exists("focal_genera", inherits = TRUE) && length(focal_genera) > 0L) {
      as.character(focal_genera)
    } else {
      as.character(highlight_genus_pattern)
    }
    focal_slug <- function(x) {
      x <- as.character(x)
      x <- gsub("[^A-Za-z0-9]+", "_", x, perl = TRUE)
      x <- gsub("^_+|_+$", "", x, perl = TRUE)
      tolower(x)
    }
    for (g in unique(focal_list)) {
      pred_g <- fit$predictions %>%
        dplyr::mutate(Time = tv) %>%
        annotate_predictions_highlight(ps_rare, ., tax_rank, g)

      out_new <- file.path(path_tbl, paste0("focal_taxa__", focal_slug(g), "__Time", tv, "__", run_tag, ".csv"))
      readr::write_csv(pred_g %>% dplyr::filter(.data$is_culici), out_new)

      # Backward-compatible legacy filename for Culicoidibacter
      if (identical(as.character(g), "Culicoidibacter")) {
        out_legacy <- file.path(path_tbl, paste0("culicoidibacter_taxa__Time", tv, "__", run_tag, ".csv"))
        readr::write_csv(pred_g %>% dplyr::filter(.data$is_culici), out_legacy)
      }
    }
  } else {
    message("Run ", run_tag, ": fit failed — ", fit$error)
  }

  message("Bootstrap (", run_tag, ") Time = ", tv, " ...")
  boot_by_time[[as.character(tv)]] <- bootstrap_sncm_stats(spp, n_boot = n_bootstrap)

  fit_stats_combined <- purrr::imap_dfr(
    fits_by_time,
    function(fit, nm) {
      if (!isTRUE(fit$success) || is.null(fit$fit_stats)) {
        return(tibble::tibble(Time = as.numeric(nm)))
      }
      fit$fit_stats %>% dplyr::mutate(Time = as.numeric(nm))
    }
  )

  predictions_bound <- dplyr::bind_rows(predictions_all)
  partitions_summary <- dplyr::bind_rows(partitions_long)

  boot_summary <- purrr::imap_dfr(
    boot_by_time,
    function(bt, nm) {
      if (nrow(bt) < 1L) {
        return(tibble::tibble(Time = as.numeric(nm)))
      }
      bt %>%
        dplyr::summarise(
          m_nls_q025 = stats::quantile(.data$m_nls, 0.025, na.rm = TRUE),
          m_nls_q975 = stats::quantile(.data$m_nls, 0.975, na.rm = TRUE),
          rsqr_q025 = stats::quantile(.data$rsqr, 0.025, na.rm = TRUE),
          rsqr_q975 = stats::quantile(.data$rsqr, 0.975, na.rm = TRUE)
        ) %>%
        dplyr::mutate(Time = as.numeric(nm))
    }
  )

  readr::write_csv(fit_stats_combined, file.path(path_tbl, paste0("fit_stats__Time", tv, "__", run_tag, ".csv")))
  readr::write_csv(predictions_bound, file.path(path_tbl, paste0("predictions_per_taxon__Time", tv, "__", run_tag, ".csv")))
  readr::write_csv(partitions_summary, file.path(path_tbl, paste0("partition_counts__Time", tv, "__", run_tag, ".csv")))
  readr::write_csv(boot_summary, file.path(path_tbl, paste0("bootstrap_quantiles__Time", tv, "__", run_tag, ".csv")))

  results_bundle <- list(
    tax_rank = tax_rank,
    run_tag = run_tag,
    target_time = target_time,
    time_stratum = tv,
    neutral_time_levels = time_levels_neutral,
    fits_by_time = fits_by_time,
    boot_by_time = boot_by_time,
    fit_stats = fit_stats_combined,
    predictions = predictions_bound,
    partitions_summary = partitions_summary,
    boot_summary = boot_summary,
    partition_equalized = partition_equalized,
    ps_rare_sample_data = meta_df,
    highlight_pattern = highlight_genus_pattern
  )
  saveRDS(results_bundle, file.path(path_rds, paste0("neutral_model_results__Time", tv, "__", run_tag, ".rds")))

  # --- Figures ---
  fit <- fits_by_time[[as.character(tv)]]
  if (isTRUE(fit$success) && !is.null(fit$predictions)) {
    pred <- predictions_bound
    N <- fit$fit_stats$N[[1]]
    m_est <- fit$fit_stats$m_nls[[1]]
    d <- fit$fit_stats$detect_d[[1]]

    p_grid <- 10^seq(log10(max(min(pred$p), 1e-8)), 0, length.out = 200)
    curve_df <- tibble::tibble(
      p = p_grid,
      freq_pred = stats::pbeta(d, N * m_est * p_grid, N * m_est * (1 - p_grid), lower.tail = FALSE)
    )
    n_samples_curve <- fit$fit_stats$n_samples[[1]]
    curve_ci_wilson <- Hmisc::binconf(
      curve_df$freq_pred * n_samples_curve,
      n_samples_curve,
      alpha = 0.05,
      method = "wilson",
      return.df = TRUE
    )
    curve_df <- curve_df %>%
      dplyr::mutate(
        freq_lwr = curve_ci_wilson$Lower,
        freq_upr = curve_ci_wilson$Upper
      )

    part_cols <- partition_colors

    p_main <- plot_neutral_fitted(
      pred,
      curve_df,
      tv,
      m_est,
      run_tag,
      highlight = FALSE,
      highlight_label = highlight_genus_pattern,
      partition_colors = part_cols
    )
    ggplot2::ggsave(
      filename = file.path(path_fig, paste0("neutral_model_Time", tv, "__", run_tag, ".png")),
      plot = p_main,
      width = 7,
      height = 5,
      dpi = 300L
    )

    rsqr_est <- fit$fit_stats$rsqr[[1]]
    focal_list <- if (exists("focal_genera", inherits = TRUE) && length(focal_genera) > 0L) {
      as.character(focal_genera)
    } else {
      as.character(highlight_genus_pattern)
    }
    focal_slug <- function(x) {
      x <- as.character(x)
      x <- gsub("[^A-Za-z0-9]+", "_", x, perl = TRUE)
      x <- gsub("^_+|_+$", "", x, perl = TRUE)
      tolower(x)
    }
    for (g in unique(focal_list)) {
      pred_g <- annotate_predictions_highlight(ps_rare, pred, tax_rank, g)
      p_hi <- plot_neutral_fitted(
        pred_g,
        curve_df,
        tv,
        m_est,
        run_tag,
        highlight = TRUE,
        highlight_label = g,
        partition_colors = part_cols,
        rsqr = rsqr_est
      )
      out_path <- file.path(path_fig, paste0("neutral_model_Time", tv, "__", run_tag, "__", focal_slug(g), ".png"))
      ggplot2::ggsave(filename = out_path, plot = p_hi, width = 8, height = 6.5, dpi = 300L)

      # Legacy filename retained for Culicoidibacter
      if (identical(as.character(g), "Culicoidibacter")) {
        ggplot2::ggsave(
          filename = file.path(path_fig, paste0("neutral_model_Time", tv, "__", run_tag, "__Culicoidibacter.png")),
          plot = p_hi,
          width = 8,
          height = 6.5,
          dpi = 300L
        )
      }
    }

    # One panel: all focal genera on the same Genus-level neutral-model axes
    if (identical(run_tag, "Genus")) {
      pred_multi <- pred %>%
        dplyr::mutate(
          genus_label = as.character(.data$taxon_id),
          is_focal = tolower(as.character(.data$taxon_id)) %in%
            tolower(as.character(unique(focal_list)))
        )
      if (sum(pred_multi$is_focal, na.rm = TRUE) >= 1L) {
        p_multi <- plot_neutral_fitted_multi_focal(
          pred_multi,
          curve_df,
          tv,
          m_est,
          run_tag,
          partition_colors = part_cols,
          rsqr = rsqr_est
        )
        p_neutral_focal_four_genus <- p_multi
        stem_multi <- paste0("neutral_model_Time", tv, "__", run_tag, "__focal_four_genera")
        ggplot2::ggsave(
          filename = file.path(path_fig, paste0(stem_multi, ".png")),
          plot = p_multi,
          width = 7,
          height = 7,
          dpi = 300L
        )
        ggplot2::ggsave(
          filename = file.path(path_fig, paste0(stem_multi, ".pdf")),
          plot = p_multi,
          width = 7,
          height = 7,
          device = "pdf"
        )
      }
    }
  }

  message("Finished run: ", run_tag, " (Time ", tv, ")")
  }
}

# --- Culicoidibacter ASVs: focal neutral plot + pair proportion by Treatment / History ----------
if (!is.null(ps_rare_asv) && !is.null(asv_fit_for_plot) && isTRUE(asv_fit_for_plot$success)) {
  # Additional focal genera (reviewer appendix): write the same “two ASVs” focal panel when possible.
  if (exists("focal_genera", inherits = TRUE) && length(focal_genera) > 0L) {
    focal_slug <- function(x) {
      x <- as.character(x)
      x <- gsub("[^A-Za-z0-9]+", "_", x, perl = TRUE)
      x <- gsub("^_+|_+$", "", x, perl = TRUE)
      tolower(x)
    }
    other_genera <- setdiff(as.character(focal_genera), "Culicoidibacter")
    for (g in other_genera) {
      focal_path <- file.path(path_tbl, paste0("focal_taxa__", focal_slug(g), "__Time", target_time, "__ASV.csv"))
      if (!file.exists(focal_path)) {
        next
      }
      focal_tbl <- readr::read_csv(focal_path, show_col_types = FALSE)
      pair <- focal_tbl %>%
        dplyr::filter(.data$partition %in% c("above", "below")) %>%
        dplyr::distinct(.data$taxon_id, .data$partition)

      taxon_above <- pair$taxon_id[pair$partition == "above"][[1]]
      taxon_below <- pair$taxon_id[pair$partition == "below"][[1]]
      if (is.null(taxon_above) || is.null(taxon_below) || !nzchar(taxon_above) || !nzchar(taxon_below)) {
        next
      }

      pred_asv <- asv_fit_for_plot$predictions %>%
        dplyr::mutate(Time = target_time) %>%
        annotate_predictions_highlight(ps_rare_asv, ., NULL, g)

      N <- asv_fit_for_plot$fit_stats$N[[1]]
      m_est <- asv_fit_for_plot$fit_stats$m_nls[[1]]
      d <- asv_fit_for_plot$fit_stats$detect_d[[1]]
      n_samples_asv <- asv_fit_for_plot$fit_stats$n_samples[[1]]
      p_grid <- 10^seq(log10(max(min(pred_asv$p), 1e-8)), 0, length.out = 200)
      curve_df_asv <- tibble::tibble(
        p = p_grid,
        freq_pred = stats::pbeta(d, N * m_est * p_grid, N * m_est * (1 - p_grid), lower.tail = FALSE)
      )
      curve_ci_asv <- Hmisc::binconf(
        curve_df_asv$freq_pred * n_samples_asv,
        n_samples_asv,
        alpha = 0.05,
        method = "wilson",
        return.df = TRUE
      )
      curve_df_asv <- curve_df_asv %>%
        dplyr::mutate(
          freq_lwr = curve_ci_asv$Lower,
          freq_upr = curve_ci_asv$Upper
        )

      p_focal <- plot_neutral_focal_two_asvs(
        pred_asv,
        curve_df_asv,
        target_time,
        m_est,
        taxon_above,
        taxon_below,
        partition_colors = partition_colors
      )
      ggplot2::ggsave(
        filename = file.path(path_fig, paste0("neutral_model_Time", target_time, "__ASV__focal_two_", focal_slug(g), ".png")),
        plot = p_focal,
        width = 8,
        height = 6.5,
        dpi = 300L
      )
    }
  }

  cul_path <- file.path(path_tbl, paste0("culicoidibacter_taxa__Time", target_time, "__ASV.csv"))
  if (file.exists(cul_path)) {
    cul_asv_tbl <- readr::read_csv(cul_path, show_col_types = FALSE)
    pair <- cul_asv_tbl %>%
      dplyr::filter(.data$partition %in% c("above", "below")) %>%
      dplyr::distinct(.data$taxon_id, .data$partition)
    taxon_above <- pair$taxon_id[pair$partition == "above"][[1]]
    taxon_below <- pair$taxon_id[pair$partition == "below"][[1]]

    pred_asv <- asv_fit_for_plot$predictions %>%
      dplyr::mutate(Time = target_time) %>%
      annotate_predictions_highlight(ps_rare_asv, ., NULL, highlight_genus_pattern)

    N <- asv_fit_for_plot$fit_stats$N[[1]]
    m_est <- asv_fit_for_plot$fit_stats$m_nls[[1]]
    d <- asv_fit_for_plot$fit_stats$detect_d[[1]]
    n_samples_asv <- asv_fit_for_plot$fit_stats$n_samples[[1]]
    p_grid <- 10^seq(log10(max(min(pred_asv$p), 1e-8)), 0, length.out = 200)
    curve_df_asv <- tibble::tibble(
      p = p_grid,
      freq_pred = stats::pbeta(d, N * m_est * p_grid, N * m_est * (1 - p_grid), lower.tail = FALSE)
    )
    curve_ci_asv <- Hmisc::binconf(
      curve_df_asv$freq_pred * n_samples_asv,
      n_samples_asv,
      alpha = 0.05,
      method = "wilson",
      return.df = TRUE
    )
    curve_df_asv <- curve_df_asv %>%
      dplyr::mutate(
        freq_lwr = curve_ci_asv$Lower,
        freq_upr = curve_ci_asv$Upper
      )

    p_focal <- plot_neutral_focal_two_asvs(
      pred_asv,
      curve_df_asv,
      target_time,
      m_est,
      taxon_above,
      taxon_below,
      partition_colors = partition_colors
    )
    ggplot2::ggsave(
      filename = file.path(path_fig, paste0("neutral_model_Time", target_time, "__ASV__focal_two_Culici.png")),
      plot = p_focal,
      width = 8,
      height = 6.5,
      dpi = 300L
    )

    if (taxon_above %in% colnames(asv_otu_mat) && taxon_below %in% colnames(asv_otu_mat)) {
      sample_ids <- rownames(asv_otu_mat)
      props <- tibble::tibble(
        .sample_name = sample_ids,
        asv_above = asv_otu_mat[, taxon_above],
        asv_below = asv_otu_mat[, taxon_below]
      ) %>%
        dplyr::mutate(
          denom = .data$asv_above + .data$asv_below,
          prop_of_above_in_pair = dplyr::if_else(
            .data$denom > 0,
            .data$asv_above / .data$denom,
            NA_real_
          )
        )

      meta_b <- asv_meta_df %>% dplyr::mutate(.sample_name = as.character(.sample_name))

      if (!"History_Label" %in% names(meta_b)) {
        meta_b <- meta_b %>%
          dplyr::mutate(
            History_Label = dplyr::case_when(
              as.character(.data$Treatment) %in% c("A- T- P-", "A- T- P+") ~ "No prior stressors",
              as.character(.data$Treatment) %in% c(
                "A+ T- P-", "A+ T- P+", "A- T+ P-", "A- T+ P+"
              ) ~ "One prior stressor",
              as.character(.data$Treatment) %in% c("A+ T+ P-", "A+ T+ P+") ~ "Two prior stressors",
              TRUE ~ NA_character_
            ),
            History_Label = factor(
              .data$History_Label,
              levels = c("No prior stressors", "One prior stressor", "Two prior stressors")
            )
          )
      }

      plot_bar_df <- meta_b %>%
        dplyr::left_join(props, by = ".sample_name")

      if (exists("treatment_order", inherits = TRUE)) {
        plot_bar_df <- plot_bar_df %>%
          dplyr::mutate(Treatment = factor(as.character(.data$Treatment), levels = treatment_order))
      }

      agg_bar <- plot_bar_df %>%
        dplyr::group_by(.data$Treatment, .data$History_Label) %>%
        dplyr::summarise(
          mean_prop_above = mean(.data$prop_of_above_in_pair, na.rm = TRUE),
          n = dplyr::n(),
          n_with_pair = sum(!is.na(.data$prop_of_above_in_pair)),
          .groups = "drop"
        ) %>%
        dplyr::mutate(mean_prop_below = 1 - .data$mean_prop_above)

      readr::write_csv(
        agg_bar,
        file.path(path_tbl, paste0("culicoidibacter_two_asv_pair_props__Time", target_time, ".csv"))
      )

      lab_above <- paste0(taxon_above, " (above-neutral)")
      lab_below <- paste0(taxon_below, " (below-neutral)")
      al <- agg_bar %>%
        tidyr::pivot_longer(
          cols = c(mean_prop_above, mean_prop_below),
          names_to = "which",
          values_to = "prop"
        ) %>%
        dplyr::mutate(
          asv_label = dplyr::if_else(.data$which == "mean_prop_above", lab_above, lab_below),
          # Stack order: first level = bottom (below = steelblue), last = top (above = firebrick)
          asv_label = factor(.data$asv_label, levels = c(lab_below, lab_above))
        )

      fill_vec <- stats::setNames(c("steelblue", "firebrick"), c(lab_below, lab_above))

      # Total rarefied reads per focal ASV (all samples at this time point; exposure regimes pooled)
      total_reads_above <- sum(asv_otu_mat[, taxon_above])
      total_reads_below <- sum(asv_otu_mat[, taxon_below])
      fmt_reads <- function(x) format(as.integer(round(x)), big.mark = ",", trim = TRUE)

      p_bar <- ggplot2::ggplot(al, ggplot2::aes(x = .data$Treatment, y = .data$prop, fill = .data$asv_label)) +
        ggplot2::geom_col(position = "stack", width = 0.85, color = "white", linewidth = 0.2) +
        ggplot2::facet_wrap(~History_Label, nrow = 1L, scales = "free_x") +
        ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
        ggplot2::scale_fill_manual(values = fill_vec, name = "ASV (neutral partition)") +
        ggplot2::labs(
          title = "Culicoidibacter ASVs: mean within-pair proportion by treatment and prior stress history",
          subtitle = paste0("Day ", target_time),
          caption = paste0(
            "Proportion = ", taxon_above, "/(", taxon_above, "+", taxon_below, ") per sample (rarefied).\n",
            "Bars: means by treatment × prior-stress history.\n",
            "Total rarefied reads (all samples, regimes pooled): ", taxon_above, " = ",
            fmt_reads(total_reads_above), "; ", taxon_below, " = ", fmt_reads(total_reads_below), "."
          ),
          x = "Treatment (exposure regime)",
          y = "Mean proportion (within ASV pair)"
        ) +
        theme_sieler2026_publication(base_size = 14, legend_position = "bottom") +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.box = "horizontal",
          plot.caption = ggplot2::element_text(size = 8.5, hjust = 0, colour = "grey25", lineheight = 1.25),
          plot.caption.position = "plot",
          plot.margin = ggplot2::margin(t = 5.5, r = 8, b = 18, l = 8, unit = "pt")
        ) +
        manuscript_facet_strip_theme()

      ggplot2::ggsave(
        filename = file.path(path_fig, paste0("culicoidibacter_two_asv_pair_by_treatment_History__Time", target_time, ".png")),
        plot = p_bar,
        width = 11,
        height = 5.95,
        dpi = 300L
      )

      # --- Three-way stack: above + below + other Culicoidibacter ASVs (neutral partition) ----------
      neutral_taxa <- pred_asv %>%
        dplyr::filter(.data$partition == "neutral", .data$is_culici %in% TRUE) %>%
        dplyr::pull(.data$taxon_id) %>%
        unique()
      neutral_taxa <- intersect(neutral_taxa, colnames(asv_otu_mat))

      other_counts <- if (length(neutral_taxa) > 0L) {
        rowSums(asv_otu_mat[, neutral_taxa, drop = FALSE])
      } else {
        rep(0, nrow(asv_otu_mat))
      }
      total_reads_other <- sum(other_counts)

      props3 <- tibble::tibble(
        .sample_name = sample_ids,
        asv_above = asv_otu_mat[, taxon_above],
        asv_below = asv_otu_mat[, taxon_below],
        other_culici = other_counts
      ) %>%
        dplyr::mutate(
          total_culici = .data$asv_above + .data$asv_below + .data$other_culici,
          prop_above = dplyr::if_else(.data$total_culici > 0, .data$asv_above / .data$total_culici, NA_real_),
          prop_below = dplyr::if_else(.data$total_culici > 0, .data$asv_below / .data$total_culici, NA_real_),
          prop_other = dplyr::if_else(.data$total_culici > 0, .data$other_culici / .data$total_culici, NA_real_)
        )

      plot_bar_df3 <- meta_b %>%
        dplyr::left_join(props3, by = ".sample_name")

      if (exists("treatment_order", inherits = TRUE)) {
        plot_bar_df3 <- plot_bar_df3 %>%
          dplyr::mutate(Treatment = factor(as.character(.data$Treatment), levels = treatment_order))
      }

      agg_bar3 <- plot_bar_df3 %>%
        dplyr::group_by(.data$Treatment, .data$History_Label) %>%
        dplyr::summarise(
          mean_prop_above = mean(.data$prop_above, na.rm = TRUE),
          mean_prop_below = mean(.data$prop_below, na.rm = TRUE),
          mean_prop_other = mean(.data$prop_other, na.rm = TRUE),
          n = dplyr::n(),
          n_with_culici = sum(.data$total_culici > 0, na.rm = TRUE),
          .groups = "drop"
        )

      readr::write_csv(
        agg_bar3,
        file.path(path_tbl, paste0("culicoidibacter_three_way_props__Time", target_time, ".csv"))
      )

      lab_other <- "Other Culicoidibacter (neutral expectation)"
      lab_a3 <- paste0(taxon_above, " (above-neutral)")
      lab_b3 <- paste0(taxon_below, " (below-neutral)")

      al3 <- agg_bar3 %>%
        tidyr::pivot_longer(
          cols = c(mean_prop_above, mean_prop_below, mean_prop_other),
          names_to = "which",
          values_to = "prop"
        ) %>%
        dplyr::mutate(
          asv_label = dplyr::case_when(
            .data$which == "mean_prop_above" ~ lab_a3,
            .data$which == "mean_prop_below" ~ lab_b3,
            TRUE ~ lab_other
          ),
          # Stack order: bottom = below (steelblue), middle = other neutral (grey35), top = above (firebrick)
          asv_label = factor(.data$asv_label, levels = c(lab_b3, lab_other, lab_a3))
        )

      fill_vec3 <- stats::setNames(
        c("steelblue", "grey35", "firebrick"),
        c(lab_b3, lab_other, lab_a3)
      )

      p_bar3 <- ggplot2::ggplot(al3, ggplot2::aes(x = .data$Treatment, y = .data$prop, fill = .data$asv_label)) +
        ggplot2::geom_col(position = "stack", width = 0.85, color = "white", linewidth = 0.2) +
        ggplot2::facet_wrap(~History_Label, nrow = 1L, scales = "free_x") +
        ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
        ggplot2::scale_fill_manual(values = fill_vec3, name = "Culicoidibacter ASVs") +
        ggplot2::labs(
          title = "Culicoidibacter: proportion of rarefied reads (above vs below vs other neutral ASVs)",
          subtitle = paste0("Day ", target_time),
          caption = paste0(
            "Per sample: category / total Culicoidibacter reads (rarefied).\n",
            "Bars: means by treatment × prior-stress history.\n",
            "Total rarefied reads (all samples, regimes pooled): ", taxon_above, " = ",
            fmt_reads(total_reads_above), "; ", taxon_below, " = ", fmt_reads(total_reads_below), ".\n",
            "Other neutral ASVs = ", fmt_reads(total_reads_other), "; total = ",
            fmt_reads(total_reads_above + total_reads_below + total_reads_other), "."
          ),
          x = "Treatment (exposure regime)",
          y = "Mean proportion (of Culicoidibacter reads)"
        ) +
        theme_sieler2026_publication(base_size = 14, legend_position = "bottom") +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.box = "horizontal",
          plot.caption = ggplot2::element_text(size = 8.5, hjust = 0, colour = "grey25", lineheight = 1.25),
          plot.caption.position = "plot",
          plot.margin = ggplot2::margin(t = 5.5, r = 8, b = 22, l = 8, unit = "pt")
        ) +
        manuscript_facet_strip_theme()

      ggplot2::ggsave(
        filename = file.path(path_fig, paste0("culicoidibacter_three_way_by_treatment_History__Time", target_time, ".png")),
        plot = p_bar3,
        width = 11,
        height = 6.15,
        dpi = 300L
      )
    }
  }
}

message("Neutral model analysis complete. Outputs -> ", out_dir)

# =============================================================================
# Inlined: 08__NeutralModel__byregime_partitions.R
# =============================================================================

# Genus-level Sloan predictions at `target_time` (same stem as main loop `write_csv`).
pred_path <- file.path(path_tbl, paste0("predictions_per_taxon__Time", target_time, "__Genus.csv"))

if (!file.exists(pred_path)) {
  stop("Run Code/01__Analysis/08__NeutralModel.R first to create ", pred_path)
}

pred <- readr::read_csv(pred_path, show_col_types = FALSE) %>%
  dplyr::distinct(.data$taxon_id, .data$partition) %>%
  dplyr::filter(.data$partition %in% c("above", "neutral", "below"))

if (!exists("ps.list", inherits = TRUE) || is.null(ps.list) || !"All" %in% names(ps.list)) {
  stop("ps.list not found. Run 04__DataPreProcess.R.")
}

ps_all <- ps.list[["All"]]
ps_t <- ps_all %>%
  microViz::ps_filter(as.numeric(.data$Time) == target_time) %>%
  # Must match Genus-level `predictions_per_taxon__Time*__Genus.csv` from the main neutral loop.
  microViz::tax_agg(rank = "Genus")

if (phyloseq::nsamples(ps_t) < 2L) {
  stop("Need >= 2 samples at Time == ", target_time)
}

set.seed(42L)
ps_rare <- phyloseq::rarefy_even_depth(ps_t, rngseed = 42L, replace = FALSE, trimOTUs = TRUE, verbose = FALSE)

otu_mat <- as(phyloseq::otu_table(ps_rare), "matrix")
if (phyloseq::taxa_are_rows(ps_rare)) {
  otu_mat <- t(otu_mat)
}
colnames(otu_mat) <- phyloseq::taxa_names(ps_rare)
sample_ids <- phyloseq::sample_names(ps_rare)
otu_mat <- otu_mat[sample_ids, , drop = FALSE]

meta <- microViz::samdat_tbl(ps_rare) %>%
  dplyr::mutate(
    .sample_name = as.character(.data$.sample_name),
    History_Label = dplyr::case_when(
      as.character(.data$Treatment) %in% c("A- T- P-", "A- T- P+") ~ "No prior stressors",
      as.character(.data$Treatment) %in% c(
        "A+ T- P-", "A+ T- P+", "A- T+ P-", "A- T+ P+"
      ) ~ "One prior stressor",
      as.character(.data$Treatment) %in% c("A+ T+ P-", "A+ T+ P+") ~ "Two prior stressors",
      TRUE ~ NA_character_
    ),
    History_Label = factor(
      .data$History_Label,
      levels = c("No prior stressors", "One prior stressor", "Two prior stressors")
    )
  )

if (exists("treatment_order", inherits = TRUE)) {
  meta <- meta %>%
    dplyr::mutate(Treatment = factor(as.character(.data$Treatment), levels = treatment_order))
}

taxa_fit <- intersect(colnames(otu_mat), pred$taxon_id)
if (length(taxa_fit) < 3L) {
  stop("Too few overlapping taxa between phyloseq and predictions.")
}

pred_fit <- pred %>% dplyr::filter(.data$taxon_id %in% taxa_fit)

# Long counts: fitted taxa only
long <- tibble::as_tibble(otu_mat[, taxa_fit, drop = FALSE], rownames = ".sample_name") %>%
  tidyr::pivot_longer(
    cols = -".sample_name",
    names_to = "taxon_id",
    values_to = "count"
  ) %>%
  dplyr::inner_join(pred_fit, by = "taxon_id")

if (!".sample_name" %in% names(meta)) {
  stop("sample_data must include `.sample_name` (microViz::samdat_tbl).")
}

meta_b <- meta %>%
  dplyr::mutate(.sample_name = as.character(.data$.sample_name))

long <- long %>%
  dplyr::inner_join(
    meta_b %>% dplyr::select(".sample_name", "Treatment", "History_Label"),
    by = ".sample_name"
  )

# --- Taxon-based: among fitted ASVs detected (>0) in at least one sample in the stratum ------
taxa_detected_regime <- long %>%
  dplyr::filter(.data$count > 0) %>%
  dplyr::distinct(.data$Treatment, .data$taxon_id, .data$partition)

tab_taxa_regime <- taxa_detected_regime %>%
  dplyr::count(.data$Treatment, .data$partition, name = "n_taxa") %>%
  tidyr::pivot_wider(
    names_from = "partition",
    values_from = "n_taxa",
    values_fill = 0L
  )

for (nm in c("above", "neutral", "below")) {
  if (!nm %in% names(tab_taxa_regime)) {
    tab_taxa_regime[[nm]] <- 0L
  }
}

tab_taxa_regime <- tab_taxa_regime %>%
  dplyr::group_by(.data$Treatment) %>%
  dplyr::mutate(n_total = sum(.data$above + .data$neutral + .data$below, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    dplyr::across(dplyr::all_of(c("above", "neutral", "below")), as.integer),
    pct_above = 100 * .data$above / .data$n_total,
    pct_neutral = 100 * .data$neutral / .data$n_total,
    pct_below = 100 * .data$below / .data$n_total
  ) %>%
  dplyr::left_join(
    meta_b %>%
      dplyr::distinct(.data$Treatment, .data$History_Label),
    by = "Treatment"
  ) %>%
  dplyr::arrange(.data$Treatment)

readr::write_csv(
  tab_taxa_regime,
  file.path(path_tbl, paste0("neutral_partition__by_regime__taxon_counts__Time", target_time, "__ASV.csv"))
)

# --- Read-based: per sample, proportion of rarefied reads in fitted taxa by partition ---------
props_sample <- long %>%
  dplyr::group_by(.data$.sample_name, .data$Treatment, .data$History_Label, .data$partition) %>%
  dplyr::summarise(reads = sum(.data$count), .groups = "drop") %>%
  dplyr::group_by(.data$.sample_name) %>%
  dplyr::mutate(
    total_fit = sum(.data$reads),
    prop = dplyr::if_else(.data$total_fit > 0, .data$reads / .data$total_fit, NA_real_)
  ) %>%
  dplyr::ungroup()

mean_prop_regime <- props_sample %>%
  dplyr::group_by(.data$Treatment, .data$History_Label, .data$partition) %>%
  dplyr::summarise(mean_prop_within_fit = mean(.data$prop, na.rm = TRUE), .groups = "drop")

wide_reads <- mean_prop_regime %>%
  tidyr::pivot_wider(
    names_from = "partition",
    values_from = "mean_prop_within_fit",
    values_fill = 0
  )

for (nm in c("above", "neutral", "below")) {
  if (!nm %in% names(wide_reads)) {
    wide_reads[[nm]] <- 0
  }
}

readr::write_csv(
  wide_reads %>%
    dplyr::arrange(.data$Treatment),
  file.path(path_tbl, paste0("neutral_partition__by_regime__read_mean_prop__Time", target_time, "__ASV.csv"))
)

# Pooled read mass per regime (descriptive)
pooled_regime <- long %>%
  dplyr::group_by(.data$Treatment, .data$partition) %>%
  dplyr::summarise(reads = sum(.data$count), .groups = "drop") %>%
  dplyr::group_by(.data$Treatment) %>%
  dplyr::mutate(pooled_prop = .data$reads / sum(.data$reads)) %>%
  dplyr::ungroup()

readr::write_csv(
  pooled_regime %>%
    tidyr::pivot_wider(names_from = "partition", values_from = c("reads", "pooled_prop"), values_fill = 0),
  file.path(path_tbl, paste0("neutral_partition__by_regime__read_pooled__Time", target_time, "__ASV.csv"))
)

# --- By prior stressor history only (pooled regimes) -----------------------------------------
taxa_detected_hist <- long %>%
  dplyr::filter(.data$count > 0) %>%
  dplyr::distinct(.data$History_Label, .data$taxon_id, .data$partition)

tab_taxa_hist <- taxa_detected_hist %>%
  dplyr::count(.data$History_Label, .data$partition, name = "n_taxa") %>%
  tidyr::pivot_wider(names_from = "partition", values_from = "n_taxa", values_fill = 0L)

for (nm in c("above", "neutral", "below")) {
  if (!nm %in% names(tab_taxa_hist)) {
    tab_taxa_hist[[nm]] <- 0L
  }
}

tab_taxa_hist <- tab_taxa_hist %>%
  dplyr::mutate(
    dplyr::across(dplyr::all_of(c("above", "neutral", "below")), as.integer),
    n_total = .data$above + .data$neutral + .data$below
  )

readr::write_csv(
  tab_taxa_hist,
  file.path(path_tbl, paste0("neutral_partition__by_history__taxon_counts__Time", target_time, "__ASV.csv"))
)

mean_prop_hist <- props_sample %>%
  dplyr::group_by(.data$History_Label, .data$partition) %>%
  dplyr::summarise(mean_prop_within_fit = mean(.data$prop, na.rm = TRUE), .groups = "drop")

readr::write_csv(
  mean_prop_hist %>%
    tidyr::pivot_wider(names_from = "partition", values_from = "mean_prop_within_fit", values_fill = 0) %>%
    dplyr::arrange(.data$History_Label),
  file.path(path_tbl, paste0("neutral_partition__by_history__read_mean_prop__Time", target_time, "__ASV.csv"))
)

# --- ggplot: stacked bars (mean read proportion within fitted taxa) -------------------------
lab_part <- c(
  above = "Above neutral expectation",
  neutral = "Neutral",
  below = "Below neutral expectation"
)

plot_df <- mean_prop_regime %>%
  dplyr::mutate(
    partition = factor(.data$partition, levels = c("below", "neutral", "above")),
    partition_lab = lab_part[as.character(.data$partition)]
  ) %>%
  dplyr::mutate(
    partition_lab = factor(.data$partition_lab, levels = unname(lab_part[c("below", "neutral", "above")]))
  )

# Stack order (ggplot: first level = bottom of stack): below (blue) -> neutral (grey) -> above (red, top)
fill_vec <- stats::setNames(
  c(partition_colors[["below"]], partition_colors[["neutral"]], partition_colors[["above"]]),
  levels(plot_df$partition_lab)
)
# Legend order: top-of-stack (above) first, for congruence with the figure
fill_breaks_rev <- rev(levels(plot_df$partition_lab))

caption_txt <- paste0(
  "Partitions from a single global Sloan neutral model fit at day ", target_time,
  " (metacommunity = all samples). Y = mean sample-wise proportion of rarefied reads ",
  "assigned to ASVs in each partition, among ASVs included in the fit.\n",
  "Stratification by regime does not re-estimate neutral expectations per regime."
)

p_stack <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$Treatment, y = .data$mean_prop_within_fit, fill = .data$partition_lab)) +
  ggplot2::geom_col(width = 0.85, color = "white", linewidth = 0.2) +
  ggplot2::facet_wrap(~History_Label, nrow = 1L, scales = "free_x") +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::scale_fill_manual(
    values = fill_vec,
    breaks = fill_breaks_rev,
    name = "Neutral model partition"
  ) +
  ggplot2::labs(
    title = "ASV neutral-model partitions: mean read proportion by exposure regime",
    subtitle = paste0("Day ", target_time, "; among fitted ASVs only"),
    caption = caption_txt,
    x = "Treatment (exposure regime)",
    y = "Mean proportion of reads (fitted ASVs)"
  ) +
  theme_sieler2026_publication(base_size = 14, legend_position = "bottom") +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.direction = "horizontal",
    strip.background = ggplot2::element_rect(fill = "white", colour = "black", linewidth = 0.4),
    strip.text = ggplot2::element_text(face = "bold", colour = "black"),
    panel.spacing.x = grid::unit(1.1, "lines"),
    plot.caption = ggplot2::element_text(size = 8.5, hjust = 0, colour = "grey25", lineheight = 1.2),
    plot.caption.position = "plot",
    plot.margin = ggplot2::margin(t = 5.5, r = 8, b = 20, l = 8, unit = "pt")
  )

ggplot2::ggsave(
  filename = file.path(path_fig, paste0("neutral_partition_read_mean_prop_by_regime__Time", target_time, "__ASV.png")),
  plot = p_stack,
  width = 12,
  height = 6.2,
  dpi = 300L
)

# Taxon fraction stacked bar (percent of detected fitted taxa per partition)
plot_taxa <- tab_taxa_regime %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(c("above", "neutral", "below")),
    names_to = "partition",
    values_to = "n_taxa"
  ) %>%
  dplyr::mutate(
    prop_taxa = .data$n_taxa / .data$n_total,
    partition_lab = lab_part[.data$partition],
    partition_lab = factor(.data$partition_lab, levels = unname(lab_part[c("below", "neutral", "above")]))
  )

p_taxa <- ggplot2::ggplot(plot_taxa, ggplot2::aes(x = .data$Treatment, y = .data$prop_taxa, fill = .data$partition_lab)) +
  ggplot2::geom_col(width = 0.85, color = "white", linewidth = 0.2) +
  ggplot2::facet_wrap(~History_Label, nrow = 1L, scales = "free_x") +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::scale_fill_manual(
    values = fill_vec,
    breaks = fill_breaks_rev,
    name = "Neutral model partition"
  ) +
  ggplot2::labs(
    title = "ASV neutral-model partitions: fraction of detected fitted taxa",
    subtitle = paste0("Day ", target_time, "; taxa with >0 rarefied reads in regime; among fitted ASVs only"),
    caption = paste0(
      "Partitions from global day ", target_time, " fit. For each regime, fraction = ",
      "n taxa in partition / n detected fitted taxa in that regime. ",
      "Stratification describes detection within regime, not separate neutral expectations per regime."
    ),
    x = "Treatment (exposure regime)",
    y = "Fraction of detected taxa"
  ) +
  theme_sieler2026_publication(base_size = 14, legend_position = "bottom") +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.direction = "horizontal",
    strip.background = ggplot2::element_rect(fill = "white", colour = "black", linewidth = 0.4),
    strip.text = ggplot2::element_text(face = "bold", colour = "black"),
    panel.spacing.x = grid::unit(1.1, "lines"),
    plot.caption = ggplot2::element_text(size = 8.5, hjust = 0, colour = "grey25", lineheight = 1.2),
    plot.margin = ggplot2::margin(t = 5.5, r = 8, b = 18, l = 8, unit = "pt")
  )

ggplot2::ggsave(
  filename = file.path(path_fig, paste0("neutral_partition_taxa_fraction_by_regime__Time", target_time, "__ASV.png")),
  plot = p_taxa,
  width = 12,
  height = 6,
  dpi = 300L
)

# --- gt table: regime + history colors + counts ---------------------------------------------
gt_tbl <- tab_taxa_regime %>%
  dplyr::select(
    "Exposure regime" = "Treatment",
    "Prior stressor history" = "History_Label",
    "n above" = "above",
    "n neutral" = "neutral",
    "n below" = "below",
    "n total (detected)" = "n_total"
  ) %>%
  dplyr::mutate(
    `pct above` = round(100 * .data$`n above` / .data$`n total (detected)`, 1),
    `pct neutral` = round(100 * .data$`n neutral` / .data$`n total (detected)`, 1),
    `pct below` = round(100 * .data$`n below` / .data$`n total (detected)`, 1)
  )

gt_out <- gt::gt(gt_tbl) %>%
  gt::tab_header(
    title = "Neutral model partitions by exposure regime (global ASV fit, day 60)",
    subtitle = "Taxon counts among fitted ASVs with >0 rarefied reads in the regime"
  ) %>%
  gt::fmt_number(columns = c("n above", "n neutral", "n below", "n total (detected)"), decimals = 0) %>%
  gt::cols_label(
    `pct above` = "% above",
    `pct neutral` = "% neutral",
    `pct below` = "% below"
  )

if (exists("treatment_color_scale", inherits = TRUE)) {
  for (trt in unique(as.character(gt_tbl$`Exposure regime`))) {
    if (trt %in% names(treatment_color_scale)) {
      row_idx <- which(as.character(gt_tbl$`Exposure regime`) == trt)
      gt_out <- gt_out %>%
        gt::tab_style(
          style = gt::cell_fill(color = grDevices::adjustcolor(treatment_color_scale[[trt]], alpha.f = 0.35)),
          locations = gt::cells_body(columns = "Exposure regime", rows = row_idx)
        )
    }
  }
}

if (exists("history_color_scale", inherits = TRUE)) {
  for (h in unique(as.character(gt_tbl$`Prior stressor history`))) {
    if (!is.na(h) && h %in% names(history_color_scale)) {
      row_idx <- which(as.character(gt_tbl$`Prior stressor history`) == h)
      gt_out <- gt_out %>%
        gt::tab_style(
          style = gt::cell_fill(color = grDevices::adjustcolor(history_color_scale[[h]], alpha.f = 0.35)),
          locations = gt::cells_body(columns = "Prior stressor history", rows = row_idx)
        )
    }
  }
}

gt::gtsave(
  gt_out,
  filename = file.path(path_tbl, paste0("neutral_partition__by_regime_summary__Time", target_time, "__ASV.html"))
)

gt_hist <- tab_taxa_hist %>%
  dplyr::select(
    "Prior stressor history" = "History_Label",
    "n above" = "above",
    "n neutral" = "neutral",
    "n below" = "below",
    "n total (detected)" = "n_total"
  ) %>%
  dplyr::mutate(
    `pct above` = round(100 * .data$`n above` / .data$`n total (detected)`, 1),
    `pct neutral` = round(100 * .data$`n neutral` / .data$`n total (detected)`, 1),
    `pct below` = round(100 * .data$`n below` / .data$`n total (detected)`, 1)
  )

gt_hist_out <- gt::gt(gt_hist) %>%
  gt::tab_header(
    title = "Neutral model partitions pooled by prior stressor history (global ASV fit, day 60)",
    subtitle = "Taxon counts among fitted ASVs with >0 reads in any sample in that history stratum"
  ) %>%
  gt::fmt_number(columns = c("n above", "n neutral", "n below", "n total (detected)"), decimals = 0)

if (exists("history_color_scale", inherits = TRUE)) {
  for (h in unique(as.character(gt_hist$`Prior stressor history`))) {
    if (!is.na(h) && h %in% names(history_color_scale)) {
      row_idx <- which(as.character(gt_hist$`Prior stressor history`) == h)
      gt_hist_out <- gt_hist_out %>%
        gt::tab_style(
          style = gt::cell_fill(color = grDevices::adjustcolor(history_color_scale[[h]], alpha.f = 0.35)),
          locations = gt::cells_body(columns = "Prior stressor history", rows = row_idx)
        )
    }
  }
}

gt::gtsave(
  gt_hist_out,
  filename = file.path(path_tbl, paste0("neutral_partition__by_history_summary__Time", target_time, "__ASV.html"))
)

# --- Combined figure: community vs focal genus (same 3 partitions, two dodged stacks) -----------
focal_slug <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x, perl = TRUE)
  x <- gsub("^_+|_+$", "", x, perl = TRUE)
  tolower(x)
}

make_comm_vs_focal_plot <- function(focal_label, focal_taxa_path, out_png, out_long_csv, legacy_title = NULL) {
  if (!file.exists(focal_taxa_path)) {
    message("Skipping combined community vs ", focal_label, " figure: missing ", basename(focal_taxa_path), ".")
    return(invisible(FALSE))
  }

  focal_taxa <- readr::read_csv(focal_taxa_path, show_col_types = FALSE)
  if (!"is_culici" %in% names(focal_taxa)) {
    message("Skipping combined community vs ", focal_label, " figure: no is_culici column in ", basename(focal_taxa_path), ".")
    return(invisible(FALSE))
  }

  focal_ids <- focal_taxa %>%
    dplyr::filter(.data$is_culici %in% TRUE) %>%
    dplyr::pull(.data$taxon_id) %>%
    unique()
  focal_ids <- intersect(focal_ids, unique(long$taxon_id))

  if (length(focal_ids) < 1L) {
    message("Skipping combined community vs ", focal_label, " figure: no focal taxa overlap fitted counts.")
    return(invisible(FALSE))
  }

  # Per sample: partition shares among rarefied reads summed over all focal ASVs (global partition labels)
  props_focal_sample <- long %>%
    dplyr::filter(.data$taxon_id %in% focal_ids) %>%
    dplyr::group_by(.data$.sample_name, .data$Treatment, .data$History_Label, .data$partition) %>%
    dplyr::summarise(reads = sum(.data$count), .groups = "drop") %>%
    dplyr::group_by(.data$.sample_name) %>%
    dplyr::mutate(
      total_focal = sum(.data$reads),
      prop = dplyr::if_else(.data$total_focal > 0, .data$reads / .data$total_focal, NA_real_)
    ) %>%
    dplyr::ungroup()

  mean_prop_focal_regime <- props_focal_sample %>%
    dplyr::group_by(.data$Treatment, .data$History_Label, .data$partition) %>%
    dplyr::summarise(mean_prop = mean(.data$prop, na.rm = TRUE), .groups = "drop")

  regime_grid <- plot_df %>%
    dplyr::distinct(.data$Treatment, .data$History_Label) %>%
    tidyr::crossing(partition = c("below", "neutral", "above"))

  mean_prop_focal_regime <- regime_grid %>%
    dplyr::left_join(mean_prop_focal_regime, by = c("Treatment", "History_Label", "partition")) %>%
    dplyr::mutate(mean_prop = tidyr::replace_na(.data$mean_prop, 0))

  focal_long <- mean_prop_focal_regime %>%
    dplyr::mutate(
      partition_lab = lab_part[as.character(.data$partition)],
      partition_lab = factor(.data$partition_lab, levels = levels(plot_df$partition_lab)),
      bar_scope = paste0(focal_label, " (genus reads)"),
      prop = .data$mean_prop
    ) %>%
    dplyr::select("Treatment", "History_Label", "bar_scope", "partition_lab", "prop")

  comm_long <- plot_df %>%
    dplyr::mutate(
      bar_scope = "Community (fitted ASVs)",
      prop = .data$mean_prop_within_fit
    ) %>%
    dplyr::select("Treatment", "History_Label", "bar_scope", "partition_lab", "prop")

  plot_combined <- dplyr::bind_rows(comm_long, focal_long) %>%
    dplyr::mutate(
      bar_scope = factor(
        .data$bar_scope,
        levels = c("Community (fitted ASVs)", paste0(focal_label, " (genus reads)"))
      ),
      Bar_facet = dplyr::if_else(
        .data$bar_scope == "Community (fitted ASVs)",
        "All",
        "Focal"
      ),
      Bar_facet = factor(.data$Bar_facet, levels = c("All", "Focal")),
      History_Label = factor(
        as.character(.data$History_Label),
        levels = c("No prior stressors", "One prior stressor", "Two prior stressors")
      )
    )

  if (exists("treatment_order", inherits = TRUE)) {
    plot_combined <- plot_combined %>%
      dplyr::mutate(Treatment = factor(as.character(.data$Treatment), levels = treatment_order))
  }

  panel_levels <- plot_combined %>%
    dplyr::distinct(.data$Treatment, .data$History_Label) %>%
    dplyr::arrange(.data$History_Label, .data$Treatment)

  lev_panels <- paste(as.character(panel_levels$History_Label), as.character(panel_levels$Treatment), sep = "\n")

  plot_pair <- plot_combined %>%
    dplyr::mutate(
      panel_fac = paste(as.character(.data$History_Label), as.character(.data$Treatment), sep = "\n"),
      panel_fac = factor(.data$panel_fac, levels = lev_panels)
    )

  pal_pastel <- RColorBrewer::brewer.pal(9L, "Pastel1")
  vals_focal <- stats::setNames(
    c(pal_pastel[[2L]], pal_pastel[[9L]], pal_pastel[[1L]]),
    levels(plot_df$partition_lab)
  )
  part_leg_labs <- c(
    `Below neutral expectation` = "Below neutral",
    `Neutral` = "Neutral",
    `Above neutral expectation` = "Above neutral"
  )

  plot_pair_all <- plot_pair %>% dplyr::filter(.data$Bar_facet == "All")
  plot_pair_focal <- plot_pair %>% dplyr::filter(.data$Bar_facet == "Focal")

  readr::write_csv(plot_pair, out_long_csv)

  caption_combo <- paste0(
    "Sloan neutral partitions from a single global ASV fit (day ", target_time, "). ",
    "All: mean proportion across rarefied reads for all fitted ASVs; ",
    "Focal: same partition labels, proportions among ", focal_label, " reads only. ",
    "Stratification by regime does not refit the model."
  )

  if (!requireNamespace("ggnewscale", quietly = TRUE)) {
    stop("Install ggnewscale for the combined All vs focal figure (see Code/00__Setup/01__Libraries.R).")
  }

  p_combo <- ggplot2::ggplot() +
    ggplot2::geom_col(
      data = plot_pair_all,
      ggplot2::aes(x = .data$Bar_facet, y = .data$prop, fill = .data$partition_lab),
      width = 0.72,
      color = "black",
      linewidth = 0.28
    ) +
    ggplot2::scale_fill_manual(
      name = "All ASVs",
      values = fill_vec,
      breaks = fill_breaks_rev,
      labels = part_leg_labs[fill_breaks_rev],
      drop = FALSE
    ) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_col(
      data = plot_pair_focal,
      ggplot2::aes(x = .data$Bar_facet, y = .data$prop, fill = .data$partition_lab),
      width = 0.72,
      color = "black",
      linewidth = 0.28
    ) +
    ggplot2::scale_fill_manual(
      name = paste0(focal_label, "-associated ASVs"),
      values = vals_focal,
      breaks = fill_breaks_rev,
      labels = part_leg_labs[fill_breaks_rev],
      drop = FALSE
    ) +
    ggplot2::facet_wrap(~panel_fac, nrow = 1L, scales = "fixed") +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::labs(
      title = if (is.null(legacy_title)) paste0("Neutral-model partitions: All vs ", focal_label, " within each regime") else legacy_title,
      subtitle = paste0(
        "Day ", target_time,
        "; one row: prior stressors on top, exposure regime below; stack bottom-to-top = below / neutral / above"
      ),
      caption = caption_combo,
      x = NULL,
      y = "Mean proportion"
    ) +
    theme_sieler2026_publication(base_size = 14, legend_position = "bottom") +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 12, face = "plain"),
      axis.title.y = ggplot2::element_text(size = 12, face = "plain"),
      axis.text.x = ggplot2::element_text(size = 11, angle = 0, hjust = 0.5),
      axis.text.y = ggplot2::element_text(size = 11),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "vertical",
      legend.spacing.y = grid::unit(0.25, "cm"),
      strip.background = ggplot2::element_rect(fill = "white", colour = "black", linewidth = 0.4),
      strip.text = ggplot2::element_text(face = "bold", colour = "black", size = 8),
      panel.spacing = grid::unit(0.65, "lines"),
      plot.caption = ggplot2::element_text(size = 8.5, hjust = 0, colour = "grey25", lineheight = 1.25),
      plot.caption.position = "plot",
      plot.margin = ggplot2::margin(t = 5.5, r = 10, b = 36, l = 10, unit = "pt")
    )

  ggplot2::ggsave(filename = out_png, plot = p_combo, width = 22, height = 7.8, dpi = 300L)
  invisible(TRUE)
}

focal_labels <- if (exists("focal_genera", inherits = TRUE) && length(focal_genera) > 0L) {
  as.character(focal_genera)
} else if (exists("highlight_genus_pattern", inherits = TRUE) && nzchar(highlight_genus_pattern)) {
  as.character(highlight_genus_pattern)
} else {
  character(0)
}

for (lab in unique(focal_labels)) {
  slug <- focal_slug(lab)
  focal_taxa_path <- file.path(path_tbl, paste0("focal_taxa__", slug, "__Time", target_time, "__ASV.csv"))
  out_png <- file.path(path_fig, paste0("neutral_partition_read_vs_", slug, "_side_by_side__Time", target_time, "__ASV.png"))
  out_long <- file.path(
    path_tbl,
    paste0("neutral_partition__community_vs_", slug, "_genus_by_partition__stacked_long__Time", target_time, "__ASV.csv")
  )
  make_comm_vs_focal_plot(lab, focal_taxa_path, out_png, out_long)

  if (identical(lab, "Culicoidibacter")) {
    # Legacy file naming for supplement continuity
    legacy_png <- file.path(path_fig, paste0("neutral_partition_read_vs_culicoidibacter_side_by_side__Time", target_time, "__ASV.png"))
    make_comm_vs_focal_plot(
      lab,
      file.path(path_tbl, paste0("culicoidibacter_taxa__Time", target_time, "__ASV.csv")),
      legacy_png,
      out_long,
      legacy_title = "Neutral-model partitions: All vs Culicoidibacter within each regime"
    )
  }
}

message(
  "Neutral partition stratification complete. Tables -> ", path_tbl, "; figures -> ", path_fig
)

# --- Bundle for 02__Results/08__NeutralModel.Rmd -----------------------------------------------
# --- Bundle for 02__Results/08__NeutralModel.Rmd -----------------------------------------------
bundle_rds <- file.path(path_rds, "neutral_model__bundle.rds")
fig_asv_focal <- file.path(
  path_fig,
  paste0("neutral_model_Time", target_time, "__ASV__focal_two_Culici.png")
)
fig_genus_hi <- file.path(
  path_fig,
  paste0("neutral_model_Time", target_time, "__Genus__Culicoidibacter.png")
)
fig_side <- file.path(
  path_fig,
  paste0("neutral_partition_read_vs_culicoidibacter_side_by_side__Time", target_time, "__ASV.png")
)
rds_asv <- file.path(path_rds, paste0("neutral_model_results__Time", target_time, "__ASV.rds"))
rds_genus <- file.path(path_rds, paste0("neutral_model_results__Time", target_time, "__Genus.rds"))

neutral_bundle <- list(
  ggplot_focal_four_genera = p_neutral_focal_four_genus,
  fig_neutral_focal_two_asv = fig_asv_focal,
  fig_neutral_genus_culicoidibacter = fig_genus_hi,
  fig_partition_side_by_side = fig_side,
  focal = list(
    genera = focal_genera,
    fig_neutral_genus_highlight = stats::setNames(
      lapply(focal_genera, function(g) {
        slug <- tolower(gsub("[^A-Za-z0-9]+", "_", g))
        file.path(path_fig, paste0("neutral_model_Time", target_time, "__Genus__", slug, ".png"))
      }),
      focal_genera
    ),
    fig_neutral_genus_focal_four = file.path(
      path_fig,
      paste0("neutral_model_Time", target_time, "__Genus__focal_four_genera.png")
    ),
    fig_partition_side_by_side = stats::setNames(
      lapply(focal_genera, function(g) {
        slug <- tolower(gsub("[^A-Za-z0-9]+", "_", g))
        file.path(path_fig, paste0("neutral_partition_read_vs_", slug, "_side_by_side__Time", target_time, "__ASV.png"))
      }),
      focal_genera
    ),
    fig_focal_two_asv = stats::setNames(
      lapply(focal_genera, function(g) {
        slug <- tolower(gsub("[^A-Za-z0-9]+", "_", g))
        file.path(path_fig, paste0("neutral_model_Time", target_time, "__ASV__focal_two_", slug, ".png"))
      }),
      focal_genera
    )
  ),
  rds_time60_asv = rds_asv,
  rds_time60_genus = rds_genus,
  tables_dir = path_tbl,
  meta = list(
    disclaimer = paste0(
      "Neutral partitions from a single global ASV/Genus fit (day ", target_time,
      "). Regime-stratified comparisons do not refit the model per regime."
    )
  )
)
saveRDS(neutral_bundle, bundle_rds)

sieler2026_sync_main_figures_from_manifest(driver_script = "08__NeutralModel.R")

message("08__NeutralModel complete. Outputs -> ", path_res, " | bundle -> ", bundle_rds)
