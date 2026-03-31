suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(readr)
  library(forcats)
  library(purrr)
  library(ggpubr)
  library(grid)
})

## ---------------------------------------------------------
## Paths
## ---------------------------------------------------------
dir_null <- file.path(base_dir, "results_null_final_updated")
dir_alt  <- file.path(base_dir, "results_alternative_final_updated")
fig_dir  <- file.path(base_dir, "figs")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

## ---------------------------------------------------------
## Settings
## ---------------------------------------------------------
my_palette <- c("#264653", "#2A9D8F", "#E9C46A", "#F4A261", "#E76F51")

keep_methods <- c("proposed_F", "wrong_ps_F", "wrong_cond_leak_F", "wrong_both_F")

## show only (i)-(iv) on x-axis
method_short_map <- c(
  "proposed_F"        = "(i)",
  "wrong_ps_F"        = "(ii)",
  "wrong_cond_leak_F" = "(iii)",
  "wrong_both_F"      = "(iv)"
)
method_levels_short <- c("(i)", "(ii)", "(iii)", "(iv)")

## legend shows full names (one per line)
method_legend_labels <- c(
  "(i)  Sample splitting and correctly specified propensity score (proposed)",
  "(ii) Sample splitting and incorrectly specified propensity score",
  "(iii) Without sample splitting and correctly specified propensity score",
  "(iv) Without sample splitting and incorrectly specified propensity score"
)
names(method_legend_labels) <- method_levels_short

## colors mapped to (i)-(iv)
color_map_short <- setNames(my_palette[1:4], method_levels_short)

## delta grid used in alternatives
deltas_power <- c(0.02, 0.04, 0.06, 0.08, 0.10)
p_levels     <- c(21, 51)
n_levels     <- c(100, 200)

## y-axis limits requested
ylims <- list(
  # FWER  = c(0, 0.75), # for n = 100
  FWER  = c(0, 1.00), # for n  = 200
  Power = c(0, 1.00),
  # FDP   = c(0, 0.70), # for n = 100
  FDP   = c(0, 1.00), # for n = 200
  FDPex = c(0, 1.00)
)

## dodge width to separate overlapping methods
pd_width <- 0.35

## ---------------------------------------------------------
## Helpers for reading raw data
## ---------------------------------------------------------
read_null_bias <- function(dir_null) {
  pat <- "results_globalnull_F_final_p\\d+_n\\d+_lam[0-9.]+_seeds\\d+_\\d+_\\d+reps_bias\\.rds$"
  files <- list.files(dir_null, pattern = pat, full.names = TRUE)
  if (!length(files)) return(NULL)
  
  files %>%
    lapply(readRDS) %>%
    bind_rows() %>%
    mutate(
      method    = as.character(method),
      fwer      = as.integer(fwer),
      mean_bias = as.numeric(mean_bias),
      mean_abs  = as.numeric(mean_abs),
      l2_mean   = as.numeric(l2_mean),
      max_abs   = as.numeric(max_abs),
      p         = as.integer(p),
      n         = as.integer(n)
    ) %>%
    filter(method %in% keep_methods)
}

read_alt <- function(dir_alt) {
  pat <- "^powerfdp_F_final_p\\d+_n\\d+_delta[0-9.]+_.*reps\\.rds$"
  files <- list.files(dir_alt, pattern = pat, full.names = TRUE)
  if (!length(files)) return(NULL)
  
  files %>%
    lapply(readRDS) %>%
    bind_rows() %>%
    mutate(
      method = as.character(method),
      p      = as.integer(p),
      n      = as.integer(n),
      delta  = round(as.numeric(delta), 3),
      power  = as.numeric(power),
      fdp    = as.numeric(fdp),
      fdx    = as.integer(fdp > 0.05)
    ) %>%
    filter(method %in% keep_methods)
}

null_df <- read_null_bias(dir_null)
alt_df  <- read_alt(dir_alt)

## ---------------------------------------------------------
## Data builders for Figure 3 panels (parameterized by n)
## ---------------------------------------------------------
build_panelA <- function(null_df, n_use) {
  if (is.null(null_df)) return(NULL)
  
  null_df %>%
    filter(n == n_use) %>%
    mutate(
      method_label = factor(recode(method, !!!method_short_map), levels = method_levels_short),
      p = factor(p, levels = p_levels)
    ) %>%
    group_by(method_label, p) %>%
    summarise(Value = mean(fwer), .groups = "drop")
}

build_panelBCD <- function(alt_df, n_use) {
  if (is.null(alt_df)) return(list(B=NULL, C=NULL, D=NULL))
  
  base <- alt_df %>%
    filter(n == n_use, delta %in% deltas_power) %>%
    mutate(
      method_label = factor(recode(method, !!!method_short_map), levels = method_levels_short),
      p = factor(p, levels = p_levels),
      delta = factor(sprintf("%.2f", delta), levels = sprintf("%.2f", deltas_power))
    )
  
  panelB <- base %>%
    group_by(method_label, p, delta) %>%
    summarise(Value = mean(power, na.rm = TRUE), .groups = "drop")
  
  panelC <- base %>%
    group_by(method_label, p, delta) %>%
    summarise(Value = mean(fdp, na.rm = TRUE), .groups = "drop")
  
  panelD <- base %>%
    group_by(method_label, p, delta) %>%
    summarise(Value = mean(fdx, na.rm = TRUE), .groups = "drop")
  
  list(B = panelB, C = panelC, D = panelD)
}

## ---------------------------------------------------------
## Legend helper: bottom, one item per line
## ---------------------------------------------------------
legend_bottom_1perline <- function() {
  guides(color = guide_legend(
    nrow = 4, byrow = TRUE,
    override.aes = list(linewidth = 0.8, shape = 16)
  ))
}

## ---------------------------------------------------------
## Plot helpers (with position_dodge to avoid overlap)
## ---------------------------------------------------------
plot_FWER_single_n <- function(df, y_lim) {
  pd <- position_dodge(width = pd_width)
  
  ggplot(df, aes(x = method_label, y = Value,
                 color = method_label,
                 group = method_label)) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey65") +
    geom_point(position = pd, size = 2.6) +
    geom_line(position = pd, linewidth = 0.7) +
    facet_wrap(
      ~ p, nrow = 1,
      labeller = labeller(p = function(x) paste0("p = ", x))
    ) +
    scale_color_manual(
      values = color_map_short,
      breaks = method_levels_short,
      labels = method_legend_labels,
      name   = "Method"
    ) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = y_lim) +
    labs(title = expression("FWER (" * delta == 0 * ")"),
         x = "Method", y = "FWER") +
    theme_minimal(base_size = 13) +
    guides(color = guide_legend(nrow = 4, byrow = TRUE)) +
    theme(
      panel.grid.minor = element_blank(),
      
      ## spacing between p=21 and p=51
      panel.spacing.x = unit(2.2, "lines"),
      
      ## facet labels
      strip.text = element_text(size = 15, face = "bold"),
      
      ## legend
      legend.position = "bottom",
      legend.title = element_text(size = 15, face = "bold"),
      legend.text  = element_text(size = 13),
      
      axis.text.x = element_text(angle = 0, hjust = 0.5)
    )
}


plot_vs_delta_single_n <- function(df, title_expr, y_lim, y_lab = "Rate") {
  pd <- position_dodge(width = pd_width)
  
  ggplot(df, aes(x = delta, y = Value,
                 color = method_label,
                 group = method_label)) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey85") +
    geom_point(position = pd, size = 2.6) +
    geom_line(position = pd, linewidth = 0.7) +
    facet_wrap(
      ~ p, nrow = 1,
      labeller = labeller(p = function(x) paste0("p = ", x))
    ) +
    scale_color_manual(
      values = color_map_short,
      breaks = method_levels_short,
      labels = method_legend_labels,
      name   = "Method"
    ) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = y_lim) +
    labs(title = title_expr, x = expression(delta), y = y_lab) +
    theme_minimal(base_size = 13) +
    guides(color = guide_legend(nrow = 4, byrow = TRUE)) +
    theme(
      panel.grid.minor = element_blank(),
      
      ## spacing between columns
      panel.spacing.x = unit(2.2, "lines"),
      
      ## facet labels
      strip.text = element_text(size = 15, face = "bold"),
      
      ## legend
      legend.position = "bottom",
      legend.title = element_text(size = 15, face = "bold"),
      legend.text  = element_text(size = 13),
      
      axis.text.x = element_text(angle = 20, hjust = 1)
    )
}


## ---------------------------------------------------------
## Build & save Figure 3 (main): n = 100 only
## ---------------------------------------------------------
n_main <- 100
panelA_main <- build_panelA(null_df, n_main)
panels_main <- build_panelBCD(alt_df, n_main)

if (!is.null(panelA_main) && !is.null(panels_main$B) && !is.null(panels_main$C) && !is.null(panels_main$D)) {
  
  gA_main <- plot_FWER_single_n(panelA_main, ylims$FWER)
  gB_main <- plot_vs_delta_single_n(panels_main$B, expression("Power (" * delta * " > 0)"),
                                    ylims$Power, y_lab = "Power")
  gC_main <- plot_vs_delta_single_n(panels_main$C, expression("FDP (" * delta * " > 0)"),
                                    ylims$FDP, y_lab = "FDP")
  gD_main <- plot_vs_delta_single_n(panels_main$D, expression("FDPex (" * delta * " > 0)"),
                                    ylims$FDPex, y_lab = "FDPex")
  
  ## Make each subfigure less tall:
  ## - reduce relative heights
  ## - reduce overall saved height
  fig3_main <- ggarrange(
    gA_main, gB_main, gC_main, gD_main,
    labels = c("(a)", "(b)", "(c)", "(d)"),
    ncol = 1,
    heights = c(2.8, 2.8, 2.8, 3.1),  # <-- smaller than before
    align = "v",
    common.legend = TRUE,
    legend = "bottom"
  )
  
  fig3_height <- 14.2  # <-- reduce total height (was ~18.5 before)
  ggsave(file.path(fig_dir, "Figure3_main_n100.pdf"),
         fig3_main, width = 11.5, height = fig3_height, device = cairo_pdf)
}

## ---------------------------------------------------------
## Build & save Figure 3 (Appendix): n = 200 only
## ---------------------------------------------------------
n_app <- 200
panelA_app <- build_panelA(null_df, n_app)
panels_app <- build_panelBCD(alt_df, n_app)

if (!is.null(panelA_app) && !is.null(panels_app$B) && !is.null(panels_app$C) && !is.null(panels_app$D)) {
  
  gA_app <- plot_FWER_single_n(panelA_app, ylims$FWER)
  gB_app <- plot_vs_delta_single_n(panels_app$B, expression("Power (" * delta * " > 0)"),
                                   ylims$Power, y_lab = "Power")
  gC_app <- plot_vs_delta_single_n(panels_app$C, expression("FDP (" * delta * " > 0)"),
                                   ylims$FDP, y_lab = "FDP")
  gD_app <- plot_vs_delta_single_n(panels_app$D, expression("FDPex (" * delta * " > 0)"),
                                   ylims$FDPex, y_lab = "FDPex")
  
  fig3_app <- ggarrange(
    gA_app, gB_app, gC_app, gD_app,
    labels = c("(a)", "(b)", "(c)", "(d)"),
    ncol = 1,
    heights = c(2.8, 2.8, 2.8, 2.8),
    align = "v",
    common.legend = TRUE,
    legend = "bottom"
  )
  
  fig3_height <- 14.2
  ggsave(file.path(fig_dir, "Appendix_Figure3_n200.pdf"),
         fig3_app, width = 11.5, height = fig3_height, device = cairo_pdf)
}

## ---------------------------------------------------------
## Bias / RMSE Tables
## ---------------------------------------------------------
make_bias_rmse_tables <- function(null_df,
                                  p_levels = c(21, 51),
                                  n_levels = c(100, 200),
                                  digits = 3) {
  
  ## 1) Summarize once
  summ <- null_df %>%
    filter(p %in% p_levels, n %in% n_levels) %>%
    mutate(
      Method = factor(recode(method, !!!method_short_map), levels = method_levels_short),
      p = factor(p, levels = p_levels),
      n = factor(n, levels = n_levels),
      RMSE = sqrt(l2_mean),
      pn = paste0("p=", p, ", n=", n)
    ) %>%
    group_by(Method, pn) %>%
    summarise(
      Bias = mean(mean_bias, na.rm = TRUE),
      RMSE = mean(RMSE, na.rm = TRUE),
      .groups = "drop"
    )
  
  pn_order <- as.vector(outer(p_levels, n_levels, FUN = function(p, n) paste0("p=", p, ", n=", n)))
  
  ## Helper to build a gt table for one metric
  build_one <- function(metric = c("Bias", "RMSE")) {
    metric <- match.arg(metric)
    
    tab <- summ %>%
      select(Method, pn, all_of(metric)) %>%
      pivot_wider(names_from = pn, values_from = all_of(metric)) %>%
      arrange(Method) %>%
      select(Method, all_of(pn_order))  # enforce column order
    
    # Column labels (short and clean)
    col_labels <- c(Method = "Method")
    for (pn in pn_order) col_labels[[pn]] <- pn
    
    gt(tab) %>%
      cols_label(.list = col_labels) %>%
      fmt_number(columns = -Method, decimals = digits) %>%
      cols_width(
        Method ~ px(60),
        everything() ~ px(120)
      ) %>%
      tab_options(
        table.font.size = 11,
        data_row.padding = px(4),
        heading.align = "left"
      )
  }
  
  list(
    bias_table = build_one("Bias"),
    rmse_table = build_one("RMSE")
  )
}

tabs <- make_bias_rmse_tables(null_df)

gtsave(tabs$bias_table, file = file.path(fig_dir, "bias_table.tex"))
gtsave(tabs$rmse_table, file = file.path(fig_dir, "rmse_table.tex"))
