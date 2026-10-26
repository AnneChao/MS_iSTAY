library(tidyverse)
library(gridExtra)
library(lmerTest)
library(mvtnorm)
library(RColorBrewer)

output_dir <- "figures"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
save_figure <- function(p, name, width, height) {
  # ggsave(file.path(output_dir, paste0(name, ".pdf")), p,
  #        width = width, height = height, bg = "white")
  ggsave(file.path(output_dir, paste0(name, ".png")), p,
         width = width, height = height, dpi = 300, bg = "white")
}
plot_theme <- function(size = 12) {
  theme_bw(base_size = size) +
    theme(panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
          panel.grid.minor = element_line(colour = "grey94", linewidth = 0.2),
          panel.border = element_rect(colour = "grey25", fill = NA, linewidth = 0.4),
          legend.position = "bottom", legend.title = element_blank(),
          legend.key = element_blank(), strip.background = element_rect(fill = "grey85"))
}

# Appendix S1: Figure S1
S_q <- function(data, q){
  data <- as.matrix(data)
  data <- data[rowSums(data) > 0, , drop = FALSE]
  T_step <- ncol(data)
  K <- nrow(data)
  if(K >= 1){
    gamma_data <- colSums(data)
    Z_gamma <- sum(gamma_data)
    if (Z_gamma == 0) {
      return(list(S_alpha=NA, S_beta=NA, max_S_beta=NA, S_gamma=NA))
    }
    p_t_gamma <- gamma_data / Z_gamma
    p_t_gamma <- p_t_gamma[p_t_gamma > 0]
    if(q != 1){
      I_gamma <- sum(p_t_gamma^q)^(1/(1-q))
    } else {
      I_gamma <- exp(-sum(p_t_gamma * log(p_t_gamma)))
    }
    if (T_step == 1) {
      S_gamma <- 1
    } else {
      S_gamma <- (I_gamma - 1) / (T_step - 1)
    }
    if (K == 1) {
      return(list(S_alpha=S_gamma, S_beta=0, max_S_beta=0, S_gamma=S_gamma))
    }
    return(list(S_gamma=S_gamma))
  }
}
calculate_stability_profile <- function(series_data, series_name, q_range) {
  data_matrix <- matrix(series_data, nrow = 1)
  profile <- data.frame(q = q_range, Stability = numeric(length(q_range)), Series = series_name)
  for (i in 1:length(q_range)) {
    q_val <- q_range[i]
    stability_result <- S_q(data_matrix, q = q_val)
    profile$Stability[i] <- stability_result$S_gamma
  }
  return(profile)
}
series_A_data <- c(0.999, 0.001)
series_B_data <- c(rep(0.3049, 3), rep(0.0122, 6), 0.0121)
series_C_data <- c(0.25, 0.3, 0.45)
series_D_data <- c(0.8, rep(0.02, 2), rep(0.01, 16))
q_sequence <- seq(0, 2, by = 0.01)
profile_A <- calculate_stability_profile(series_A_data, "A", q_sequence)
profile_B <- calculate_stability_profile(series_B_data, "B", q_sequence)
profile_C <- calculate_stability_profile(series_C_data, "C", q_sequence)
profile_D <- calculate_stability_profile(series_D_data, "D", q_sequence)
data_for_plot_a <- rbind(profile_A, profile_B)
data_for_plot_b <- rbind(profile_C, profile_D)
profile_style <- plot_theme(10) + theme(legend.position = "none")
plot_a <- ggplot(data_for_plot_a, aes(q, Stability, colour = Series, linetype = Series)) +
  geom_line(linewidth = 0.65) +
  scale_colour_manual(values = c(A = "darkorange", B = "#99CC00")) +
  scale_linetype_manual(values = c(A = "solid", B = "dashed")) +
  scale_x_continuous(breaks = seq(0, 2, 0.5)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  coord_cartesian(ylim = c(0, 1.02)) +
  annotate("text", x = 1.18, y = 0.09, label = "Series A", colour = "darkorange", size = 5) +
  annotate("text", x = 0.87, y = 0.54, label = "Series B", colour = "#99CC00", size = 5) +
  labs(title = "(a)", x = "Order q", y = "Stability") + profile_style
plot_b <- ggplot(data_for_plot_b, aes(q, Stability, colour = Series, linetype = Series)) +
  geom_line(linewidth = 0.65) +
  scale_colour_manual(values = c(C = "red", D = "#337AB7")) +
  scale_linetype_manual(values = c(C = "solid", D = "dashed")) +
  scale_x_continuous(breaks = seq(0, 2, 0.5)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  coord_cartesian(ylim = c(0, 1.02)) +
  annotate("text", x = 0.97, y = 0.87, label = "Series C", colour = "red", size = 5) +
  annotate("text", x = 1.10, y = 0.28, label = "Series D", colour = "#337AB7", size = 5) +
  labs(title = "(b)", x = "Order q", y = "Stability") + profile_style
fig_S1 <- arrangeGrob(plot_a, plot_b, ncol = 2)
save_figure(fig_S1, "Figure_S1_Stability_profiles", 9.5, 4.2)

# Appendix S3: Figure S1
ratios <- c(1, 2, 4, 5.8, 8, 16, 32)
rho_seq <- seq(-1, 1, length.out = 500)
ratio_colors <- c("1" = "#D55E00", "2" = "#0072B2", "4" = "#009E73",
                  "5.8" = "#CC79A7", "8" = "#E69F00", "16" = "#56B4E9", "32" = "black")
calc_ldm <- function(rho, R) (1 + R^2 + 2 * rho * R) / (1 + R)^2
plot_data <- expand.grid(rho = rho_seq, ratio = ratios) |>
  mutate(LDM = calc_ldm(rho, ratio), ratio_f = factor(ratio, levels = ratios))
dot_data <- data.frame(ratio = ratios) |>
  mutate(rho_cross = (0.5 * (1 + ratio)^2 - (1 + ratio^2)) / (2 * ratio), LDM = 0.5) |>
  filter(between(rho_cross, -1, 1))
p1 <- ggplot(plot_data, aes(rho, LDM, colour = ratio_f)) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.55) +
  geom_vline(xintercept = 0, colour = "black", linewidth = 0.55) +
  geom_hline(yintercept = 0.5, colour = "grey45", linetype = "dashed", linewidth = 0.6) +
  geom_line(linewidth = 0.85) +
  geom_point(data = dot_data, aes(x = rho_cross, y = LDM), inherit.aes = FALSE,
             colour = "red", size = 3) +
  scale_colour_manual(values = ratio_colors, breaks = as.character(ratios), labels = paste("R =", ratios)) +
  scale_x_continuous(breaks = seq(-1, 1, 0.5), limits = c(-1.04, 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = function(x) sprintf("%.2f", x),
                     limits = c(-0.05, 1.05), expand = c(0, 0)) +
  labs(x = "Correlation", y = expression(phi[LM] ~ "(Loreau & de Mazancourt)")) +
  plot_theme(12) +
  theme(legend.position = "right", panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey80", linewidth = 0.45),
        axis.title = element_text(face = "bold", size = 13), aspect.ratio = 1,
        legend.key.width = grid::unit(1, "cm"))
save_figure(p1, "Figure_S3_LDM_correlation", 7.4, 5.5)

# Appendix S3: Figures S2 and S3
# Equal-weight information-based stability components; q = 1 and q = 2 are retained.
S_q_calc_equal_weight_b <- function(data, q) {
  data <- as.matrix(data)
  data <- data[rowSums(data) > 0, , drop = FALSE]
  T <- ncol(data); K <- nrow(data)
  if (K == 0 || T < 2) return(list(S_beta = 0, max_S_beta = 0, S_alpha = 0, S_gamma = 0))
  p <- data / rowSums(data)
  w <- rep(1 / K, K)
  p_gamma <- colMeans(p)
  info <- function(x) {
    x <- x[x > 0]
    if (q == 1) exp(-sum(x * log(x))) else sum(x^q)^(1 / (1 - q))
  }
  I_gamma <- info(p_gamma)
  I_k <- apply(p, 1, info)
  I_alpha <- if (q == 1) exp(sum(w * log(I_k))) else sum(w * I_k^(1 - q))^(1 / (1 - q))
  S_gamma <- (I_gamma - 1) / (T - 1)
  S_alpha <- (I_alpha - 1) / (T - 1)
  list(S_beta = max(0, S_gamma - S_alpha),
       max_S_beta = (1 - 1 / K) * (S_gamma + 1 / (T - 1)),
       S_alpha = S_alpha, S_gamma = S_gamma)
}

phi_LM_calc <- function(data) {
  data <- as.matrix(data)
  data <- data[rowSums(data) > 0, , drop = FALSE]
  denom <- sum(sqrt(apply(data, 1, var)))
  if (denom > 0) var(colSums(data)) / denom^2 else NA_real_
}

# Simulate LDM, q = 1, and q = 2; return all three numerical results.
run_synchrony_simulation <- function(
    n_time = 20, mu_x = 10, mu_y = 10, n_reps = 100,
    sd_x_vals = c(0, 0.1, 0.5, 1, 1.5, 2, 3, 4, 5, 5.5, 6),
    sd_y_vals = c(0, 0.1, 0.5, 1, 1.5, 2, 3, 4, 5, 5.5, 6),
    rho_vals = c(-1, -.7, -.5, -.3, 0, .3, .5, .7, 1), seed = 10) {
  if (!is.null(seed)) set.seed(seed)
  params <- expand.grid(sd_x = sd_x_vals, sd_y = sd_y_vals, rho = rho_vals)
  ans <- matrix(NA_real_, nrow(params), 3)
  for (i in seq_len(nrow(params))) {
    sx <- params$sd_x[i]; sy <- params$sd_y[i]; rho <- params$rho[i]
    sigma <- matrix(c(sx^2, rho * sx * sy, rho * sx * sy, sy^2), 2)
    sims <- matrix(NA_real_, n_reps, 3)
    for (j in seq_len(n_reps)) {
      sims[j, ] <- tryCatch({
        valid <- matrix(numeric(0), ncol = 2)
        attempts <- 0L
        while (nrow(valid) < n_time) {
          attempts <- attempts + 1L
          if (attempts > 1000L) stop("Unable to generate positive observations")
          z <- mvtnorm::rmvnorm(n_time, mean = c(mu_x, mu_y), sigma = sigma)
          valid <- rbind(valid, z[z[, 1] > 1e-9 & z[, 2] > 1e-9, , drop = FALSE])
        }
        community <- t(valid[seq_len(n_time), , drop = FALSE])
        synchrony <- function(q) {
          s <- S_q_calc_equal_weight_b(community, q)
          if (s$max_S_beta > 1e-10) 1 - s$S_beta / s$max_S_beta else 1
        }
        c(phi_LM_calc(community), synchrony(1), synchrony(2))
      }, error = function(e) rep(NA_real_, 3))
    }
    ans[i, ] <- colMeans(sims, na.rm = TRUE)
    if (sx == 0 && sy == 0) ans[i, 1] <- NA_real_
  }
  names(params) <- c("sd_x", "sd_y", "rho")
  cbind(params, setNames(as.data.frame(ans), c("ldm_mean", "q1_mean", "q2_mean")))
}

# Draw a heatmap using 10 bins and the same limits as Appendix S3.
make_heatmap <- function(data, metric, lower) {
  breaks <- seq(lower, 1, length.out = 11)
  labels <- vapply(seq_len(10), function(i)
    paste0("[", format(round(breaks[i], 2), trim = TRUE), ",",
           format(round(breaks[i + 1], 2), trim = TRUE), if (i == 10) "]" else ")"), character(1))
  colors <- rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlBu"))(10))
  index <- cut(pmax(lower, pmin(1, data[[metric]])), breaks = breaks,
               include.lowest = TRUE, right = TRUE, labels = FALSE)
  df <- data |>
    mutate(bin = factor(ifelse(is.na(index), "NA", labels[index]), levels = c("NA", rev(labels))),
           rho_f = factor(rho, levels = unique(rho), labels = paste("ρ =", unique(rho))))
  ggplot(df, aes(factor(sd_x, levels = unique(sd_x)), factor(sd_y, levels = unique(sd_y)), fill = bin)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    facet_wrap(~rho_f, nrow = 3) +
    scale_fill_manual(name = "Synchrony", values = c("NA" = "grey80", setNames(colors, labels)), drop = FALSE) +
    guides(fill = guide_legend(title.position = "top", title.hjust = 0.5,
                               keywidth = grid::unit(0.7, "cm"), keyheight = grid::unit(0.45, "cm"))) +
    labs(x = "SD of the first dataset", y = "SD of the second dataset") +
    coord_fixed() + theme_minimal(base_size = 12) +
    theme(axis.title = element_text(face = "bold", size = 16),
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_rect(fill = "grey90", colour = NA),
          strip.text = element_text(face = "bold", size = 13),
          panel.border = element_rect(colour = "grey80", fill = NA),
          panel.grid = element_blank(), legend.position = "right",
          legend.title = element_text(face = "bold"))
}


heatmap_results <- run_synchrony_simulation()
fig_S3_S2 <- make_heatmap(heatmap_results, "ldm_mean", 0)
fig_S3_S3 <- make_heatmap(heatmap_results, "q2_mean", 0.7)
save_figure(fig_S3_S2, "Appendix_S3_Figure_S2", 12, 10)
save_figure(fig_S3_S3, "Appendix_S3_Figure_S3", 12, 10)

# Appendix S3: Figures S2 and S3
T_sim <- 50
K_values <- seq(2, 50, 2)
n_sim <- 1000

generate_data <- function(K, T_sim) {
  list(Gamma = matrix(rgamma(K*T_sim, shape = 0.25, scale = 4), K, T_sim),
       Exponential = matrix(rexp(K*T_sim), K, T_sim),
       `Log-normal` = matrix(rlnorm(K*T_sim), K, T_sim),
       Uniform = matrix(runif(K*T_sim), K, T_sim))
}

calc_Syn3 <- function(Z) {
  K <- nrow(Z)
  P <- Z / rowSums(Z)
  D_alpha <- 1 / mean(rowSums(P^2))
  D_gamma <- 1 / sum(colMeans(P)^2)
  D_joint <- 1 / sum((P/K)^2)
  1 - (D_gamma - D_alpha) / (D_gamma - D_alpha*D_gamma/D_joint)
}

CV <- c(Gamma = 2, Exponential = 1, `Log-normal` = sqrt(exp(1)-1), Uniform = sqrt(1/3))
summary_q2 <- map_dfr(K_values, function(K) {
  sims <- replicate(n_sim, lapply(generate_data(K, T_sim), calc_Syn3), simplify = FALSE)
  tibble(Distribution = names(CV), K = K,
         Syn3_mean = map_dbl(names(CV), ~mean(vapply(sims, function(x) x[[.x]], numeric(1)))),
         Syn3_theory = 1 - ((1-1/T_sim)*CV^2)/(1+(1-1/T_sim)*CV^2))
})

dist_order <- c("Uniform", "Exponential", "Log-normal", "Gamma")
dist_labels <- c(Uniform = "Uniform (0,1)", Exponential = "Exponential (1)",
                 `Log-normal` = "Lognormal (0,1)", Gamma = "Gamma (0.25,4)")
plot_q2_Syn3 <- summary_q2 |>
  pivot_longer(c(Syn3_mean, Syn3_theory), names_to = "Type", values_to = "Value") |>
  mutate(Type = recode(Type, Syn3_mean = "Simulation", Syn3_theory = "Theory")) |>
  bind_rows(summary_q2 |> distinct(Distribution, K) |>
              mutate(Type = "phi_LM", Value = 1/K)) |>
  mutate(Distribution = factor(Distribution, levels = dist_order))

p_q2_Syn3 <- ggplot(plot_q2_Syn3, aes(K, Value, color = Type, linetype = Type)) +
  geom_line(linewidth = 1.5) +
  facet_wrap(~Distribution, nrow = 1, labeller = labeller(Distribution = dist_labels)) +
  scale_color_manual(breaks = c("Simulation", "Theory", "phi_LM"),
                     values = c(Simulation = "black", Theory = "red", phi_LM = "blue"),
                     labels = c("Simulated mean", "Limiting value", "ϕ<sub>LM</sub> = 1/<i>K</i>")) +
  scale_linetype_manual(breaks = c("Simulation", "Theory", "phi_LM"),
                        values = c(Simulation = "solid", Theory = "dashed", phi_LM = "dotdash"),
                        labels = c("Simulated mean", "Limiting value", "ϕ<sub>LM</sub> = 1/<i>K</i>")) +
  labs(x = "Number of datasets (K)", y = "Synchrony", color = NULL, linetype = NULL) +
  theme_bw(base_size = 20) +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22, face = "bold"),
        strip.text = element_text(size = 18, face = "bold"),
        legend.text = ggtext::element_markdown(size = 17, family = "Cambria Math"),
        legend.position = "bottom", legend.key.width = grid::unit(1.8, "cm"),
        panel.spacing = grid::unit(0.12, "cm"))

print(p_q2_Syn3)
save_figure(p_q2_Syn3, "Figure_S3_S4", 12, 6)

# Appendix S4: Figure S1
iSTAY_Multiple = function (data, order.q = c(1, 2), equal_weights = T, Alltime = TRUE, start_T = NULL, end_T = NULL)
{
  NA_num <- sum(is.na(data))
  if (NA_num != 0) {
    stop("There are some NA in the data.")
  }
  if (sum(which(order.q <= 0))) {
    stop("Order q need larger than 0.")
  }
  if (Alltime != TRUE) {
    if (is.null(start_T) | is.null(end_T)) {
      stop("Need to set the length of time series for calculating.")
    }
    if ((start_T >= end_T) | length(start_T) != 1 | length(end_T) !=
        1) {
      stop("Starting and ending time need to be a number, and ending time needs larger than starting time.")
    }
  }
  Stabillity_Multiple <- function(ZZ, q){
    z_iplus <- apply(ZZ, 1, sum)
    if (length(which(z_iplus == 0)) != 0) {
      ZZ <- ZZ[-which(z_iplus == 0), ]
    }
    ZZ <- as.matrix(ZZ)
    z_iplus <- apply(ZZ, 1, sum)
    Times <- ncol(ZZ); K <- nrow(ZZ)
    p_i <- as.data.frame(apply(ZZ, 2, function(w) w / z_iplus))
    if(dim(p_i)[2] == 1) p_i = t(p_i)
    if(equal_weights){
      p_pool <-  colMeans(p_i)
      ww <- rep(1/K,K)
    }else{
      z_plusk <- apply(ZZ, 2, sum)
      ww <- z_iplus/sum(ZZ)
      p_pool <- ww %*% as.matrix(p_i)
    }
    if(q == 1){
      alpha <- (exp(sum(ww * apply(p_i, 1, function(p) -sum(p[p > 0] * log(p[p > 0]))))) - 1)/(Times - 1)
      gamma <- (exp(-sum(p_pool[p_pool > 0] * log(p_pool[p_pool > 0]))) - 1)/(Times - 1)
    }else{
      alpha <- (sum(ww * rowSums(p_i^q))^(1/(1-q)) - 1)/(Times - 1)
      gamma <- (sum(p_pool^q)^(1/(1 - q)) - 1)/(Times - 1)
    }
    return( c(gamma, alpha, (gamma/alpha), (gamma - alpha)) )
  }
  Synchrony <- function(ZZ, q){
    z_iplus <- apply(ZZ, 1, sum)
    if (length(which(z_iplus == 0)) != 0) {
      ZZ <- ZZ[-which(z_iplus == 0), ]
    }
    ZZ <- as.matrix(ZZ)
    Times <- ncol(ZZ); K <- nrow(ZZ)
    z_iplus <- apply(ZZ, 1, sum)
    gamma_diversity = Stabillity_Multiple(ZZ,q)[1] * (Times - 1) + 1
    alpha_diversity = Stabillity_Multiple(ZZ,q)[2] * (Times - 1) + 1
    if(sum(rowSums(ZZ)>0)!=1){
      if(dim(ZZ)[1] > dim(ZZ)[2]){
        1 - (1 - alpha_diversity/gamma_diversity)/(1 - 1/gamma_diversity)
      }else{
        p_i <- as.data.frame(apply(ZZ, 2, function(w) w / z_iplus))
        if(equal_weights){
          ww <- rep(1/K, K)
          p_joint <- diag(ww, K) %*% as.matrix(p_i)
        }else{
          ww <- z_iplus/sum(ZZ)
          p_joint <- diag(ww, K) %*% as.matrix(p_i)
        }
        if(q == 1){
          joint_diversity <- exp(-sum(p_joint[p_joint > 0] * log(p_joint[p_joint > 0])) )
        }else{
          joint_diversity <- sum((p_joint)^q)^(1/(1 - q))
        }
        1 - (1 - alpha_diversity/gamma_diversity)/(1 - alpha_diversity/joint_diversity)
      }
    }else{
      1
    }
  }
  if (is.data.frame(data) | is.matrix(data)) {
    if (Alltime == TRUE) {
      subdata <- data
    }else {
      subdata <- data[, c(start_T:end_T)]
    }
    out <- as.matrix(sapply(order.q, function(qq) c(Stabillity_Multiple(subdata,
                                                                        q = qq), Synchrony(subdata, q = qq))))
    result <- data.frame(Site = rep(1, length(order.q)),
                         Order_q = order.q, t(out))
    colnames(result)[3:7] <- c("Stab_Gamma", "Stab_Alpha",
                               "Stab_Beta_multiplicative", "Stab_Beta_additive",
                               "Synchrony")
    result <- result[, -5]
  }else {
    out <- lapply(order.q, function(qq) {
      cal <- lapply(data, function(ZZ) {
        if (Alltime == TRUE) {
          subZZ <- ZZ
        }else {
          subZZ <- ZZ[, c(start_T:end_T)]
        }
        outout <- c(Stabillity_Multiple(subZZ, q = qq),
                    Synchrony(subZZ, q = qq))
        result <- data.frame(Order_q = qq, t(outout))
        colnames(result)[2:6] <- c("Stab_Gamma", "Stab_Alpha",
                                   "Stab_Beta_multiplicative", "Stab_Beta_additive",
                                   "Synchrony")
        result <- result[, -4]
        return(result)
      })
      cal2 <- do.call(rbind, cal)
      calcal <- data.frame(Site = names(data), cal2)
      return(calcal)
    })
    result <- do.call(rbind, out)
  }
  colnames(result) <- c("Dataset", "Order_q", "Gamma", "Alpha",
                        "Beta", "Synchrony")
  return(result)
}
load("data/Data_Jena_76_community_populations.rda")
fit_jena <- function(data, metric, block_data) {
  computed <- iSTAY_Multiple(data, order.q = c(0.5, 1, 2), equal_weights = T)
  stopifnot(nrow(computed) == 3 * nrow(block_data))
  full <- data.frame(block = rep(block_data$block, 3),
                     sowndiv = rep(block_data$sowndiv, 3),
                     value = computed[[metric]],
                     version = factor(paste0("q = ", computed$Order_q),
                                      levels = c("q = 0.5", "q = 1", "q = 2")))
  bind_rows(lapply(split(full, full$version), function(df) {
    fit <- lmerTest::lmer(value ~ sowndiv + (1 + sowndiv | block), data = df)
    df$predicted <- predict(fit, newdata = df, re.form = NA)
    df$significance <- if (summary(fit)$coefficients[2, 5] < 0.05) "significant" else "non-significant"
    df$slope <- unname(lme4::fixef(fit)[2])
    df
  }))
}
plot_names <- names(Data_Jena_76_community_populations)
name_parts <- str_split(plot_names, "_")
block_data <- data.frame(block = map_chr(name_parts, 2),
                         sowndiv = log2(as.numeric(map_chr(name_parts, 3))))
keep <- block_data$sowndiv != 0
windows <- lapply(c(2005, 2010, 2015), function(year) as.character(year:(year + 9)))
window_labels <- vapply(windows, function(x) paste0(x[1], "-", tail(x, 1)), character(1))
window_colors <- setNames(c("#EA0000", "#0066CC", "purple"), window_labels)
run_windows <- function(metric) {
  bind_rows(lapply(seq_along(windows), function(i) {
    dat <- lapply(Data_Jena_76_community_populations, function(x) {
      x <- x[, windows[[i]], drop = FALSE]
      x[x == 0] <- 1e-15
      x
    })
    if (metric %in% c("Synchrony", "Asynchrony")) {
      dat <- dat[keep]
      structure <- block_data[keep, , drop = FALSE]
    } else structure <- block_data
    fit <- fit_jena(dat, if (metric == "Asynchrony") "Synchrony" else metric, structure)
    if (metric == "Asynchrony") {
      fit$value <- 1 - fit$value
      fit$predicted <- 1 - fit$predicted
      fit$slope <- -fit$slope
    }
    fit$year <- window_labels[i]
    fit
  }))
}
plot_jena <- function(metric, tag, y_label) {
  dat <- run_windows(metric) |>
    mutate(year = factor(year, levels = window_labels),
           significance = factor(significance, levels = c("significant", "non-significant")))
  slopes <- dat |>
    group_by(year, version) |>
    summarise(slope = first(slope), .groups = "drop") |>
    mutate(label = sprintf("Slope = %.4f", slope),
           label_x = if (metric %in% c("Synchrony", "Asynchrony")) 3.45 else 3.92,
           label_y = if (metric == "Alpha") 0.90 else if (metric == "Gamma") 0.34 else 0.19)
  ggplot() +
    geom_point(data = dat, aes(sowndiv, value, colour = year),
               shape = 1, size = 1.3, alpha = 0.35, show.legend = FALSE) +
    geom_line(data = distinct(dat, version, year, sowndiv, predicted, significance),
              aes(sowndiv, predicted, colour = year, linetype = significance, group = year),
              linewidth = 0.75) +
    geom_text(data = slopes, aes(label_x, label_y - (as.integer(year) - 1) * 0.065,
                                 label = label, colour = year), size = 3, hjust = 1,
              show.legend = FALSE) +
    scale_colour_manual(values = window_colors) +
    scale_linetype_manual(values = c(significant = "solid", `non-significant` = "dashed"),
                          labels = c("Significant", "Non-significant"), drop = FALSE) +
    scale_x_continuous(breaks = if (metric %in% c("Synchrony", "Asynchrony")) 1:4 else 0:4,
                       labels = function(x) 2^x) +
    scale_y_continuous(breaks = seq(0, 1, 0.25)) +
    coord_cartesian(ylim = c(if (metric == "Gamma") 0.15 else if (metric == "Alpha") 0.05 else 0, 1.03)) +
    facet_grid(~version) +
    labs(x = "Number of species (log2 scale)", y = y_label, tag = tag) +
    plot_theme(12)
}
fig_S4a <- plot_jena("Gamma", "(a)", "Community stability")
fig_S4b <- plot_jena("Alpha", "(b)", "Population stability")
fig_S4c <- plot_jena("Synchrony", "(c)", "Population synchrony")
fig_S4d <- plot_jena("Asynchrony", "(d)", "Population asynchrony")
save_figure(fig_S4a, "Figure_S4a", 9.5, 4)
save_figure(fig_S4b, "Figure_S4b", 9.5, 4)
save_figure(fig_S4c, "Figure_S4c", 9.5, 4)
save_figure(fig_S4d, "Figure_S4d", 9.5, 4)
# fig_S4 <- arrangeGrob(fig_S4a, fig_S4b, fig_S4c, fig_S4d, ncol = 1)
# save_figure(fig_S4, "Figure_S4_Combined", 9.5, 15)

# Appendix S5: Figure S1

iSTAY_Single = function (data, order.q = c(1, 2), Alltime = TRUE, start_T = NULL, end_T = NULL)
{
  if (is.vector(data)) {
    data <- matrix(data, nrow = 1)
  }
  NA_num <- sum(is.na(data))
  if (NA_num != 0) {
    stop("There are some NA in the data.")
  }
  if (sum(which(order.q <= 0))) {
    stop("Order q need larger than 0.")
  }
  if (Alltime != TRUE) {
    if (is.null(start_T) | is.null(end_T)) {
      stop("Need to set the length of time series for calculating.")
    }
    if ((start_T >= end_T) | length(start_T) != 1 | length(end_T) !=
        1) {
      stop("Starting and ending time need to be a number, and ending time needs larger than starting time.")
    }
  }
  stability <- function(vector, q) {
    if (sum(vector != 0) == 0) {
      out <- NA
    }
    else {
      K <- length(vector)
      if (q == 1) {
        H <- sum((vector[vector > 0]/sum(vector)) *
                   log(vector[vector > 0]/sum(vector))) * (-1)
        out <- (exp(H) - 1)/(K - 1)
      }
      else {
        up <- 1 - ((sum((vector[vector > 0]/sum(vector))^q))^(1/(1 - q)))
        out <- up/(1 - K)
      }
    }
    return(out)
  }
  if (Alltime == TRUE) {
    stab <- as.matrix(apply(data, 1, function(vec) sapply(order.q,
                                                          function(qq) stability(vec, q = qq))))
  }
  else {
    subdata <- data[, c(start_T:end_T)]
    stab <- as.matrix(apply(subdata, 1, function(vec) sapply(order.q,
                                                             function(qq) stability(vec, q = qq))))
  }
  result <- data.frame(Assemblage = rep(rownames(as.data.frame(data)),
                                        length(order.q)), Order_q = rep(order.q, each = nrow(data)),
                       Stability = as.vector(t(stab)))
  colnames(result)[1] <- c("Dataset")
  return(result)
}

load("data/Data_Jena_462_populations.rda")
Jena_tidy <- Data_Jena_462_populations |>
  mutate(rowname = row.names(Data_Jena_462_populations),
         block = str_sub(rowname, 1, 2),
         number_of_species = as.numeric(sub(".*_([0-9]+)_.*", "\\1", rowname)),
         plot = str_sub(sub(".*_(B[1234]A[0-9]{2})_.*", "\\1", rowname), 3, 5),
         species = str_sub(rowname, -7, -1),
         label = paste(block, plot, number_of_species, sep = "_"))
make_s5 <- function(q) {
  temp <- Jena_tidy |> filter(number_of_species != 1)
  year_columns <- grep("^[0-9]{4}$", names(temp), value = TRUE)
  biomass <- as.matrix(temp[, year_columns])
  temp$total_biomass <- rowSums(biomass)
  single_result <- iSTAY_Single(biomass, order.q = q)
  stability_column <- intersect(c("Stab", "Stability"), names(single_result))[1]
  if (is.na(stability_column)) stop("Cannot identify iSTAY single stability column.")
  temp$Stab <- single_result[[stability_column]]
  temp <- temp |> group_by(label) |>
    mutate(weight = total_biomass / sum(total_biomass)) |> ungroup() |>
    mutate(number_of_species = factor(number_of_species, levels = c(2, 4, 8, 16)),
           species_key = factor(if_else(species == "Ono.vic", "Ono.vic", "Species"),
                                levels = c("Species", "Ono.vic")))
  fits <- temp |> group_by(number_of_species) |> group_modify(function(x, key) {
    fit <- lm(Stab ~ weight, data = x)
    tibble(weight = sort(unique(x$weight)),
           predicted = predict(fit, newdata = data.frame(weight = sort(unique(x$weight)))))
  }) |> ungroup()
  slope_labels <- temp |> group_by(number_of_species) |>
    summarise(Slope = unname(coef(lm(Stab ~ weight))[2]), .groups = "drop")
  ono <- temp |> filter(species == "Ono.vic") |> group_by(number_of_species) |>
    summarise(mean_ono = mean(Stab), .groups = "drop")
  plot <- ggplot(temp, aes(weight, Stab)) +
    geom_hline(data = ono, aes(yintercept = mean_ono), inherit.aes = FALSE,
               linetype = "dashed", linewidth = 0.7, colour = "#E88080") +
    geom_point(size = 0.45, colour = "black", alpha = 0.6, na.rm = TRUE) +
    geom_point(aes(colour = species_key, size = species_key), alpha = 0.5, na.rm = TRUE) +
    geom_line(data = fits, aes(weight, predicted), inherit.aes = FALSE,
              colour = "#20B2AA", linewidth = 0.75) +
    geom_text(data = slope_labels, aes(x = 1.01, y = 0.085,
                                       label = sprintf("Slope = %.4f", Slope)),
              inherit.aes = FALSE, colour = "#20B2AA", hjust = 1, size = 3.1) +
    geom_text(data = ono, aes(x = 1.01, y = 0.025,
                              label = sprintf("Mean Ono.Vic = %.4f", mean_ono)),
              inherit.aes = FALSE, colour = "#E88080", hjust = 1, size = 3.1) +
    scale_colour_manual(values = c(Species = "black", Ono.vic = "#F2A09E")) +
    scale_size_manual(values = c(Species = 0.7, Ono.vic = 4)) +
    scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1.05), expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0, 0.9, 0.3), expand = expansion(mult = 0)) +
    coord_cartesian(ylim = c(-0.05, 1.15)) +
    facet_wrap(~number_of_species, nrow = 1,
               labeller = labeller(number_of_species = function(x) paste(x, "species"))) +
    labs(title = paste0("q = ", q), x = "Relative biomass", y = "Stability") +
    plot_theme(11) + theme(plot.tag.position = "topleft", panel.spacing = grid::unit(0.5, "cm")) +
    guides(colour = guide_legend(override.aes = list(size = c(1.2, 4), alpha = c(1, 0.6))),
           size = "none")
  list(plot = plot, summary = left_join(slope_labels, ono, by = "number_of_species") |> mutate(q = q))
}
# Appendix S5: Figure S1, q = 1
s5_q1 <- make_s5(1)
fig_S5_q1 <- s5_q1$plot
print(fig_S5_q1)
save_figure(fig_S5_q1, "Figure_S5_q1", 10.5, 4.5)

# Appendix S5: Figure S1, q = 2
s5_q2 <- make_s5(2)
fig_S5_q2 <- s5_q2$plot
print(fig_S5_q2)
save_figure(fig_S5_q2, "Figure_S5_q2", 10.5, 4.5)
