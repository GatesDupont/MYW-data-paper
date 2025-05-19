library(ggplot2)
library(cowplot)
library(dplyr)

out <- rbind(
  read.csv("output/out_df_1.csv", header = TRUE),
  read.csv("output/out_df_2.csv", header = TRUE),
  read.csv("output/out_df_3.csv", header = TRUE),
  read.csv("output/out_df_4.csv", header = TRUE),
  read.csv("output/out_df_5.csv", header = TRUE),
  read.csv("output/out_df_6.csv", header = TRUE))

out <- subset(out,species != "test")
out$design <- factor(out$design, levels = c("emp60","reg60","emp100","reg100"))
out$d0_bias <- (out$d0_est - out$d0_tru) / out$d0_tru 
out$g0_bias <- (out$g0_est - out$g0_tru) / out$g0_tru 
out$sig_bias <- (out$sig_est - out$sig_tru) / out$sig_tru 



out_all <- rbind(
  out %>%
    group_by(species,design) %>%
    summarise(parameter = "Density",
              bias = mean(d0_bias),
              lwr = quantile(d0_bias,0.05),
              upr = quantile(d0_bias,0.95)),
  out %>%
    group_by(species,design) %>%
    summarise(parameter = "Detection",
              bias = mean(g0_bias),
              lwr = quantile(g0_bias,0.05),
              upr = quantile(g0_bias,0.95)),
  out %>%
    group_by(species,design) %>%
    summarise(parameter = "Sigma",
              bias = mean(sig_bias),
              lwr = quantile(sig_bias,0.05),
              upr = quantile(sig_bias,0.95)))
gd <- ggplot(out_all, aes(x = design, y = bias, color = design, fill = design)) +
    geom_hline(yintercept = 0, color = "grey") +
    geom_errorbar(aes(ymin=lwr,ymax=upr),width=0) +
    geom_point(size=3) +
    ylab("Bias in density") +
    xlab("") +
    facet_grid(parameter~species, scales = "free") +
    scale_x_discrete(guide=guide_axis(n.dodge=2)) +
    scale_fill_brewer(palette = "Dark2") +
    scale_color_brewer(palette = "Dark2") +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

gd



out_all <- out_all %>%
  mutate(IW = upr - lwr)

giw <- ggplot(out_all, aes(x = design, y = IW, color = design, fill = design)) +
  geom_point(size=3) +
  ylab("Interval width (bias)") +
  xlab("") +
  facet_grid(parameter~species, scales = "free") +
  scale_x_discrete(guide=guide_axis(n.dodge=2)) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

giw
