#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Post processing (null models, all species)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(oSCR)
library(tidyverse)
load("output/null_models.RData")

species_list <- c("cheetah", "genet", "hyena", "leopard", "lion", "serval")
designs <- c("Leopard", "Leopard+small")

ndf <- data.frame(Session = 1, effort_sc = 0)

pred_df <- data.frame(
  Species = rep(species_list,each = 3*2),
  Parameter = rep(c("Density", "Detection", "Sigma"),times = 12),
  Design = rep(rep(designs, each=3), times = 6),
  Estimate = c(rbind(
    unlist(lapply(m0_out_list,function(x) get.real(x, type = "dens", newdata = ndf)$estimate)),
    unlist(lapply(m0_out_list,function(x) get.real(x, type = "det", newdata = ndf)$estimate)),
    unlist(lapply(m0_out_list,function(x) get.real(x, type = "sig", newdata = ndf)$estimate)))),
  Lower_95 = c(rbind(
    unlist(lapply(m0_out_list,function(x) get.real(x, type = "dens", newdata = ndf)$lwr)),
    unlist(lapply(m0_out_list,function(x) get.real(x, type = "det", newdata = ndf)$lwr)),
    unlist(lapply(m0_out_list,function(x) get.real(x, type = "sig", newdata = ndf)$lwr)))),
  Upper_95 = c(rbind(
    unlist(lapply(m0_out_list,function(x) get.real(x, type = "dens", newdata = ndf)$upr)),
    unlist(lapply(m0_out_list,function(x) get.real(x, type = "det", newdata = ndf)$upr)),
    unlist(lapply(m0_out_list,function(x) get.real(x, type = "sig", newdata = ndf)$upr))))
)


ggplot(data = pred_df, aes(x = Design, y = Estimate)) +
  geom_errorbar(aes(ymin = Lower_95, ymax = Upper_95), width = 0) + 
  geom_point() +
  facet_grid(Parameter ~ Species, scales = "free_y") + 
  theme_bw()
  
