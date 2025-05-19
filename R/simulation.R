library(sf)
library(oSCR)
library(terra)
library(foreach)
library(tidyverse)
library(doParallel)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Simulation study
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

load("output/single_species_test_designs.RData")
load("output/null_models.RData")

fence <- st_read("spatial/Full_Phinda_2018_UTM.shp")
reg_traps_100 <- st_sample(st_buffer(fence,dist = -50), size = 100, type="regular")
reg_traps_60 <- st_sample(st_buffer(fence,dist = -50), size = 60, type="regular")

trap_labels <- read.csv("data/MYW_grid_cam_allocation.csv", header = TRUE)

emp_traps <- read.csv("data/MYW2021_cam_op_site_effort_UTM.csv", header = TRUE) %>%
  select(-Brand, -Grid_Allocation) %>%
  left_join(trap_labels, by = join_by(Site)) %>%
  select(Site, UTM_x, UTM_y, Grid_type) 

emp_traps_100 <- st_as_sf(emp_traps,coords = c(2,3))
emp_traps_60 <- st_as_sf(dplyr::filter(emp_traps, Grid_type == "Leopard"),coords = c(2,3))

par(mfrow=c(2,2), oma=c(0,0,0,0), mar=c(1,1,1,1))
plot(fence$geometry, lwd=2); plot(emp_traps_60$geometry, add = TRUE)
plot(fence$geometry, lwd=2); plot(emp_traps_100$geometry, add = TRUE)
plot(fence$geometry, lwd=2); plot(reg_traps_60, add = TRUE)
plot(fence$geometry, lwd=2); plot(reg_traps_100, add = TRUE)

trap_df <- data.frame(
  X = c(st_coordinates(emp_traps_60)[,1],st_coordinates(emp_traps_100)[,1],
        st_coordinates(reg_traps_60)[,1],st_coordinates(reg_traps_100)[,1]),
  Y = c(st_coordinates(emp_traps_60)[,2],st_coordinates(emp_traps_100)[,2],
        st_coordinates(reg_traps_60)[,2],st_coordinates(reg_traps_100)[,2]),
  Type = c(rep("Optimal",60), rep("Optimal",100), 
           rep("Regular",60), rep("Regular",100)),
  Number = factor((c(rep(60,60), rep(100,100), rep(60,60),rep(100,100))))
)

ggplot() +
  geom_sf(data = fence, fill = "white") +
  geom_point(data = trap_df, aes(x = X, y = Y)) + 
  geom_point() +
  facet_grid(Type ~ Number) +
  theme_bw() +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank())


trap_list <- list(
  emp60 = st_coordinates(emp_traps_60),
  emp100 = st_coordinates(emp_traps_100),
  reg60 = st_coordinates(reg_traps_60),
  reg100 = st_coordinates(reg_traps_60)
)

species_vector <- c("cheetah", "genet", "hyena", "leopard", "lion", "serval")

nsim <- 100

out_df <- data.frame(
  iteration = 0,
  species = "test",
  design = "test",
  d0_est = NA,
  g0_est = NA,
  sig_est = NA,
  d0_tru = NA,
  g0_tru = NA,
  sig_tru = NA,
  seed_i = NA,
  seed_j = NA
)


# Set up parallel computing
cl <- makeCluster(6)
registerDoParallel(cl)


out <- foreach(i=1:length(species_vector),
               .packages = c("sf","terra", "oSCR","tidyverse"),
               .combine='rbind') %dopar% {

  set.seed(i)
  
  dens <- exp(m0_out_list[[i]]$coef.mle$mle[4]) 
  g0 <- exp(m0_out_list[[i]]$coef.mle$mle[1]) 
  sig <- exp(m0_out_list[[i]]$coef.mle$mle[2])
  ss <- m0_out_list[[i]]$ssDF
  ssr <- as.polygons(rast(as.matrix(ss[[1]]), type="xyz"), dissolve=TRUE)
  N <- dens * nrow(ss[[1]])
  acs <- spatSample(ssr, size = N, method="random")
  acs <- crds(acs)
  
  for(j in 1:length(trap_list)){
    set.seed(j)
    trp <- data.frame(trapID = paste(1:nrow(trap_list[[j]])),
                      X = trap_list[[j]][,1]/1000,
                      Y = trap_list[[j]][,2]/1000)
                      
    dmat <- e2dist(acs,trp[,2:3])
    pmat <- g0 * exp(-dmat^2 / (2*sig^2))
    y <- matrix(0,nrow=nrow(acs),ncol=nrow(trp)) 
      
    for(s in 1:nsim){

      for(nn in 1:nrow(acs)){
        for(jj in 1:nrow(trp)){
        
          y[nn,jj] <- rpois(n = 1, lambda = pmat[nn,jj])
        
        }
      }
      
      edf <- data.frame(Session = 1,
                        ID = paste(rep(1:nrow(acs),times = nrow(trp))),
                        trapID = paste(rep(1:nrow(trp),each = nrow(acs))),
                        Occasion = 1,
                        Detected = c(y))
      edf <- edf[rep(1:nrow(edf),edf$Detected),]
      
      sf <- data2oscr(edf = edf,
                      sess.col = 1,
                      id.col = 2,
                      occ.col = 4,
                      trap.col = 3,
                      tdf = list(trp),
                      K = 1,
                      ntraps = nrow(trp),
                      remove.zeros = TRUE,
                      remove.extracaps = FALSE)$scrFrame
      
      
      mod <- oSCR.fit(scrFrame = sf, ssDF = ss, encmod = "P", se = FALSE)

      tmp_df <- data.frame(
        iteration = s,
        species = species_vector[i],
        design = names(trap_list)[j],
        d0_est = mod$coef.mle$mle[3],
        g0_est = mod$coef.mle$mle[1],
        sig_est = mod$coef.mle$mle[2],
        d0_tru = log(dens),
        g0_tru = log(g0),
        sig_tru = log(sig),
        seed_i = i,
        seed_j = j)

      out_df <- rbind(out_df,tmp_df)
      write.csv(out_df, file=paste0("output/out_df_",i,".csv"), append = FALSE, 
                row.names = FALSE, col.names = TRUE)
    }
  }
}

stopCluster(cl)

out_df <- read.csv("output/out_df.csv", header = TRUE)
