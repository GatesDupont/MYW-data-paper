library(oSCR)
library(lubridate)
library(sf)
library(ggplot2)
library(dplyr)
library(patchwork)


# ---- Basics ----

species_list <- c("cheetah", "genet", "hyena", "leopard", "lion", "serval")
designs <- c("Leopard", "Leopard+Small")

trap_designs <- read.csv("data/MYW_grid_cam_allocation.csv")
fence <- st_read("spatial/Full_Phinda_2018_UTM.shp")

n_designs <- length(designs)
n_spp <- length(species_list)

m0_out_list <- list()

s <- 1


# ---- Loop over tests ----

# Loop through species
for(i in 1:n_spp){
    
  for(j in 1:2){
    
    cat("\n\n\nspecies: ", species_list[i], "\ndesign:  ", designs[j], "\n\n", sep = "")

    traps_to_keep <- list(c("Leopard"), c("Leopard", "Small"))[[j]]
    
    # ---- Encounter file ----
    
    # Generate species file name
    edf_file_sp <- paste0("data/", species_list[i], "_detections_grid21.csv")
    
    # Load species encounter file, remove traps from other design if necessary
    edf_sp <- read.csv(file = edf_file_sp) %>%
    left_join(trap_designs) %>%
    select(Grid_type, everything()) %>%
    filter(Grid_type %in% traps_to_keep) %>%
    select(-Grid_type)
  
    # Set first day date
    dayone <- mdy("9-17-2021")
    edf_sp$Occasion <- as.numeric(dmy(edf_sp$Date)-dayone) + 1
   
    # Set session
    edf_sp$Session <- 1
  
    
    # ---- Trap data ----- 
 
    # Load file 
    tdf <- read.csv(file = "data/MYW2021_cam_op_site_effort_UTM.csv") %>%
      left_join(trap_designs) %>%
      select(Grid_type, everything()) %>%
      filter(Grid_type %in% traps_to_keep) %>%
      select(-Grid_type)
     
    # Keeping only trap operation data
    tmp_op <- tdf[,-(1:4)]
    
    # Rename columns as Op.1...Op.n
    colnames(tmp_op) <- paste0("Op.",1:ncol(tmp_op))
    
    # Reformat operated = NA to operated = 0
    tmp_op[is.na(tmp_op)] <- 0
  
    # Convert to df and from m to km
    tdf_df <- data.frame(TrapID = tdf[,1],
                         X = tdf[,2]/1000,
                         Y = tdf[,3]/1000,
                         tmp_op)
     
    # Convert to list object
    tdf_ls <- list(tdf_df)
     
    
    # ---- Combine to oSCR data object (scrFrame) ----
    
    # Convert to scrFrame
    sp_sf <- data2oscr(edf = edf_sp,
                       sess.col = which(colnames(edf_sp) %in% "Session"),
                       id.col = which(colnames(edf_sp) %in% "Animal_ID"),
                       occ.col = which(colnames(edf_sp) %in% "Occasion"),
                       trap.col = which(colnames(edf_sp) %in% "Site"),
                       tdf = tdf_ls,
                       K = ncol(tmp_op),
                       ntraps = nrow(tmp_op),
                       remove.zeros = TRUE,
                       remove.extracaps = TRUE)$scrFrame
    
    # Collapse across ocassions
    sp_sf_collapsed <- collapse.k(sp_sf)
    sp_sf_collapsed$trapCovs[[1]][[1]]$effort_sc <- c(scale(sp_sf_collapsed$trapCovs[[1]][[1]]$effort))
  
    # ---- Make state-space data frame ----
  
    # Make state-space data frame

    buffres <- matrix(c(0, 1,    #cheetah
                        2, 0.5,  #genet
                        5, 1,    #hyena
                        5, 1,    #leopard
                        0, 1,    #lion
                        5, 1),   #serval
                      nrow=6,ncol=2,byrow = TRUE)
  
    bfence <- st_buffer(fence, dist = buffres[i,1]*1000)
    sz <- as.numeric(round(st_area(bfence) / (buffres[i,2]*1000)^2))
    sp_ss <- list(data.frame(st_coordinates(st_sample(bfence, size = sz, type="regular"))/1000))
    class(sp_ss) <- "ssDF"
    plot(fence$geometry)
    points(sp_ss[[1]][,1:2]*1000)  

    # Fit SCR model
    m <- oSCR.fit(
      model = list(~1, ~effort_sc, ~1),   
      scrFrame = sp_sf_collapsed, 
      ssDF = sp_ss, 
      encmod = "P",
      print.level = 0)
    
    m0_out_list[[s]] <- m
    s <- s+1

  }
}

save(m0_out_list, file = "output/null_models.RData")

