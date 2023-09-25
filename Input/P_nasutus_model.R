##########################
#     Load libraries     #
##########################

library(spocc)
library(spThin)
library(dismo)
library(rgeos)
library(ENMeval)
library(dplyr)
library(rgdal)
library(rJava)
library(maptools)
library(ggplot2)
library(devtools)
library(ENMGadgets)
library(corrplot)
library(ntbox)
library(wallace)
library(usdm)
library(raster)
library(ENMTools)
library(doParallel)
library(devtools)
library(remotes)

########################
#      Load data       #
########################

setwd("~/Documentos/Projects/Promops/")

source("R/thresh.R")

##########################################
#LoadPoints

occ <- read.csv("Data/Promops nasutus R.csv", header = T)

#########################################
#Load M polygon

Mpol <- readOGR("Data/shp/Promops nasutus_mcp.shp") %>% rgeos::gBuffer(width = 2.5)

#########################################
#Read WorldClim variables
bio1 <- raster("/media/omar/HD710 PRO/Omar/QGIS/WorldClim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_1.tif")
bio2 <- raster("/media/omar/HD710 PRO/Omar/QGIS/WorldClim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_2.tif")
bio3 <- raster("/media/omar/HD710 PRO/Omar/QGIS/WorldClim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_3.tif")
bio4 <- raster("/media/omar/HD710 PRO/Omar/QGIS/WorldClim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_4.tif")
bio5 <- raster("/media/omar/HD710 PRO/Omar/QGIS/WorldClim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_5.tif")
bio6 <- raster("/media/omar/HD710 PRO/Omar/QGIS/WorldClim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_6.tif")
bio7 <- raster("/media/omar/HD710 PRO/Omar/QGIS/WorldClim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_7.tif")
bio8 <- raster("/media/omar/HD710 PRO/Omar/QGIS/WorldClim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_8.tif")
bio9 <- raster("/media/omar/HD710 PRO/Omar/QGIS/WorldClim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_9.tif")
bio10 <- raster("/media/omar/HD710 PRO/Omar/QGIS/WorldClim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_10.tif")
bio11 <- raster("/media/omar/HD710 PRO/Omar/QGIS/WorldClim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_11.tif")
bio12 <- raster("/media/omar/HD710 PRO/Omar/QGIS/WorldClim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_12.tif")
bio13 <- raster("/media/omar/HD710 PRO/Omar/QGIS/WorldClim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_13.tif")
bio14 <- raster("/media/omar/HD710 PRO/Omar/QGIS/WorldClim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_14.tif")
bio15 <- raster("/media/omar/HD710 PRO/Omar/QGIS/WorldClim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_15.tif")
bio16 <- raster("/media/omar/HD710 PRO/Omar/QGIS/WorldClim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_16.tif")
bio17 <- raster("/media/omar/HD710 PRO/Omar/QGIS/WorldClim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_17.tif")
bio18 <- raster("/media/omar/HD710 PRO/Omar/QGIS/WorldClim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_18.tif")
bio19 <- raster("/media/omar/HD710 PRO/Omar/QGIS/WorldClim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_19.tif")

envs <- stack(bio1,bio2,bio3,bio4,bio5,bio6,bio7,bio8,bio9,
              bio10,bio11,bio12,bio13,bio14,bio15,bio16,bio17,bio18,bio19)

#########################################################

# crop the environmental rasters by the background extent shape

envsBgCrop <- raster::crop(envs, Mpol)
envsBgMsk <- raster::mask(envsBgCrop, Mpol)

# sample random background points
#bg.xy <- as.data.frame(dismo::randomPoints(envsBgMsk, 50000))
bg.xy <- read.csv("Data/BackGround_P_Nasutus.csv")[,c(2,3)]
#write.csv(bg.xy, "Data/BackGround_P_Nasutus.csv")

plot(Mpol)
points(bg.xy$x, bg.xy$y)
points(occ$Lon, occ$Lat, col="red")

#########################################################

envsBgMsk <- usdm::exclude(envsBgMsk,vifstep(envsBgMsk))
names(envsBgMsk)

envs <- usdm::exclude(envs,vifstep(envsBgMsk))
names(envs)

#########################################################

t <- spThin::thin(occ, long.col = "Lon", lat.col = "Lat", spec.col = "Sp", thin.par = 10, reps = 100, locs.thinned.list.return = T, write.files = F)
t[[1]]

occs.xy <- t[[1]]
colnames(occs.xy) <- colnames(bg.xy)

#Create group partition
group.data <- ENMeval::get.block(occs = occs.xy, bg = bg.xy)

occs.grp <- group.data[[1]]
bg.grp <- group.data[[2]]

par(mfrow=c(1,2))
plot(Mpol, main="Original Occ")
points(occs.xy, col=occs.grp)
plot(Mpol, main="Background Occ")
points(bg.xy, col=bg.grp)
par(mfrow=c(1,1))

############################################################

##########################################################
#                 Build and Evaluate ENM                 #
##########################################################

# define the vector of regularization multipliers to test
rms <- seq(1, 10, 1)

# iterate model building over all chosen parameter settings
e <- ENMeval::ENMevaluate(occs.xy, envsBgMsk, bg.coords = bg.xy, 
                          RMvalues = rms, fc = c('L','Q','LQ','H','LH','LQH','LQHP'),
                          partitions = 'user', user.grp = group.data,clamp = TRUE, 
                          algorithm = "maxent.jar")

#load("Maxent/P_nasutus/e.rda")
save(e, file = "Maxent/P_nasutus/2023/e.rda")
##########################################################
#            Extract results from the analysis           #
##########################################################

evalMods <- e@models
evalTbl <- e@results
evalPreds <- e@predictions

write.csv(evalTbl, file = "Maxent/P_nasutus/2023/Table_Models.csv", quote = F, row.names = F)

##########################################################
#                  See model plots                       #
##########################################################

ENMeval::evalplot.stats(e, stats = ("auc.val.avg"), color.var = "fc", x.var = "rm", error.bars = FALSE)
ENMeval::evalplot.stats(e, stats = ("or.mtp.avg"), color.var = "fc", x.var = "rm", error.bars = FALSE)
ENMeval::evalplot.stats(e, stats = ("or.10p.avg"), color.var = "fc", x.var = "rm", error.bars = FALSE)
ENMeval::evalplot.stats(e, stats = ("AICc"), color.var = "fc", x.var = "rm", error.bars = FALSE)

##########################################################
#               Model Selection Rules                    #
##########################################################

MTPtable <- evalTbl[which(evalTbl$or.mtp.avg<0.05),]
MTPtable

AIC <- MTPtable[which(MTPtable$delta.AICc == min(MTPtable$delta.AICc, na.rm=T)), c(1,2)]
AIC

AUC <- MTPtable[which(MTPtable$auc.val.avg == max(MTPtable$auc.val.avg)), c(1,2)]
AUC

#AIC #Se guardan los outputs de MaxEnt para cada modelo (Comando solo para Linux, Windows a manita)
system(paste0("cp -R ",evalMods$rm.4_fc.LH@path," ", getwd(),"/Maxent/P_nasutus/2023/","LQHP_6")) 
#AUC
system(paste0("cp -R ",evalMods$rm.6_fc.L@path," ", getwd(),"/Maxent/P_nasutus/2023/","L_6"))

#Guardar Info Importante de Cada Modelo

png("Maxent/P_nasutus/2023/LH_4/Var_Response.png"); dismo::response(evalMods$rm.4_fc.LH); dev.off()
png("Maxent/P_nasutus/L_6/Var_Response.png");dismo::response(evalMods$rm.6_fc.L); dev.off()

#########################################################
#                     Partial ROC                       #
#########################################################

projPoly <- readOGR("Data/shp/G2_Area.shp")
predsProj <- raster::crop(envs, projPoly)
predsProj <- raster::mask(predsProj, projPoly)

predM <- dismo::predict(evalMods$rm.4_fc.LH, envsBgMsk, args="outputformat=cloglog")
occPredVals <- raster::extract(predM, occs.xy) #Cambiar M por G y viceversa
thr <- thresh(modOccVals = occPredVals, type = "p10") #mtp o p10
thr
pred.bin <- predM > thr

back <- raster::extract(pred.bin, bg.xy)
back <- bg.xy[which(back==1),]

partial_roc <- pROC(continuous_mod=predM,
                    test_data = back,
                    n_iter=1000,E_percent=5,
                    boost_percent=50,
                    parallel=FALSE)

partial_roc$pROC_summary

#PartialROC Resultas
#LH4=    Mean_AUC     Mean_pAUC_ratio_at_5%           P_value 
#       0.6001534             1.1357897             0.0000000
#L_6 =  Mean_AUC     Mean_pAUC_ratio_at_5%           P_value 
#       0.5978467          1.1332998                0.0000000 

##########################################################
#                  Predict over G                        #
##########################################################

#projPoly <- readOGR("Data/shp/G2_Area.shp")
#predsProj <- raster::crop(envs, projPoly)
#predsProj <- raster::mask(predsProj, projPoly)

predM <- dismo::predict(evalMods$rm.4_fc.LH, envsBgMsk, args="outputformat=cloglog")
plot(predM)
points(occs.xy)

predG <- dismo::predict(evalMods$rm.4_fc.LH, predsProj, args="outputformat=cloglog")
plot(predG)

#mop.out <- ntbox::mop(envsBgMsk, predsProj)
#plot(mop.out)

#########################################################
#                  Threshold and binary                 #
#########################################################

occPredVals <- raster::extract(predG, occs.xy) #Cambiar M por G y viceversa

thr <- thresh(modOccVals = occPredVals, type = "p10") #mtp o p10

thr

pred.binM <- predM > thr
plot(pred.binM)
points(occs.xy$x,occs.xy$y)

pred.binG <- predG > thr
plot(pred.binG)
points(occs.xy$x,occs.xy$y)

#########################################################

writeRaster(predM, filename = "Maxent/P_nasutus/2023/LH_4/predM_LH4", format = "GTiff")
writeRaster(predG, filename = "Maxent/P_nasutus/2023/LH_4/predG_LH4", format = "GTiff")
writeRaster(pred.binM, filename = "Maxent/P_nasutus/2023/LH_4/pred_binM_LH4", format = "GTiff")
writeRaster(pred.binG, filename = "Maxent/P_nasutus/2023/LH_4/pred_binG_LH4", format = "GTiff")
