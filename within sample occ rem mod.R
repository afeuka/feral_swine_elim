### Title: Within-sample validation of occupancy-abundance model
### Author: Abbey Feuka
### Date: 11JAN2024
### Notes: single or multi-chain inference

#setup------------------------
library(tidyverse)
library(pROC)
library(nimble)

#load samples -----------
# load("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/nimble/Model outputs/ziBinMod_area_04SEP24_logit_det.RData")

##load data --------------------------------
# load("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/Model Ready Data/sysbait_10day_season_nlcd_neweff.RData")

if(class(samples)=="list"){
  nChains <- length(samples)
  nmcmc <- nrow(samples[[1]])
} else {
  nChains <- 1
  nmcmc <- nrow(samples)
}

if(nbeta==3){
  subfolder<-"No NFSP"
} else {
  subfolder<-"NFSP2"
}

#occupancy ----------------------
if(nChains>1){
  psi <- 1- do.call("rbind",lapply(1:nChains,function(i){samples[[i]][,grep("pabs",colnames(samples[[i]]))]}))
} else {
  psi <- 1- samples[,grep("pabs",colnames(samples))]
}
psi_mn <- colMeans(psi)
psi_md <- sapply(1:ncol(psi),function(i)quantile(psi[,i],prob=0.5))

psi_sum <- data.frame(idx=names(psi_mn),mn=psi_mn)
psi_sum$site_idx <- sapply(1:nrow(psi_sum),function(x){
  as.numeric(gsub("\\[","",gsub("pabs","",unlist(str_split(psi_sum[x,"idx"],","))[1])))})

psi_sum$period_idx <- sapply(1:nrow(psi_sum),function(x){
  as.numeric(gsub("\\]","",unlist(str_split((psi_sum[x,"idx"]),","))[2]))})

## AUC for prevalance process ----------------
if(nChains>1){
  pocc <- do.call("rbind",lapply(1:nChains,function(i){samples[[i]][grep("pocc",colnames(samples[[i]]))]}))
} else {
  pocc <- samples[,grep("pocc",colnames(samples))]
}
pocc_mn <- colMeans(pocc)
pocc_md <- sapply(1:ncol(pocc),function(i)quantile(pocc[,i],prob=0.5))

pocc_sum <- data.frame(idx=names(pocc_mn),mn=pocc_mn,md=pocc_md)

ypred_vec <- rbinom(nrow(dat_occ),dat_occ$nweeks,pocc_sum$mn)
ypred_vec_bin <- ifelse(ypred_vec>0,1,0)
ybin <- ifelse(dat_occ$detections>0,1,0)
auc_prev_bin <- as.numeric(pROC::roc(ybin,ypred_vec_bin)$auc)
auc_prev_pocc <- as.numeric(pROC::roc(ybin,pocc_sum$mn)$auc)

## AUC for occupancy process ---------------------
psi_mn_vec <- numeric(nrow(dat_occ))
for(i in 1:nrow(dat_occ)){
  psi_mn_vec[i] <- psi_sum$mn[psi_sum$site_idx==dat_occ$site_idx[i] & 
                                psi_sum$period_idx==dat_occ$period_idx[i]]
}
auc_occ <- suppressMessages(as.numeric(pROC::roc(ybin,psi_mn_vec)$auc))
auc_prev_bin
auc_prev_pocc
auc_occ

#bayesian p-value ----------------------
pocc <- as.matrix(pocc)
occ_ll <- dev_ypred_occ <- matrix(NA,nrow(dat_occ),nmcmc*nChains)

if(nChains>1){
  ypred_occ <- as.matrix(do.call("rbind",lapply(1:nChains,function(i){samples[[i]][grep("yocc_pred",colnames(samples[[i]]))]})))
} else {
  ypred_occ <- samples[,grep("yocc_pred",colnames(samples))]
}

for(i in 1:nrow(dat_occ)){
  occ_ll[i,] <- dbinom(dat_occ$detections[i],dat_occ$nweeks[i],pocc[,i],log=T)
  dev_ypred_occ[i,] <- dbinom(ypred_occ[,i],dat_occ$nweeks[i],pocc[,i],log=T)
}

if(nChains>1){
  N_latent <- as.matrix(do.call("rbind",lapply(1:nChains,function(i){
      samples[[i]][,grepl("N_latent",colnames(samples[[i]]))]})))
  pip_a <- as.matrix(do.call("rbind",lapply(1:nChains,function(i){
    samples[[i]][,grepl("pip_a",colnames(samples[[i]]))]})))
  ypred_rem_a <- as.matrix(do.call("rbind",lapply(1:nChains,function(i){samples[[i]][,grep("yrem_pred_aerial",colnames(samples[[i]]))]})))
  pip_t <- as.matrix(do.call("rbind",lapply(1:nChains,function(i){
    samples[[i]][,grepl("pip_t",colnames(samples[[i]]))]})))
  ypred_rem_t <- as.matrix(do.call("rbind",lapply(1:nChains,function(i){samples[[i]][,grep("yrem_pred_trap",colnames(samples[[i]]))]})))
  pip_g <- as.matrix(do.call("rbind",lapply(1:nChains,function(i){
    samples[[i]][,grepl("pip_g",colnames(samples[[i]]))]})))
  ypred_rem_g <- as.matrix(do.call("rbind",lapply(1:nChains,function(i){samples[[i]][,grep("yrem_pred_ground",colnames(samples[[i]]))]})))
} else {
  N_latent <- samples[,grepl("N_latent",colnames(samples))]
  pip_a <- samples[,grepl("pip_a",colnames(samples))]
  ypred_rem_a <- samples[,grepl("yrem_pred_aerial",colnames(samples))]
  pip_t <- samples[,grepl("pip_t",colnames(samples))]
  ypred_rem_t <- samples[,grepl("yrem_pred_trap",colnames(samples))]
  pip_g <- samples[,grepl("pip_g",colnames(samples))]
  ypred_rem_g <- samples[,grepl("yrem_pred_ground",colnames(samples))]
}

samp_N_a <- rem_ll_a <-dev_ypred_rem_a <- matrix(NA,nrow(dat_aerial),nmcmc*nChains)
samp_N_t <- rem_ll_t <-dev_ypred_rem_t <- matrix(NA,nrow(dat_trap),nmcmc*nChains)
for(i in 1:nrow(dat_aerial)){
  # if(dat_aerial$elim_area_idx[i]%in%c(5,7)){
    samp_N_a[i,] <- N_latent[,paste0("N_latent[",dat_aerial$elim_area_idx[i],", ",
                                     dat_aerial$period_idx[i],"]")]
    rem_ll_a[i,] <- dbinom(x=dat_aerial$tot_rem[i],size=samp_N_a[i,],
                           prob=pip_a[,i],log=T)
    dev_ypred_rem_a[i,] <- dbinom(x=ypred_rem_a[,i],size=samp_N_a[i,],
                                  prob=pip_a[,i],log=T)
  # }
}
for(i in 1:nrow(dat_trap)){
  # if(dat_trap$elim_area_idx[i]%in%c(5,7)){
    samp_N_t[i,] <- N_latent[,paste0("N_latent[",dat_trap$elim_area_idx[i],", ",
                                     dat_trap$period_idx[i],"]")]
    rem_ll_t[i,] <- dbinom(x=dat_trap$tot_rem[i],size=samp_N_t[i,],
                           prob=pip_t[,i],log=T)
    dev_ypred_rem_t[i,] <- dbinom(x=ypred_rem_t[,i],size=samp_N_t[i,],
                                  prob=pip_t[,i],log=T)
  # }
}

rem_ll <- rbind(rem_ll_a,rem_ll_t)
dev_ypred_rem <- rbind(dev_ypred_rem_a,dev_ypred_rem_t)

dev_y <- colSums(occ_ll) + colSums(rem_ll,na.rm=T)
dev_ypred <- colSums(dev_ypred_occ) + colSums(dev_ypred_rem,na.rm = T)

pVal_ll <- sum(dev_ypred>dev_y)/(nmcmc*nChains)

pVal_occ <- sum(colSums(dev_ypred_occ)>colSums(occ_ll))/(nmcmc*nChains)
pVal_rem <- sum(colSums(dev_ypred_rem,na.rm=T)>colSums(rem_ll,na.rm=T))/(nmcmc*nChains)

pVal_ll
pVal_occ
pVal_rem

ggplot()+
  geom_histogram(aes(x=colSums(dev_ypred_occ),fill="predictions"),alpha=0.5)+
  geom_histogram(aes(x=colSums(occ_ll),fill="data"),alpha=0.5)+
  xlab('loglikelihood - occupancy')+
  ggtitle(paste("pVal =",round(pVal_occ,2)))
ggsave(filename = paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/nimble/Model validation/Within sample/Plots/",subfolder,"/pval_occ_ll.jpeg"),
       width=7,height=5,units="in",device="jpeg")

ggplot()+
  geom_histogram(aes(x=colSums(dev_ypred_rem,na.rm=T),fill="predictions"),alpha=0.5)+
  geom_histogram(aes(x=colSums(rem_ll,na.rm=T),fill="data"),alpha=0.5)+
  xlab('loglikelihood - abundance')+
  ggtitle(paste("pVal =",round(pVal_rem,2)))
ggsave(filename = paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/nimble/Model validation/Within sample/Plots/",subfolder,"/pval_rem_ll.jpeg"),
       width=7,height=5,units="in",device="jpeg")

#mean squared error -----------------------------
ypred_occ <- as.matrix(ypred_occ)
mse_occ <- mse_pred_occ <- rep(NA,nmcmc*nChains)
for(i in 1:(nmcmc*nChains)){
  mse_occ[i] <- mean((dat_occ$detections - dat_occ$nweeks*pocc[i,])^2)
  mse_pred_occ[i] <- mean((ypred_occ[i,] - dat_occ$nweeks*pocc[i,])^2)
}

ypred_rem_a <- as.matrix(ypred_rem_a)
ypred_rem_t <- as.matrix(ypred_rem_t)
N_latent <- as.matrix(N_latent)
mse_rem_a <- mse_pred_rem_a <- matrix(NA,nrow(dat_aerial),nmcmc*nChains)
mse_rem_t <- mse_pred_rem_t <- matrix(NA,nrow(dat_trap),nmcmc*nChains)
for(i in 1:nrow(dat_aerial)){
  # if(dat_aerial$elim_area_idx[i]%in%c(5,7)){
    exp_N <- samp_N_a[i,]*pip_a[,i]
    mse_rem_a[i,] <- (dat_aerial$tot_rem[i] - exp_N)^2
    mse_pred_rem_a[i,] <- (ypred_rem_a[,i] - exp_N)^2
  # }
}
for(i in 1:nrow(dat_trap)){
  # if(dat_trap$elim_area_idx[i]%in%c(5,7)){
    exp_N <- samp_N_t[i,]*pip_t[,i]
    mse_rem_t[i,] <- (dat_trap$tot_rem[i] - exp_N)^2
    mse_pred_rem_t[i,] <- (ypred_rem_t[,i] - exp_N)^2
  # }
 }

mse_rem <- colMeans(rbind(mse_rem_a,mse_rem_t),na.rm=T)
mse_pred_rem <- colMeans(rbind(mse_pred_rem_a,mse_pred_rem_t),na.rm=T)

pVal_occ_mse <- sum(mse_pred_occ>mse_occ)/(nmcmc*nChains)
pVal_rem_mse <- sum(mse_pred_rem>mse_rem)/(nmcmc*nChains)

ggplot()+
  geom_histogram(aes(x=mse_pred_occ,fill="predictions"),alpha=0.5)+
  geom_histogram(aes(x=mse_occ,fill="data"),alpha=0.5)+
  xlab('mse - occupancy')+
  ggtitle(paste("pVal =",round(pVal_occ_mse,2)))
ggsave(filename = paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/nimble/Model validation/Within sample/Plots/",subfolder,"/pval_occ_mse.jpeg"),
       width=7,height=5,units="in",device="jpeg")

ggplot()+
  geom_histogram(aes(x=mse_pred_rem,fill="predictions"),alpha=0.5)+
  geom_histogram(aes(x=mse_rem,fill="data"),alpha=0.5)+
  xlab('mse - abundance')+
  ggtitle(paste("pVal =",round(pVal_rem_mse,2)))
ggsave(filename = paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/nimble/Model validation/Within sample/Plots/",subfolder,"/pval_rem_mse.jpeg"),
       width=7,height=5,units="in",device="jpeg")


#p-value mean and var stats -------------------
## mean ----------------
# ypred_rem_a_46 <- ypred_rem_a[,dat_aerial$elim_area_idx%in%c(5,7)]
# ypred_rem_t_46 <- ypred_rem_t[,dat_trap$elim_area_idx%in%c(5,7)]

ypred_rem_mn <- rowMeans(cbind(ypred_rem_a,ypred_rem_t))
dat_rem_mn <- mean(c(dat_aerial$tot_rem,dat_trap$tot_rem))
pVal_rem_mn <- sum(ypred_rem_mn>dat_rem_mn)/(nmcmc*nChains)
ggplot()+geom_histogram(aes(x=ypred_rem_mn))+
  geom_vline(xintercept=dat_rem_mn,col="red",lwd=1.5)+
  xlab("mean - removal")+
  ggtitle(paste("pVal =",round(pVal_rem_mn,2)))
ggsave(filename = paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/nimble/Model validation/Within sample/Plots/",subfolder,"/pval_rem_mn.jpeg"),
       width=7,height=5,units="in",device="jpeg")

ypred_occ_mn <- rowMeans(ypred_occ)
dat_occ_mn <- mean(dat_occ$detections)
pVal_occ_mn <- sum(ypred_occ_mn>dat_occ_mn)/(nmcmc*nChains)
ggplot()+geom_histogram(aes(x=ypred_occ_mn))+
  geom_vline(xintercept=dat_occ_mn,col="red",lwd=1.5)+
  xlab("mean - occupancy")+
  ggtitle(paste("pVal =",round(pVal_occ_mn,2)))
ggsave(filename = paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/nimble/Model validation/Within sample/Plots/",subfolder,"/pval_occ_mn.jpeg"),
       width=7,height=5,units="in",device="jpeg")

## standard deviation ------------------
ypred_rem_var <- sqrt(apply(cbind(ypred_rem_a,ypred_rem_t),1,var))
dat_rem_var <- sqrt(var(c(dat_aerial$tot_rem,dat_trap$tot_rem)))
pVal_rem_var <- sum(ypred_rem_var>dat_rem_var)/(nmcmc*nChains)
ggplot()+geom_histogram(aes(x=ypred_rem_var))+
  geom_vline(xintercept=dat_rem_var,col="red",lwd=1.5)+
  xlab("sd - removal")+
  ggtitle(paste("pVal =",round(pVal_rem_var,2)))
ggsave(filename = paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/nimble/Model validation/Within sample/Plots/",subfolder,"/pval_rem_sd.jpeg"),
       width=7,height=5,units="in",device="jpeg")

ypred_occ_var <- sqrt(apply(ypred_occ,1,var))
dat_occ_var <- sqrt(var(dat_occ$detections))
pVal_occ_var <- sum(ypred_occ_var>dat_occ_var)/(nmcmc*nChains)
ggplot()+geom_histogram(aes(x=ypred_occ_var))+
  geom_vline(xintercept=dat_occ_var,col="red",lwd=1.5)+
  xlab("sd - occupancy")+
  ggtitle(paste("pVal =",round(pVal_occ_var,2)))
ggsave(filename = paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/nimble/Model validation/Within sample/Plots/",subfolder,"/pval_occ_sd.jpeg"),
       width=7,height=5,units="in",device="jpeg")

pVal_rem_mn
pVal_occ_mn
pVal_rem_var
pVal_occ_var

save(auc_prev_bin,auc_prev_pocc,auc_occ,
     pVal_ll,pVal_occ,pVal_rem,
     pVal_occ_mse,pVal_rem_mse,
     pVal_rem_mn,pVal_occ_mn,pVal_rem_var,pVal_occ_var,
     file = paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/nimble/Model validation/Within sample/Plots/",subfolder,"/modval_wi_occ_rem_04SEP24.Rdata"))

# load(paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/nimble/Model validation/Within sample/Plots/",subfolder,"/modval_wi_occ_rem_04SEP24.Rdata"))

# #posterior predictive checks ---------------------
# samp_rem <- dat_aerial %>% 
#   group_by(elim_area_idx)%>% 
#   summarise(total=sum(tot_rem)) %>% 
#   filter(total>100)
# rem_idx <- which(dat_aerial$elim_area_idx%in%samp_rem$elim_area_idx & dat_aerial$tot_rem>0)
# ypred_rem_a_samp <- ypred_rem_a[,rem_idx]
# dat_rem_samp <- dat_aerial[rem_idx,]
# 
# ypred_rem_sum <- data.frame(date = dat_rem_samp$Date,
#                             site_idx = dat_rem_samp$elim_area_idx,
#                             elim_area_idx = dat_rem_samp$period_idx,
#                             tot_rem = dat_rem_samp$tot_rem,
#                             mn=as.numeric(colMeans(ypred_rem_a_samp)),
#                             md=as.numeric(sapply(1:ncol(ypred_rem_a_samp),function(i)quantile(ypred_rem_a_samp[,i],probs=0.5))),
#                             lci=as.numeric(sapply(1:ncol(ypred_rem_a_samp),function(i)quantile(ypred_rem_a_samp[,i],probs=0.025))),
#                             uci=as.numeric(sapply(1:ncol(ypred_rem_a_samp),function(i)quantile(ypred_rem_a_samp[,i],probs=0.975))))
# 
# g<-list()
# for(i in 1:length(unique(ypred_rem_sum$elim_area_idx))){
#   g[[i]] <- ggplot(ypred_rem_sum %>% filter(elim_area_idx==unique(ypred_rem_sum$elim_area_idx)[i]))+
#     geom_errorbar(aes(x=date,ymin=lci,ymax=uci),col="grey")+
#     geom_point(aes(x=date,y=tot_rem,col="black"),alpha=0.5)+
#     geom_point(aes(x=date,y=mn,col="red"),alpha=0.5)+
#     scale_color_manual(values=c("black","red"),name="",labels=c("Observed","Predicted"))+
#     ggtitle(unique(dat_rem_samp$elim_area_idx)[i])
# }
# g[[1]]
# g[[2]]
# g[[3]]
# g[[4]]
# 
