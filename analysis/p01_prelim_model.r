## p01_prelim_model.r
## Fit Poisson models with constant sex ratio

## required packages
req <- c("rstan", "ggplot2", "ggrepel", "gridExtra", "Cairo")
lapply(req, library, character.only = TRUE)

set.seed(893749)
sessionInfo()
##############################################
load("./data/d00_ci5.Rdta")

registry_key <- data.frame(registry_id=as.numeric(ci5_long$registry_cde),
                           registry_cde=ci5_long$registry_cde)
registry_key <- registry_key[!duplicated(registry_key), ]
registry_key <- registry_key[order(registry_key$registry_id), ]
rownames(registry_key) <- NULL

age_key <- data.frame(age_id=as.numeric(ci5_long$age_group),
                      age_group=ci5_long$age_group)
age_key <- age_key[!duplicated(age_key), ]
age_key <- age_key[order(age_key$age_id), ]
rownames(age_key) <- NULL

year_key <- data.frame(year_id=as.numeric(ci5_long$year_cde),
                   year=ci5_long$year_cde)
year_key <- year_key[!duplicated(year_key), ]
year_key <- year_key[order(year_key$year_id), ]
rownames(year_key) <- NULL

## try with kidney
kid <- ci5_long[ci5_long$label=="Kidney etc.",]
lung <- ci5_long[ci5_long$label=="Lung",]
thyroid <- ci5_long[ci5_long$label=="Thyroid",]
panc <- ci5_long[ci5_long$label=="Pancreas",]
crc <- ci5_long[ci5_long$label=="Colon",]
kid_datlist <- list(N=nrow(kid),
                    n_coef=1,
                    y=kid$n_ca,
                    x=as.matrix(kid$male),
                    reg=as.numeric(kid$registry_cde),
                    n_reg=max(as.numeric(kid$registry_cde)),
                    age=as.numeric(kid$age_group),
                    n_age=max(as.numeric(kid$age_group)),
                    year=as.numeric(factor(kid$year)),
                    n_year=max(as.numeric(factor(kid$year))),
                    log_offset=log(kid$n_pop))
lung_datlist <- list(N=nrow(lung),
                     n_coef=1,
                     y=lung$n_ca,
                     x=as.matrix(lung$male),
                     reg=as.numeric(lung$registry_cde),
                     n_reg=max(as.numeric(lung$registry_cde)),
                     age=as.numeric(lung$age_group),
                     n_age=max(as.numeric(lung$age_group)),
                     year=as.numeric(factor(lung$year)),
                     n_year=max(as.numeric(factor(lung$year))),
                     log_offset=log(lung$n_pop))
thyroid_datlist <- list(N=nrow(thyroid),
                        n_coef=1,
                        y=thyroid$n_ca,
                        x=as.matrix(thyroid$male),
                        reg=as.numeric(thyroid$registry_cde),
                        n_reg=max(as.numeric(thyroid$registry_cde)),
                        age=as.numeric(thyroid$age_group),
                        n_age=max(as.numeric(thyroid$age_group)),
                        year=as.numeric(factor(thyroid$year)),
                        n_year=max(as.numeric(factor(thyroid$year))),
                        log_offset=log(thyroid$n_pop))
panc_datlist <- list(N=nrow(panc),
                     n_coef=1,
                     y=panc$n_ca,
                     x=as.matrix(panc$male),
                     reg=as.numeric(panc$registry_cde),
                     n_reg=max(as.numeric(panc$registry_cde)),
                     age=as.numeric(panc$age_group),
                     n_age=max(as.numeric(panc$age_group)),
                     year=as.numeric(factor(panc$year)),
                     n_year=max(as.numeric(factor(panc$year))),
                     log_offset=log(panc$n_pop))
crc_datlist <- list(N=nrow(crc),
                    n_coef=1,
                    y=crc$n_ca,
                    x=as.matrix(crc$male),
                    reg=as.numeric(crc$registry_cde),
                    n_reg=max(as.numeric(crc$registry_cde)),
                    age=as.numeric(crc$age_group),
                    n_age=max(as.numeric(crc$age_group)),
                    year=as.numeric(factor(crc$year)),
                    n_year=max(as.numeric(factor(crc$year))),
                    log_offset=log(crc$n_pop))

kmod <- stan("./analysis/m01_poisson_ods_full_centred.stan", data=kid_datlist,
             iter=500, chains=3, cores=3, control=list(max_treedepth=15))
lmod <- stan("./analysis/m01_poisson_ods_full_centred.stan", data=lung_datlist,
             iter=500, chains=3, cores=3, control=list(max_treedepth=15))
tmod <- stan("./analysis/m01_poisson_ods_full_centred.stan", data=thyroid_datlist,
             iter=500, chains=3, cores=3, control=list(max_treedepth=15))
pmod <- stan("./analysis/m01_poisson_ods_full_centred.stan", data=panc_datlist,
             iter=500, chains=3, cores=3, control=list(max_treedepth=15))
cmod <- stan("./analysis/m01_poisson_ods_full_centred.stan", data=crc_datlist,
             iter=500, chains=3, cores=3, control=list(max_treedepth=15))

kmod
kmod_samp <- extract(kmod)
save(kmod_samp, kmod, file="./prelim_kid_samples3.Rdta")
lmod
lmod_samp <- extract(lmod)
save(lmod_samp, lmod, file="./prelim_lung_samples3.Rdta")
tmod
tmod_samp <- extract(tmod)
save(tmod_samp, tmod, file="./prelim_thyroid_samples3.Rdta")
pmod
pmod_samp <- extract(pmod)
save(pmod_samp, pmod, file="./prelim_pancreas_samples3.Rdta")
cmod
cmod_samp <- extract(cmod)
save(cmod_samp, cmod, file="./prelim_colorectal_samples3.Rdta")

s_list <- list(kidmod=kmod_samp, lmod=lmod_samp,
               tmod=tmod_samp, pmod=pmod_samp,
               cmod=cmod_samp)
c_names <- c("Kidney", "Lung", "Thyroid", "Pancreas", "Colon")
res_age <- vector(length=length(s_list), mode="list")
res_year <- vector(length=length(s_list), mode="list")
res_registry <- vector(length=length(s_list), mode="list")
res_avg <- vector(length=length(s_list), mode="list")
res_sd <- vector(length=length(s_list), mode="list")
res_sd_ratio <- vector(length=length(s_list), mode="list")
for (j in 1:length(s_list)) {
  log_rr_age <- matrix(s_list[[j]]$log_rr_age, ncol=nrow(age_key))
  log_rr_year <- matrix(s_list[[j]]$log_rr_year, ncol=nrow(year_key))
  log_rr_reg <- matrix(s_list[[j]]$log_rr_reg, ncol=nrow(registry_key))

  log_rr_mean <- s_list[[j]]$mu_b

  log_sigma_b_age <- as.vector(log(s_list[[j]]$sigma_b_age))
  log_sigma_b_year <- as.vector(log(s_list[[j]]$sigma_b_year))
  log_sigma_b_reg <- as.vector(log(s_list[[j]]$sigma_b_reg))

  log_sigma_age <- log(s_list[[j]]$sigma_log_lambda0)
  log_sigma_year <- log(s_list[[j]]$sigma_year)
  log_sigma_reg <- log(s_list[[j]]$sigma_reg)

  log_sigma_ratio_age <- log_sigma_b_age - log_sigma_age
  log_sigma_ratio_year <- log_sigma_b_year - log_sigma_year
  log_sigma_ratio_reg <- log_sigma_b_reg - log_sigma_reg

  res_summary_age <- matrix(data=0, ncol=nrow(age_key), nrow=8)
  for (i in 1L:nrow(age_key)) {
    m <- mean(log_rr_age[,i])
    pc <- quantile(log_rr_age[,i], probs=c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95))
    res_summary_age[,i] <- c(m, pc)
  }
  res_summary_age <- t(res_summary_age)
  res_summary_age_exp <- exp(res_summary_age)
  colnames(res_summary_age_exp) <- c("mean", "p05", "p10", "p25", "p50", "p75", "p90", "p95")
  res_summary_age_exp <- cbind(age_key, res_summary_age_exp)
  res_summary_age_exp$cancer <- c_names[i]

  res_summary_year <- matrix(data=0, ncol=nrow(year_key), nrow=8)
  for (i in 1L:nrow(year_key)) {
    m <- mean(log_rr_year[,i])
    pc <- quantile(log_rr_year[,i], probs=c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95))
    res_summary_year[,i] <- c(m, pc)
  }
  res_summary_year <- t(res_summary_year)
  res_summary_year_exp <- exp(res_summary_year)
  colnames(res_summary_year_exp) <- c("mean", "p05", "p10", "p25", "p50", "p75", "p90", "p95")
  res_summary_year_exp <- cbind(year_key, res_summary_year_exp)
  res_summary_year_exp$cancer <- c_names[i]

  res_summary_reg <- matrix(data=0, ncol=nrow(registry_key), nrow=8)
  for (i in 1L:nrow(registry_key)) {
    m <- mean(log_rr_reg[,i])
    pc <- quantile(log_rr_reg[,i], probs=c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95))
    res_summary_reg[,i] <- c(m, pc)
  }
  res_summary_reg <- t(res_summary_reg)
  res_summary_reg_exp <- exp(res_summary_reg)
  colnames(res_summary_reg_exp) <- c("mean", "p05", "p10", "p25", "p50", "p75", "p90", "p95")
  res_summary_reg_exp <- cbind(registry_key, res_summary_reg_exp)

  res_summary_avg <- matrix(data=0, ncol=1, nrow=8)
  m <- mean(log_rr_mean)
  pc <- quantile(log_rr_mean, probs=c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95))
  res_summary_avg[,1] <- c(m, pc)
  res_summary_avg <- t(res_summary_avg)
  res_summary_avg_exp <- exp(res_summary_avg)
  colnames(res_summary_avg_exp) <- c("mean", "p05", "p10", "p25", "p50", "p75", "p90", "p95")

  res_summary_sd <- matrix(data=0, ncol=3, nrow=8)
  m <- mean(log_sigma_b_age)
  pc <- quantile(log_sigma_b_age, probs=c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95))
  res_summary_sd[,1] <- c(m, pc)
  m <- mean(log_sigma_b_year)
  pc <- quantile(log_sigma_b_year, probs=c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95))
  res_summary_sd[,2] <- c(m, pc)
  m <- mean(log_sigma_b_reg)
  pc <- quantile(log_sigma_b_reg, probs=c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95))
  res_summary_sd[,3] <- c(m, pc)
  res_summary_sd <- t(res_summary_sd)
  res_summary_sd <- exp(res_summary_sd)
  colnames(res_summary_sd) <- c("mean", "p05", "p10", "p25", "p50", "p75", "p90", "p95")

  res_summary_sd_ratio <- matrix(data=0, ncol=3, nrow=8)
  m <- mean(log_sigma_ratio_age)
  pc <- quantile(log_sigma_ratio_age, probs=c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95))
  res_summary_sd_ratio[,1] <- c(m, pc)
  m <- mean(log_sigma_ratio_year)
  pc <- quantile(log_sigma_ratio_year, probs=c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95))
  res_summary_sd_ratio[,2] <- c(m, pc)
  m <- mean(log_sigma_ratio_reg)
  pc <- quantile(log_sigma_ratio_reg, probs=c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95))
  res_summary_sd_ratio[,3] <- c(m, pc)
  res_summary_sd_ratio <- t(res_summary_sd_ratio)
  res_summary_sd_ratio <- exp(res_summary_sd_ratio)
  colnames(res_summary_sd_ratio) <- c("mean", "p05", "p10", "p25", "p50", "p75", "p90", "p95")

  res_summary_reg_exp$cancer <- c_names[j]
  res_summary_age_exp$cancer <- c_names[j]
  res_summary_year_exp$cancer <- c_names[j]
  res_summary_avg_exp <- data.frame(res_summary_avg_exp, cancer=c_names[j])
  res_summary_sd <- data.frame(res_summary_sd, cancer=c_names[j])
  res_summary_sd$sd_by <- c("by age", "by year", "by registry")
  res_summary_sd_ratio <- data.frame(res_summary_sd_ratio, cancer=c_names[j])
  res_summary_sd_ratio$sd_by <- c("by age", "by year", "by registry")

  res_age[[j]] <- res_summary_age_exp
  res_year[[j]] <- res_summary_year_exp
  res_registry[[j]] <- res_summary_reg_exp
  res_avg[[j]] <- res_summary_avg_exp
  res_sd[[j]] <- res_summary_sd
  res_sd_ratio[[j]] <- res_summary_sd_ratio
}
res_age <- do.call(rbind.data.frame, res_age)
res_year <- do.call(rbind.data.frame, res_year)
res_registry <- do.call(rbind.data.frame, res_registry)
res_avg <- do.call(rbind.data.frame, res_avg)
res_sd <- do.call(rbind.data.frame, res_sd)
res_sd_ratio <- do.call(rbind.data.frame, res_sd_ratio)

res_sd$sd_by <- factor(res_sd$sd_by, levels=c("by registry", "by year", "by age"))
sd_relevel <- res_sd[res_sd$sd_by=="by registry",]
sd_relevel <- as.character(sd_relevel$cancer[order(sd_relevel$mean)])
res_sd$cancer <- factor(res_sd$cancer, levels=sd_relevel)

res_sd_ratio$sd_by <- factor(res_sd_ratio$sd_by, levels=c("by registry", "by year", "by age"))
sd_relevel <- res_sd[res_sd_ratio$sd_by=="by registry",]
sd_relevel <- as.character(sd_relevel$cancer[order(sd_relevel$mean)])
res_sd_ratio$cancer <- factor(res_sd$cancer, levels=sd_relevel)


registry_relevel <- res_registry[res_registry$cancer=="Lung",]
registry_relevel <- as.character(registry_relevel$registry_cde[order(registry_relevel$mean)])
res_registry$registry_cde <- factor(res_registry$registry_cde,
                                    levels=registry_relevel)

res_avg$cancer <- as.character(res_avg$cancer)
res_avg$age_group <- NA
res_avg$year <- NA
res_avg$registry_cde <- res_registry$registry_cde[1]

sdratioplot <- ggplot(data=res_sd_ratio, aes(y=mean, x=cancer, group=sd_by)) +
  geom_hline(yintercept=1, size=.5) +
  geom_linerange(aes(ymin=p10, ymax=p90), position=position_dodge(width=.2)) +
  geom_point(aes(colour=sd_by), position=position_dodge(width=.2)) +
  scale_y_log10("Ratio of SD's", breaks=c(.01, .05, .1, 5, 1, 2)) +
  scale_x_discrete("Cancer") +
  coord_flip() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line.x=element_line(size=.3),
        axis.line.y=element_line(size=.3),
        legend.title=element_blank(),
        legend.position=c(0.8,0.1),
        legend.background=element_blank(),
        legend.key=element_blank())

sdplot <- ggplot(data=res_sd, aes(y=mean, x=cancer, group=sd_by)) +
  geom_linerange(aes(ymin=p10, ymax=p90), position=position_dodge(width=.2)) +
  geom_point(aes(colour=sd_by), position=position_dodge(width=.2)) +
  scale_y_continuous("Between group SD of rate ratios ") +
  scale_x_discrete("Cancer") +
  coord_flip() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line.x=element_line(size=.3),
        axis.line.y=element_line(size=.3),
        legend.title=element_blank(),
        legend.position=c(0.8,0.1),
        legend.background=element_blank(),
        legend.key=element_blank())

aplot <- ggplot(data=res_age, aes(y=mean,x=age_group, colour=cancer, fill=cancer, group=cancer)) +
  geom_rect(data=res_avg, aes(ymax=p90, ymin=p10, xmin=-Inf, xmax=Inf), colour=NA,  alpha=.1) +
#  geom_hline(data=res_avg, aes(yintercept=mean, colour=cancer), size=.5, alpha=.4) +
  geom_point() +
  geom_line(size=.3, alpha=.5) +
  geom_linerange(aes(ymin=p10, ymax=p90)) +
  geom_text_repel(data=res_age[as.numeric(res_age$age_group)==max(as.numeric(res_age$age_group)), ],
                  aes(label=cancer),
                  size=3.3,
                  colour="black",
                  segment.color=NA,
                  nudge_x=0.3,
                  box.padding = unit(0.1, "lines"),
                  point.padding = unit(0.2, "lines")) +
  scale_y_log10("Rate ratio, 80% credible interval",
                breaks=c(.25,.5,1,2,4,8)) +
  scale_x_discrete("Age", expand=c(.15, 0)) +
#  coord_cartesian(xlim=c(min(as.numeric(res_age$age_group))-1,
#                         max(as.numeric(res_age$age_group))*1.15)) +
#  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line.x=element_line(size=.3),
        axis.line.y=element_line(size=.3),
        axis.text.x=element_text(angle=45, hjust=1),
        legend.position="none")

yplot <- ggplot(data=res_year, aes(y=mean,x=year, colour=cancer, group=cancer)) +
  geom_rect(data=res_avg, aes(ymax=p90, ymin=p10, xmin=-Inf, xmax=Inf,
                              fill=cancer), colour=NA,  alpha=.1) +
#  geom_hline(data=res_avg, aes(yintercept=mean, colour=cancer), size=.5, alpha=.4) +
  geom_point() +
  geom_line(size=.3, alpha=.5) +
  geom_linerange(aes(ymin=p10, ymax=p90)) +
  geom_text_repel(data=res_year[as.numeric(res_year$year)==max(as.numeric(res_year$year)), ],
                  aes(label=cancer),
                  size=3.2,
                  colour="black",
                  segment.color=NA,
                  nudge_x=0.3,
                  box.padding = unit(0.1, "lines"),
                  point.padding = unit(0.2, "lines")) +
  scale_y_log10("Rate ratio, 80% credible interval",
                breaks=c(.25,.5,1,2,4,8)) +
  scale_x_discrete("Year", expand=c(.15, 0)) +
#  coord_cartesian(xlim=c(min(as.numeric(res_year$year)-2),
#                         max(as.numeric(res_year$year))*1.15)) +
#  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line.x=element_line(size=.3),
        axis.line.y=element_line(size=.3),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="none")

rplot <- ggplot(data=res_registry, aes(y=mean,x=registry_cde, colour=cancer, group=cancer)) +
  geom_rect(data=res_avg, aes(x=NULL, ymax=p90, ymin=p10, xmin=-Inf, xmax=Inf,
                              fill=cancer), colour=NA,  alpha=.1) +
#  geom_hline(data=res_avg, aes(yintercept=mean, colour=cancer), size=.5, alpha=.4) +
  geom_point() +
  geom_line(size=.3, alpha=.5) +
  geom_linerange(aes(ymin=p10, ymax=p90)) +
  geom_text_repel(data=res_registry[as.numeric(res_registry$registry_cde)==max(as.numeric(res_registry$registry_cde)), ],
                  aes(label=cancer),
                  size=3.2,
                  colour="black",
                  segment.color=NA,
                  nudge_x=0.3,
                  box.padding = unit(0.1, "lines"),
                  point.padding = unit(0.2, "lines")) +
  scale_y_log10("Rate ratio, 80% credible interval",
                breaks=c(.25,.5,1,2,4,8)) +
  scale_x_discrete("Registry", expand=c(.15, 0)) +
#  coord_cartesian(xlim=c(min(as.numeric(res_registry$registry_cde))-2,
#                         max(as.numeric(res_registry$registry_cde))*1.15)) +
#  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line.x=element_line(size=.3),
        axis.line.y=element_line(size=.3),
        axis.text.x=element_text(angle = 45, hjust = 1),
        legend.position="none")

CairoPDF(file="./analysis/output/g01_prelim_age.pdf",
         width=6, height=4)
aplot
dev.off()
CairoPDF(file="./analysis/output/g01_prelim_year.pdf",
         width=8, height=4)
yplot
dev.off()
CairoPDF(file="./analysis/output/g01_prelim_registry.pdf",
         width=9, height=5)
rplot
dev.off()
CairoPDF(file="./analysis/output/g01_prelim_sdplot.pdf",
         width=4, height=5)
sdplot
dev.off()


CairoPDF(file="./analysis/output/g01_prelim_sdratioplot.pdf",
         width=6, height=5)
sdratioplot
dev.off()
