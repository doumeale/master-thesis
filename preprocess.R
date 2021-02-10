#######load packages#########
library(tidyverse)
library(dlookr)
library(ggplot2)
library(colorblindr)
library(ReplicationSuccess)
##########read data############
effect_size<- read.csv("~/Desktop/master thesis/Confirmation and replication studies - aggregated effect sizes.csv", header=TRUE)
Decline_effect_dat<- read.csv("~/Desktop/master thesis/Decline Effect Data final verified.csv", header=TRUE)
Decline_effect_dat["lab_match"] <-
        case_when(
                Decline_effect_dat$confirmation == 1 ~ "Confirmatory test",
                Decline_effect_dat$origlab == Decline_effect_dat$replab ~ "Self-replication",
                TRUE ~ "Independent replication"
        )
###########organize data##########
######set inital dataset
replication_proj <- effect_size[,1:2]
replication_proj[,c("d_orig","d_self_rep","d_rep1", "d_rep2","d_rep3",
                    "se_orig","se_self_rep","se_rep1","se_rep2","se_rep3",
                    "p_orig","p_self_rep","p_rep1","p_rep2",
                    "p_rep3","n_orig","n_self_rep","n_rep1","n_rep2",
                    "n_rep3")] <- 0
######create a function to import data, SE was calauated using the same method in the paper
importdat <- function(dat,labeq,se,p,d,n){
     confirmation_st <- Decline_effect_dat[eval(parse(text = paste("which(abs(Decline_effect_dat$",labeq,")",sep = ""))),] 
     ord <- confirmation_st[,1] %>% match(replication_proj$study)
     dat[,c(se,p)] <- confirmation_st[ord, ] %>%
                mutate(
                se_1 = sqrt((1 / n1750e + 1 / n1750c + d^2 / (2 * (n1750 - 2))+
                1 / n2750e + 1 / n2750c + d^2 / (2 * (n2750 - 2)))/4),
                se_new = if_else(is.na(sefull), se_1, sefull),
                pval = 2 * pnorm(abs(d / se_new), lower.tail = FALSE))%>%
       select(c("se_new","pval"))
     dat[c(d,n)] <- confirmation_st[ord,c("d","N")]
     dat["rep3lab"] <- confirmation_st[ord, "replab"]
     return(dat)
     }                                                      
###### import data from self confirmatory, self replication, and 3 independent replication studies
###### columns of rep1lab, rep2lab, rep3lab represent which lab corresponds to the data of rep1, rep2, rep3.
replication_proj <- importdat(replication_proj, labeq = "confirmation)==1",
                              "se_orig","p_orig","d_orig","n_orig")                                                         
replication_proj <- importdat(replication_proj, labeq = "selfrep)==1",
                              "se_self_rep","p_self_rep","d_self_rep","n_self_rep")
replication_proj <- importdat(replication_proj, labeq = "origlab + Decline_effect_dat$replab) ==5",
                              "se_rep1","p_rep1","d_rep1","n_rep1")
replication_proj <- importdat(replication_proj, labeq = "origlab - Decline_effect_dat$replab)==2",
                              "se_rep2","p_rep2","d_rep2","n_rep2")
replication_proj <- importdat(replication_proj, labeq = "origlab + Decline_effect_dat$replab)==3|(Decline_effect_dat$origlab+Decline_effect_dat$replab==7)",
                              "se_rep3","p_rep3","d_rep3","n_rep3")
######add relative sample size#######
replication_proj["relative_sz_selfr"] <- replication_proj$n_self_rep/replication_proj$n_orig
replication_proj["relative_sz_rep1"] <- replication_proj$n_rep1/replication_proj$n_orig
replication_proj["relative_sz_rep2"] <- replication_proj$n_rep2/replication_proj$n_orig 
replication_proj["relative_sz_rep3"] <- replication_proj$n_rep3/replication_proj$n_orig

####another format
replication_proj_dat <- data.frame(c(apply(replication_proj["study"], 1,function(x) {rep(x,5)})))
names(replication_proj_dat) <- "study"
replication_proj_dat["con_lab"] <- factor(rep(c(1,-1,0,0,0),16), levels = -1:1, 
                                        labels = c("Self-replication","Independent replication","Confirmation"))
replication_proj_dat["lab"] <- c(rbind(replication_proj[,1],replication_proj[,1],
                                       paste("Lab", replication_proj[,23]),paste("Lab",replication_proj[,24]),paste("Lab",replication_proj[,25])))

replication_proj_dat["d"] <- c(rbind(replication_proj[,3],replication_proj[,4],
                                       replication_proj[,5],replication_proj[,6],replication_proj[,7]))
replication_proj_dat["se"] <- c(rbind(replication_proj[,8],replication_proj[,9],
                                     replication_proj[,10],replication_proj[,11],replication_proj[,12]))
replication_proj_dat["p"] <- c(rbind(replication_proj[,13],replication_proj[,14],
                                     replication_proj[,15],replication_proj[,16],replication_proj[,17]))
replication_proj_dat["n"] <- c(rbind(replication_proj[,18],replication_proj[,19],
                                     replication_proj[,20],replication_proj[,21],replication_proj[,22]))

save(replication_proj,file = "replication_study.rda")
save(replication_proj_dat,file = "replication_study_dat.rda")
#######EDA########
summary(replication_proj)
with(replication_proj,{plot(p_orig~ p_self_rep,xlab = "original p", ylab = "replication p", main = "Pvalue")
                       points(p_orig~ p_rep1,pch =2,col = 2)
                       points(p_orig~ p_rep2,pch =3, col = 3)
                       points(p_orig~ p_rep3,pch =4, col = 4)
                       legend("topright",c("self replication","replication 1","replication 2",
                                           "replication 3"),col = c(1:3), pch = c(1:4))})

with(replication_proj_dat,
        interaction.plot(study, con_lab, d,
                         trace.label = "replications", xlab = "study", ylab = "mean of effect"))
with(replication_proj_dat,
     interaction.plot(study, con_lab, n,
                      trace.label = "replications", xlab = "study", ylab = "mean of study size"))

par(mfrow = c(2, 2), las = 1, mai = rep(0.65, 4))
pkind <- c("p_self_rep","p_rep1","p_rep2","p_rep3")
for (p in unique(pkind)) {
        data_project <- select(replication_proj,c("lab","study", "p_orig", p) )
        significant <- ifelse(eval(parse(text = paste("data_project$",p,"< 0.05"))), "#DF536BCC", "#00000080")
        shape <- ifelse(data_project$p_orig < 0.05, 1, 3)
        plot(eval(parse(text =paste(p," ~ p_orig"))), data = data_project, ylim = c(-0.5, 1), col = significant,
             xlim = c(-0.1, 1), main = p, xlab = expression(italic(p)[o]),
             cex = 0.7, ylab = expression(italic(p)[r]))
        legend("topleft", legend = "significant", pch = 20, col = "#DF536BCC", bty = "n")
        abline(h = 0, lty = 2)
        abline(a = 0, b = 1, col = "grey")
}

replication_proj %>% 
        select(d_orig:d_rep3) %>%
        plot_correlate(d_orig:d_rep3)

replication_proj %>% 
        select(n_orig:n_rep3) %>%
        plot_correlate(n_orig:n_rep3)
######reproduce plot#######
study_order <- 
        Decline_effect_dat %>%
        select(origlab, study) %>%
        arrange(origlab) %>%
        distinct() %>%
        mutate(
                study_lab = paste0(study, " (Lab ", origlab,")")
        )
lab_shapes <- c("circle filled","square filled","diamond filled", "triangle filled")

Decline_effect_dat["lab_study"] <- factor(Decline_effect_dat$study, 
                                          levels = study_order$study, 
                                          labels = study_order$study_lab)
Decline_effect_dat["con_lab"] <- factor(Decline_effect_dat$confirmation, levels = 0:1, 
                                         labels = c("Replications","Confirmation"))
Decline_effect_dat["rep_lab"] <- paste("Lab", Decline_effect_dat$replab)
Decline_effect_dat$d95.CIlower <- as.numeric(Decline_effect_dat$d95.CIlower)
Decline_effect_dat$d95.CIupper <- as.numeric(Decline_effect_dat$d95.CIupper)
Decline_effect_dat <- Decline_effect_dat %>% mutate(
        se_basic = sqrt((1 / n1750e + 1 / n1750c + d^2 / (2 * (n1750 - 2))+
                      1 / n2750e + 1 / n2750c + d^2 / (2 * (n2750 - 2)))/4),
        se = if_else(is.na(sefull), sqrt(se_basic), sefull)
        ) %>%mutate(
          cil = if_else(is.na(sefull),d-se*qnorm(0.975),d95.CIlower),
          ciu = if_else(is.na(sefull),d+se*qnorm(0.975),d95.CIupper)
        )

ggplot(Decline_effect_dat) + 
        geom_hline(yintercept = 0) + 
        geom_pointrange(
                aes(x = con_lab, y = d, ymin = cil, ymax = ciu, 
                    color = lab_match, fill = lab_match, 
                    shape = rep_lab, group = 4 - wave),
                position = position_dodge(width = 0.8)
        ) + 
        theme_minimal() + 
        scale_shape_manual(values = lab_shapes) + 
        scale_color_OkabeIto() + 
        scale_fill_OkabeIto() + 
        scale_x_discrete(labels = NULL, breaks = NULL) + 
        scale_y_continuous(breaks = seq(0,0.6,0.2), minor_breaks = NULL) + 
        coord_flip() + 
        facet_wrap(~ lab_study) + 
        labs(shape = "", color = "", fill = "",
             x = "", y = "Effect size (Standardized mean difference)") + 
        guides(
                color = guide_legend(override.aes = list(size = 0.6)),
                shape = "none"
        ) + 
        theme(
                legend.position = "top",
                strip.text = element_text(hjust = 0),
                axis.text.x = element_text(margin = margin(5, b = 10))
        )
#######sceptical p############
#####nomial
names(replication_proj)[26] <- "relative_sz_self_rep"
scept_p <- replication_proj[,1:2]
kind <- c("self_rep","rep1","rep2","rep3")
for (p in unique(kind)) {
scept_p[paste("orig.vs",p,sep = ".")] <- pSceptical(zo = p2z(replication_proj$p_orig),
                                         zr = p2z(eval(parse(text=paste("replication_proj$p_",p, sep = "")))), 
                                         c = eval(parse(text=paste("replication_proj$relative_sz_",p,sep = ""))),
                                         alternative = "one.sided",type = "nominal")
                          }

for (p in unique(kind[-1])) {
        scept_p[paste("selfrep.vs",p,sep = ".")] <- pSceptical(zo = p2z(replication_proj$p_self_rep),
                                                            zr = p2z(eval(parse(text=paste("replication_proj$p_",p, sep = "")))), 
                                                            c = eval(parse(text=paste("replication_proj$n_",p,"/replication_proj$n_self_rep",sep = ""))),
                                                            alternative = "one.sided",type = "nominal")
}
######golden
gold_p <- replication_proj[,1:2]
for (p in unique(kind)) {
        gold_p[paste("orig.vs",p,sep = ".")] <- pSceptical(zo = p2z(replication_proj$p_orig),
                                                            zr = p2z(eval(parse(text=paste("replication_proj$p_",p, sep = "")))), 
                                                            c = eval(parse(text=paste("replication_proj$relative_sz_",p,sep = ""))),
                                                            alternative = "one.sided")
        gold_p[paste("selfrep.vs",p,sep = ".")] <- pSceptical(zo = p2z(replication_proj$p_self_rep),
                                                               zr = p2z(eval(parse(text=paste("replication_proj$p_",p, sep = "")))), 
                                                               c = eval(parse(text=paste("replication_proj$n_",p,"/replication_proj$n_self_rep",sep = ""))),
                                                               alternative = "one.sided")

}
gold_p<- gold_p[,-4]
scept_p_dat <- data.frame(x = rep(scept_p$lab,4), y = rep(scept_p$study,4))
names(scept_p_dat) <- c("lab",'study')
scept_p_dat["(nominal)orig.vs"] <- c(scept_p[,3],scept_p[,4],scept_p[,5],scept_p[,6])
scept_p_dat["contrast1"] <- c(rep("selfrep",16),rep("rep1",16),rep("rep2",16),rep("rep3",16))
scept_p_dat["(gold)orig.vs"] <- c(gold_p[,3],gold_p[,4],scept_p[,6],scept_p[,8])

scept_p_dat2 <- data.frame(x = rep(scept_p$lab,3), y = rep(scept_p$study,3))
names(scept_p_dat2) <- c("lab",'study')
scept_p_dat2["(nominal)selfrep.vs"] <- c(scept_p[,7],scept_p[,8],scept_p[,9])
scept_p_dat2["contrast1"] <- c(rep("rep1",16),rep("rep2",16),rep("rep3",16))
scept_p_dat2["(gold)selfrep.vs"] <- c(gold_p[,5],gold_p[,7],gold_p[,9])

####plot sceptical pv
alpha <- 0.05
#levelSceptical(level = 0.025, type = "nominal")

boxplot(`(nominal)orig.vs`~ contrast1, data = scept_p_dat, las = 1, cex.axis = 0.7, ylim = c(0, 1),
        xlab = "Orig study vs", ylab = expression(italic(tilde(p))[S]), outline = FALSE,
        col = "#0000000D")
abline(h = alpha/2, lty = 1, col = 4)
axis(side = 4, at = alpha/2, col.axis = 4, col = 4, las = 1, cex.axis = 0.5)
stripchart(`(nominal)orig.vs`~ contrast1, data = scept_p_dat, vertical = TRUE, add = TRUE,
           pch = 20, method = "jitter", jitter = 0.3, cex = 1, col = "#00000099")

boxplot(`(nominal)selfrep.vs`~ contrast1, data = scept_p_dat2, las = 1, cex.axis = 0.7, ylim = c(0, 1),
        xlab = "Self replication study vs", ylab = expression(italic(tilde(p))[S]), outline = FALSE,
        col = "#0000000D")
abline(h = alpha/2, lty = 1, col = 4)
axis(side = 4, at = alpha/2, col.axis = 4, col = 4, las = 1, cex.axis = 0.5)
stripchart(`(nominal)selfrep.vs`~ contrast1, data = scept_p_dat2, vertical = TRUE, add = TRUE,
           pch = 20, method = "jitter", jitter = 0.3, cex = 1, col = "#00000099")
#########meata analysis#########    
library(metafor)
m <- metagen(d,
             se,
             data=replication_proj_dat[which(replication_proj_dat$study=="Tumor"),],
             studlab=study,
             comb.fixed = TRUE,
             comb.random = TRUE,
             prediction=TRUE,
             sm="SMD")
forest(m)

m <- metagen(d,
             se,
             data=replication_proj_dat[which(replication_proj_dat$study=="FSD"),],
             studlab=study,
             comb.fixed = TRUE,
             comb.random = TRUE,
             prediction=TRUE,
             sm="SMD")
forest(m)




                                                       