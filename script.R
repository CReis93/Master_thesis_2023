#########################################
#         
# Analyses statistiques des donnees de terrain / laboratoire
# Travail de master 2022-2023
# Auteur : Christophe Reis
# Version du 06 juin 2023
#
#########################################

# sommaire

# L. 26-58 : preparation du script
# L. 60-254 : variances des indices de condition
# L. 256-310 : variances des metallothioneines
# L. 312-1288 : variances metaux (chair et coquille)
# L. 1290-1361 : fonctions diverses (pour les analyses de correlations)
# L. 1363-1395 : analyse de correlation IC-IC
# L. 1397-1448 : analyse de correlation MTs-IC
# L. 1450-1474 : analyse de correlation MTs-metaux
# L. 1476-1618 : analyse de correlation IC-metaux
# L. 1623-1817 : statistiques multivariees (Analyse discriminante linéaire)
# L. 1818-1932 : statistiques multivariees (Analyse en composante principales)
# L. 1934-2015 : analyses metaux sediments - metaux moules
# L. 2020-2049 : analyse de redondance (test)
# L. 2050-2067 : regressions lineaires

#~~~~~~~~~~~~~~~~~~~~~~~~~ telechargement des packages ~~~~~~~~~~~~~~~~~~~~~~~~~

library("faraway")
library("ggpubr")
library("ISwR")
library("rstatix")
library("tidyverse")
library("dplyr")
library("vegan")
library("car")
library("cowplot")
library("hrbrthemes")
library("viridis")
library("ggplot2")
library("MASS")
library("FactoMineR")
library("factoextra")
library("wesanderson")

#~~~~~~~~~~~~~~~~~~~~~~~~~~ telechargement des donnees ~~~~~~~~~~~~~~~~~~~~~~~~~

# metaux chair
data.Met.ch <- na.omit(read.csv("1_Metaux_chair.csv", sep = ";")) # charge les donnees (ignore les manquantes)
# metaux coquille
data.Met.coq <- na.omit(read.csv("2_Metaux_coquille.csv", sep = ";")) # charge les donnees (ignore les manquantes)
# metallothionenes
data.MTs <- na.omit(read.csv("3_Metallothioneines.csv", sep = ";")) # charge les donnees (ignore les manquantes)
# indices de conditions
data.morpho <- na.omit(read.csv("4_IndicesConditions.csv", sep = ";")) # charge les donnees (ignore les manquantes)
# metaux dans les sediments
sedi <- read.csv("5_sediment.csv",sep=";")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~ indice conditions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~ nettoyage des donnees ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#IDENTIFIE LES OUTLIERS EXTREMES ET LES GARDES DANS UN DF
out <- as.data.frame(data.morpho %>% group_by(comp_site) %>% identify_outliers(longeur..cm.))
ext <- subset(out,out$is.extreme=="TRUE")
out <- data.morpho %>% group_by(comp_site) %>% identify_outliers(largeur..cm.)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.morpho %>% group_by(comp_site) %>% identify_outliers(epaisseur..cm.)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.morpho %>% group_by(comp_site) %>% identify_outliers(poid.totale..mg.)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.morpho %>% group_by(comp_site) %>% identify_outliers(poids.eaux..mg.)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.morpho %>% group_by(comp_site) %>% identify_outliers(poids.tissus..mg.)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.morpho %>% group_by(comp_site) %>% identify_outliers(poids.coquille..mg.)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))

ext <- ext[!duplicated(ext$individu), ] # supprime les doublons

ind <- matrix(ncol=1,nrow=length(ext[,1])) # trouve les index
for( i in 1:length(ext[,1])){
  ind[i] <- which(data.morpho$individu==ext$individu[i])
}

rm(ext,out,i) # nettoie l'environnement

#STATISTIQUES DESCRIPTIVES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dt <-data.morpho[-c(ind),] %>% group_by(comp_site) %>% get_summary_stats(longeur..cm., type = "mean_sd")
mean <- t(round(dt$mean,2))
meanSD <- t(round(dt$sd,2))
dt <-data.morpho[-c(ind),] %>% group_by(comp_site) %>% get_summary_stats(largeur..cm., type = "mean_sd")
mean <- rbind(mean,round(t(dt$mean),2))
meanSD <- rbind(meanSD,round(t(dt$sd),2))
dt <-data.morpho[-c(ind),] %>% group_by(comp_site) %>% get_summary_stats(epaisseur..cm., type = "mean_sd")
mean <- rbind(mean,round(t(dt$mean),2))
meanSD <- rbind(meanSD,round(t(dt$sd),2))
dt <-data.morpho[-c(ind),] %>% group_by(comp_site) %>% get_summary_stats(poid.totale..mg., type = "mean_sd")
mean <- rbind(mean,round(t(dt$mean),2))
meanSD <- rbind(meanSD,round(t(dt$sd),2))
dt <-data.morpho[-c(ind),] %>% group_by(comp_site) %>% get_summary_stats(poids.eaux..mg., type = "mean_sd")
mean <- rbind(mean,round(t(dt$mean),2))
meanSD <- rbind(meanSD,round(t(dt$sd),2))
dt <-data.morpho[-c(ind),] %>% group_by(comp_site) %>% get_summary_stats(poids.tissus..mg., type = "mean_sd")
mean <- rbind(mean,round(t(dt$mean),2))
meanSD <- rbind(meanSD,round(t(dt$sd),2))
dt <-data.morpho[-c(ind),] %>% group_by(comp_site) %>% get_summary_stats(poids.coquille..mg., type = "mean_sd")
mean <- rbind(mean,round(t(dt$mean),2))
meanSD <- rbind(meanSD,round(t(dt$sd),2))
colnames(mean) <- c("VD_1","SP_1","VD_2")
colnames(meanSD) <- c("VD_1","SP_1","VD_2")
rownames(mean) <- c("Longeur [cm]", "Largeur [cm]","Épaisseur [cm]","Poids total [mg]",
                    "Poids eau [mg]","Poids tissu [mg]","Poids coquille [mg]")
rownames(meanSD) <- c("Longeur [cm]", "Largeur [cm]","Épaisseur [cm]","Poids total [mg]",
                    "Poids eau [mg]","Poids tissu [mg]","Poids coquille [mg]")
write.csv(mean,file="moyenne_IC.csv")
write.csv(meanSD,file="moyenneSD_IC.csv")
rm(dt,mean,meanSD)


#TEST DE NORMALITE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
# H0 : les donnees ont une distributions normale
# H1 : les donnes n'ont pas une distribution normale
#la p-value doit etre supérieure à 0.05
data.morpho[-c(ind),] %>% group_by(comp_site) %>% shapiro_test(longeur..cm.)
data.morpho[-c(ind),] %>% group_by(comp_site) %>% shapiro_test(largeur..cm.)
data.morpho[-c(ind),] %>% group_by(comp_site) %>% shapiro_test(epaisseur..cm.)
data.morpho[-c(ind),] %>% group_by(comp_site) %>% shapiro_test(poid.totale..mg.)
data.morpho[-c(ind),] %>% group_by(comp_site) %>% shapiro_test(poids.eaux..mg.)
data.morpho[-c(ind),] %>% group_by(comp_site) %>% shapiro_test(poids.tissus..mg.)
data.morpho[-c(ind),] %>% group_by(comp_site) %>% shapiro_test(poids.coquille..mg.)

# La normalite n'est pas assuree pour toutes les donnees. Il faut donc proceder 
# a un test de Kruskal-Wallis (alternative non-parametrique de l'ANOVA).

#KRUSKAL-WILLIS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#H₀ : toutes les medianes sont egales
#H₁ : au moins une mediane est differente
#seuil de significtion : 0.05

res.kruskal1 <- data.morpho[-c(ind),] %>% kruskal_test(longeur..cm. ~ comp_site)
res.kruskal2 <- data.morpho[-c(ind),] %>% kruskal_test(largeur..cm. ~ comp_site)
res.kruskal3 <- data.morpho[-c(ind),] %>% kruskal_test(epaisseur..cm. ~ comp_site)
res.kruskal4 <- data.morpho[-c(ind),] %>% kruskal_test(poid.totale..mg. ~ comp_site)
res.kruskal5 <- data.morpho[-c(ind),] %>% kruskal_test(poids.eaux..mg.~ comp_site)
res.kruskal6 <- data.morpho[-c(ind),] %>% kruskal_test(poids.tissus..mg. ~ comp_site)
res.kruskal7 <- data.morpho[-c(ind),] %>% kruskal_test(poids.coquille..mg.~ comp_site)

#le test de Kruskal-Willis nous indique q'au moins un groupe est significativement
# differente des autres. Procedons a un test de Dunn 

#DUNN~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data.morpho <- data.morpho %>% # ordone les donnees d'une certaine maniere
  reorder_levels(comp_site, order = c("VD_A", "SP_A", "VD_H"))

pwc1 <- data.morpho[-c(ind),] %>% dunn_test(longeur..cm. ~ campagne, p.adjust.method = "bonferroni") 
pwc1 <- pwc1 %>% add_xy_position(x = "campagne")
gg1 <- ggboxplot(data.morpho[-c(ind),], x = "campagne", y = "longeur..cm.", fill="comp_site") +
  stat_pvalue_manual(pwc1, hide.ns = TRUE)  +
  labs(subtitle = get_test_label(res.kruskal1, detailed = TRUE),caption = get_pwc_label(pwc1)) +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold")) +
  ylab("Longeur individu [cm]") + ggtitle("Longeur") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc2 <- data.morpho[-c(ind),] %>% dunn_test(largeur..cm. ~ campagne, p.adjust.method = "bonferroni") 
pwc2 <- pwc2 %>% add_xy_position(x = "campagne")
gg2 <- ggboxplot(data.morpho[-c(ind),], x = "campagne", y = "largeur..cm.", fill="comp_site") +
  stat_pvalue_manual(pwc2, hide.ns = TRUE)  +
  labs(subtitle = get_test_label(res.kruskal2, detailed = TRUE),caption = get_pwc_label(pwc2)) +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold")) +
  ylab("Largeur individu [cm]") + ggtitle("Largeur") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc3 <- data.morpho[-c(ind),] %>% dunn_test(epaisseur..cm. ~ campagne, p.adjust.method = "bonferroni") 
pwc3 <- pwc3 %>% add_xy_position(x = "campagne")
gg3 <- ggboxplot(data.morpho[-c(ind),], x = "campagne", y = "epaisseur..cm.", fill="comp_site") +
  stat_pvalue_manual(pwc3, hide.ns = TRUE)  +
  labs(subtitle = get_test_label(res.kruskal3, detailed = TRUE),caption = get_pwc_label(pwc3)) +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold")) +
  ylab("Épaisseur individu [cm]") + ggtitle("Épaisseur") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc4 <- data.morpho[-c(ind),] %>% dunn_test(poid.totale..mg. ~ campagne, p.adjust.method = "bonferroni") 
pwc4 <- pwc4 %>% add_xy_position(x = "campagne")
gg4 <- ggboxplot(data.morpho[-c(ind),], x = "campagne", y = "poid.totale..mg.", fill="comp_site") +
  stat_pvalue_manual(pwc4, hide.ns = TRUE)  +
  labs(subtitle = get_test_label(res.kruskal4, detailed = TRUE),caption = get_pwc_label(pwc4)) +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold")) +
  ylab("Poids total individu [mg]") + ggtitle("Masse total") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc5 <- data.morpho[-c(ind),] %>% dunn_test(poids.eaux..mg. ~ campagne, p.adjust.method = "bonferroni") 
pwc5 <- pwc5 %>% add_xy_position(x = "campagne")
gg5 <- ggboxplot(data.morpho[-c(ind),], x = "campagne", y = "poids.eaux..mg.", fill="comp_site") +
  stat_pvalue_manual(pwc5, hide.ns = TRUE)  +
  labs(subtitle = get_test_label(res.kruskal5, detailed = TRUE),caption = get_pwc_label(pwc5)) +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold")) +
  ylab("Poids eau [mg]") + ggtitle("Masse de l'eau") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc6 <- data.morpho[-c(ind),] %>% dunn_test(poids.tissus..mg. ~ campagne, p.adjust.method = "bonferroni") 
pwc6 <- pwc6 %>% add_xy_position(x = "campagne")
gg6 <- ggboxplot(data.morpho[-c(ind),], x = "campagne", y = "poids.tissus..mg.", fill="comp_site") +
  stat_pvalue_manual(pwc6, hide.ns = TRUE)  +
  labs(subtitle = get_test_label(res.kruskal6, detailed = TRUE),caption = get_pwc_label(pwc6)) +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold")) +
  ylab("Poids chair [mg]") + ggtitle("Masse chair") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc7 <- data.morpho[-c(ind),] %>% dunn_test(poids.coquille..mg. ~ campagne, p.adjust.method = "bonferroni") 
pwc7 <- pwc7 %>% add_xy_position(x = "campagne")
gg7 <- ggboxplot(data.morpho[-c(ind),], x = "campagne", y = "poids.coquille..mg.", fill="comp_site") +
  stat_pvalue_manual(pwc7, hide.ns = TRUE)  +
  labs(subtitle = get_test_label(res.kruskal7, detailed = TRUE),caption = get_pwc_label(pwc7)) +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5),
        plot.subtitle = element_text(size = 14),
        plot.caption = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold")) +
  ylab("Poids coquille [mg]") + ggtitle("Masse coquille") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

tit <- ggdraw() + draw_label("Indices de conditions", fontface='bold', size = 25)
tot <- plot_grid(gg1,gg2,gg3,gg4,gg6,gg7,gg5, nrow = 4, ncol = 2, labels = "AUTO", 
                 label_size = 12, align = "v")
plot_grid(tit,tot, ncol=1, rel_heights=c(0.1, 1))
rm(gg1,gg2,gg3,gg4,gg5,gg6,gg7,pwc1,pwc2,pwc3,pwc4,pwc5,pwc6,pwc7,tit,tot,
   res.kruskal1,res.kruskal2,res.kruskal3,res.kruskal4,res.kruskal5,res.kruskal6,res.kruskal7)


data.morpho[-c(ind),] %>% kruskal_effsize(longeur..cm. ~ comp_site)
data.morpho[-c(ind),] %>% kruskal_effsize(largeur..cm. ~ comp_site)
data.morpho[-c(ind),] %>% kruskal_effsize(epaisseur..cm. ~ comp_site)
data.morpho[-c(ind),] %>% kruskal_effsize(poid.totale..mg. ~ comp_site)
data.morpho[-c(ind),] %>% kruskal_effsize(poids.eaux..mg.~ comp_site)
data.morpho[-c(ind),] %>% kruskal_effsize(poids.tissus..mg. ~ comp_site)
data.morpho[-c(ind),] %>% kruskal_effsize(poids.coquille..mg.~ comp_site)

rm(ind)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MTs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~ nettoyage des donnees ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#IDENTIFIE LES OUTLIERS EXTREMES ET LES SUPPRIME
out <- data.MTs %>% group_by(comp_site) %>% identify_outliers(MTs..ug.g.)
ext <- subset(out,out$is.extreme=="TRUE")

ext <- ext[!duplicated(ext$individu), ] # supprime les doublons

ind <- which(data.MTs$individu==ext$individu)

rm(ext,out) # nettoie l'environnement

data.MTs$comp_site <- as.factor(data.MTs$comp_site) # petite adaptation pour le test levene

#TEST DE NORMALITE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
#la p-value doit etre supérieure à 0.05
data.MTs[-c(ind),] %>% group_by(comp_site) %>% shapiro_test(MTs..ug.g.)

# les données sont normalement distribuées (p-value > 0.05)

#TEST D'HOMOGENEITE DES VARIANCES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
#la p-value doit etre supérieure à 0.05
data.MTs[-c(ind),] %>% levene_test(MTs..ug.g. ~ comp_site)

# Il n'y a pas de differences significatives entre les variances d’un groupe a l’autre
# Par conséquent, nous pouvons supposer l’homogénéité des variances dans les différents groupes de traitement.

#TEST ANOVA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

data.MTs <- data.MTs %>% # ordone les donnees d'une certaine maniere
  reorder_levels(comp_site, order = c("VD_A", "SP_A", "VD_H"))

# TEST ANOVA
res.aov <- data.MTs[-c(ind),] %>% anova_test(MTs..ug.g. ~ comp_site)
res.aov

# Test post-hoc
pwc <- data.MTs[-c(ind),] %>% tukey_hsd(MTs..ug.g. ~ comp_site)
pwc

# plot
pwc <- pwc %>% add_xy_position(x = "comp_site")
ggboxplot(data.MTs[-c(ind),], x = "comp_site", y = "MTs..ug.g.", fill="comp_site") +
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(subtitle = get_test_label(res.aov, detailed = TRUE),caption = get_pwc_label(pwc)) +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold")) +
  ylab("Métallothionéines [μg/g]") + ggtitle("Métallothionéine") + xlab("Campagne")

rm(pwc,res.aov,ind) # nettoye l'environnement

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Metaux ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# visualisation des donnees

Met.Sums.ch <- colSums(data.Met.ch[,14:42])
Met.SUM.ch <- sum(Met.Sums.ch)
Class.ch <- sort((Met.Sums.ch/Met.SUM.ch)*100, decreasing = TRUE)

Met.Sums.coq <- colSums(data.Met.coq[,14:42])
Met.SUM.coq <- sum(Met.Sums.coq)
Class.coq <- sort((Met.Sums.coq/Met.SUM.coq)*100, decreasing = TRUE)

par(mfrow=c(2,1))
bar1 <- barplot(Class.ch[1:29], main="Pourcentages des métaux dans la chair",
        ylab = "Concentration  [%]", ylim = c(0,40))
text(x=bar1[,1],
     y=Class.ch+2,
     labels = round(Class.ch,2))
bar2 <- barplot(Class.coq[1:29], main="Pourcentages des métaux dans la coquille",
        ylab = "Concentration  [%]", ylim = c(0,119))
text(x=bar2[,1],
     y=Class.coq+6,
     labels = round(Class.coq,2))

rm(Class.ch,Class.coq,Met.SUM.ch,Met.SUM.coq,Met.Sums.ch,Met.Sums.coq,bar1,bar2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~ metaux independant chair ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~ nettoyage des donnees ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#IDENTIFIE LES OUTLIERS EXTREMES ET LES SUPPRIME
out <- data.Met.ch %>% group_by(comp_site) %>% identify_outliers(Ag)
ext <- subset(out,out$is.extreme=="TRUE")
out <- data.Met.ch %>% group_by(comp_site) %>% identify_outliers(Al)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.ch %>% group_by(comp_site) %>% identify_outliers(As)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.ch %>% group_by(comp_site) %>% identify_outliers(Ba)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.ch %>% group_by(comp_site) %>% identify_outliers(Be)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.ch %>% group_by(comp_site) %>% identify_outliers(Ca)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.ch %>% group_by(comp_site) %>% identify_outliers(Cd)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.ch %>% group_by(comp_site) %>% identify_outliers(Co)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.ch %>% group_by(comp_site) %>% identify_outliers(Cr)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.ch %>% group_by(comp_site) %>% identify_outliers(Cu)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.ch %>% group_by(comp_site) %>% identify_outliers(Fe)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.ch %>% group_by(comp_site) %>% identify_outliers(Hg)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.ch %>% group_by(comp_site) %>% identify_outliers(K)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.ch %>% group_by(comp_site) %>% identify_outliers(Mg)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.ch %>% group_by(comp_site) %>% identify_outliers(Mn)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.ch %>% group_by(comp_site) %>% identify_outliers(Mo)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.ch %>% group_by(comp_site) %>% identify_outliers(Na)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.ch %>% group_by(comp_site) %>% identify_outliers(Ni)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.ch %>% group_by(comp_site) %>% identify_outliers(P)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.ch %>% group_by(comp_site) %>% identify_outliers(Pb)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.ch %>% group_by(comp_site) %>% identify_outliers(S)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.ch %>% group_by(comp_site) %>% identify_outliers(Sb)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.ch %>% group_by(comp_site) %>% identify_outliers(Se)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.ch %>% group_by(comp_site) %>% identify_outliers(Si)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.ch %>% group_by(comp_site) %>% identify_outliers(Sn)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.ch %>% group_by(comp_site) %>% identify_outliers(Sr)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.ch %>% group_by(comp_site) %>% identify_outliers(Ti)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.ch %>% group_by(comp_site) %>% identify_outliers(V)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.ch %>% group_by(comp_site) %>% identify_outliers(Zn)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))

ext <- ext[!duplicated(ext$pool), ] # supprime les doublons

ind <- matrix(ncol=1,nrow=10) # trouve les index
for( i in 1:10){
  ind[i] <- which(data.Met.ch$pool==ext$pool[i])
}

rm(ext,out,i) # nettoie l'environnement


#KRUSKAL-WILLIS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#H₀ : toutes les medianes sont egales
#H₁ : au moins une mediane est differente
#seuil de significtion : 0.05

data.Met.ch[-c(ind),] %>% kruskal_test(Ag ~ comp_site)
data.Met.ch[-c(ind),] %>% kruskal_test(Al ~ comp_site)
data.Met.ch[-c(ind),] %>% kruskal_test(As ~ comp_site)
data.Met.ch[-c(ind),] %>% kruskal_test(Ba ~ comp_site)
data.Met.ch[-c(ind),] %>% kruskal_test(Be ~ comp_site)
data.Met.ch[-c(ind),] %>% kruskal_test(Ca ~ comp_site) # pas de diff
data.Met.ch[-c(ind),] %>% kruskal_test(Cd ~ comp_site)
data.Met.ch[-c(ind),] %>% kruskal_test(Co ~ comp_site)
data.Met.ch[-c(ind),] %>% kruskal_test(Cr ~ comp_site)
data.Met.ch[-c(ind),] %>% kruskal_test(Cu ~ comp_site)
data.Met.ch[-c(ind),] %>% kruskal_test(Fe ~ comp_site)
data.Met.ch[-c(ind),] %>% kruskal_test(Hg ~ comp_site)
data.Met.ch[-c(ind),] %>% kruskal_test(K ~ comp_site)
data.Met.ch[-c(ind),] %>% kruskal_test(Mg ~ comp_site)
data.Met.ch[-c(ind),] %>% kruskal_test(Mn ~ comp_site)
data.Met.ch[-c(ind),] %>% kruskal_test(Mo ~ comp_site)
data.Met.ch[-c(ind),] %>% kruskal_test(Na ~ comp_site)  # pas de diff
data.Met.ch[-c(ind),] %>% kruskal_test(Ni ~ comp_site)
data.Met.ch[-c(ind),] %>% kruskal_test(Pb ~ comp_site)
data.Met.ch[-c(ind),] %>% kruskal_test(S ~ comp_site)
data.Met.ch[-c(ind),] %>% kruskal_test(Sb ~ comp_site) # pas de diff
data.Met.ch[-c(ind),] %>% kruskal_test(Se ~ comp_site) # pas de diff
data.Met.ch[-c(ind),] %>% kruskal_test(Si ~ comp_site)
data.Met.ch[-c(ind),] %>% kruskal_test(Sn ~ comp_site)
data.Met.ch[-c(ind),] %>% kruskal_test(Sr ~ comp_site)
data.Met.ch[-c(ind),] %>% kruskal_test(Ti ~ comp_site)
data.Met.ch[-c(ind),]%>% kruskal_test(V ~ comp_site)
data.Met.ch[-c(ind),] %>% kruskal_test(Zn ~ comp_site)

#le test de Kruskal-Willis nous indique q'au moins un groupe est significativement
# differente des autres. Procedons a un test de Dunn et wilcox pour trouver ces differences.

#DUNN~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

data.Met.ch <- data.Met.ch %>% # ordone les donnees d'une certaine maniere
  reorder_levels(comp_site, order = c("VD_A", "SP_A", "VD_H"))

pwc1 <- data.Met.ch[-c(ind),] %>% dunn_test(Ag ~ campagne, p.adjust.method = "bonferroni") 
pwc1 <- pwc1 %>% add_xy_position(x = "campagne")
gg1 <- ggboxplot(data.Met.ch[-c(ind),], x = "campagne", y = "Ag", fill="comp_site") +
  stat_pvalue_manual(pwc1, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Argent") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc2 <- data.Met.ch[-c(ind),] %>% dunn_test( Al ~ campagne, p.adjust.method = "bonferroni") 
pwc2 <- pwc2 %>% add_xy_position(x = "campagne")
gg2 <- ggboxplot(data.Met.ch[-c(ind),], x = "campagne", y = "Al", fill="comp_site") +
  stat_pvalue_manual(pwc2, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Aluminium") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc3 <- data.Met.ch[-c(ind),] %>% dunn_test( As ~ campagne, p.adjust.method = "bonferroni") 
pwc3 <- pwc3 %>% add_xy_position(x = "campagne")
gg3 <- ggboxplot(data.Met.ch[-c(ind),], x = "campagne", y = "As", fill="comp_site") +
  stat_pvalue_manual(pwc3, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Arsenic") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc4 <- data.Met.ch[-c(ind),] %>% dunn_test( Ba ~ campagne, p.adjust.method = "bonferroni") 
pwc4 <- pwc4 %>% add_xy_position(x = "campagne")
gg4 <- ggboxplot(data.Met.ch[-c(ind),], x = "campagne", y = "Ba", fill="comp_site") +
  stat_pvalue_manual(pwc4, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Barium") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc5 <- data.Met.ch[-c(ind),] %>% dunn_test( Be ~ campagne, p.adjust.method = "bonferroni") 
pwc5 <- pwc5 %>% add_xy_position(x = "campagne")
gg5 <- ggboxplot(data.Met.ch[-c(ind),], x = "campagne", y = "Be", fill="comp_site") +
  stat_pvalue_manual(pwc5, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Bérillyum") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc6 <- data.Met.ch[-c(ind),] %>% dunn_test( Ca ~ campagne, p.adjust.method = "bonferroni") 
pwc6 <- pwc6 %>% add_xy_position(x = "campagne")
gg6 <- ggboxplot(data.Met.ch[-c(ind),], x = "campagne", y = "Ca", fill="comp_site") +
  stat_pvalue_manual(pwc6, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Calcium") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc7 <- data.Met.ch[-c(ind),] %>% dunn_test( Cd ~ campagne, p.adjust.method = "bonferroni") 
pwc7 <- pwc7 %>% add_xy_position(x = "campagne")
gg7 <- ggboxplot(data.Met.ch[-c(ind),], x = "campagne", y = "Cd", fill="comp_site") +
  stat_pvalue_manual(pwc7, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Cadmium") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc8 <- data.Met.ch[-c(ind),] %>% dunn_test( Co ~ campagne, p.adjust.method = "bonferroni") 
pwc8 <- pwc8 %>% add_xy_position(x = "campagne")
gg8 <- ggboxplot(data.Met.ch[-c(ind),], x = "campagne", y = "Co", fill="comp_site") +
  stat_pvalue_manual(pwc8, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Cobalt") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc9 <- data.Met.ch[-c(ind),] %>% dunn_test( Cr ~ campagne, p.adjust.method = "bonferroni") 
pwc9 <- pwc9 %>% add_xy_position(x = "campagne")
gg9 <- ggboxplot(data.Met.ch[-c(ind),], x = "campagne", y = "Cr", fill="comp_site") +
  stat_pvalue_manual(pwc9, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Chrome") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc10 <- data.Met.ch[-c(ind),] %>% dunn_test( Cu ~ campagne, p.adjust.method = "bonferroni") 
pwc10 <- pwc10 %>% add_xy_position(x = "campagne")
gg10 <- ggboxplot(data.Met.ch[-c(ind),], x = "campagne", y = "Cu", fill="comp_site") +
  stat_pvalue_manual(pwc10, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Cuivre") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc11 <- data.Met.ch[-c(ind),] %>% dunn_test( Fe ~ campagne, p.adjust.method = "bonferroni") 
pwc11 <- pwc11 %>% add_xy_position(x = "campagne")
gg11 <- ggboxplot(data.Met.ch[-c(ind),], x = "campagne", y = "Fe", fill="comp_site") +
  stat_pvalue_manual(pwc11, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Fer") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc12 <- data.Met.ch[-c(ind),] %>% dunn_test( Hg ~ campagne, p.adjust.method = "bonferroni") 
pwc12 <- pwc12 %>% add_xy_position(x = "campagne")
gg12 <- ggboxplot(data.Met.ch[-c(ind),], x = "campagne", y = "Hg", fill="comp_site") +
  stat_pvalue_manual(pwc12, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Mercure") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc13 <- data.Met.ch[-c(ind),] %>% dunn_test( K ~ campagne, p.adjust.method = "bonferroni") 
pwc13 <- pwc13 %>% add_xy_position(x = "campagne")
gg13 <- ggboxplot(data.Met.ch[-c(ind),], x = "campagne", y = "K", fill="comp_site") +
  stat_pvalue_manual(pwc13, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Potassium") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc14 <- data.Met.ch[-c(ind),] %>% dunn_test( Mg ~ campagne, p.adjust.method = "bonferroni") 
pwc14 <- pwc14 %>% add_xy_position(x = "campagne")
gg14 <- ggboxplot(data.Met.ch[-c(ind),], x = "campagne", y = "Mg", fill="comp_site") +
  stat_pvalue_manual(pwc14, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Magnésium") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc15 <- data.Met.ch[-c(ind),] %>% dunn_test( Mn ~ campagne, p.adjust.method = "bonferroni") 
pwc15 <- pwc15 %>% add_xy_position(x = "campagne")
gg15 <- ggboxplot(data.Met.ch[-c(ind),], x = "campagne", y = "Mn", fill="comp_site") +
  stat_pvalue_manual(pwc15, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Manganèse") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc16 <- data.Met.ch[-c(ind),] %>% dunn_test( Mo ~ campagne, p.adjust.method = "bonferroni") 
pwc16 <- pwc16 %>% add_xy_position(x = "campagne")
gg16 <- ggboxplot(data.Met.ch[-c(ind),], x = "campagne", y = "Mo", fill="comp_site") +
  stat_pvalue_manual(pwc16, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Molybdène") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc17 <- data.Met.ch[-c(ind),] %>% dunn_test(Na ~ campagne, p.adjust.method = "bonferroni") 
pwc17 <- pwc17 %>% add_xy_position(x = "campagne")
gg17 <- ggboxplot(data.Met.ch[-c(ind),], x = "campagne", y = "Na", fill="comp_site") +
  stat_pvalue_manual(pwc17, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Sodium") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc18 <- data.Met.ch[-c(ind),] %>% dunn_test(Ni ~ campagne, p.adjust.method = "bonferroni") 
pwc18 <- pwc18 %>% add_xy_position(x = "campagne")
gg18 <- ggboxplot(data.Met.ch[-c(ind),], x = "campagne", y = "Ni", fill="comp_site") +
  stat_pvalue_manual(pwc18, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Nickel") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc19 <- data.Met.ch[-c(ind),] %>% dunn_test(P ~ campagne, p.adjust.method = "bonferroni") 
pwc19 <- pwc19 %>% add_xy_position(x = "campagne")
gg19 <- ggboxplot(data.Met.ch[-c(ind),], x = "campagne", y = "P", fill="comp_site") +
  stat_pvalue_manual(pwc19, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Phosphore") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc20 <- data.Met.ch[-c(ind),] %>% dunn_test(Pb ~ campagne, p.adjust.method = "bonferroni") 
pwc20 <- pwc20 %>% add_xy_position(x = "campagne")
gg20 <- ggboxplot(data.Met.ch[-c(ind),], x = "campagne", y = "Pb", fill="comp_site") +
  stat_pvalue_manual(pwc20, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Plomb") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc21 <- data.Met.ch[-c(ind),] %>% dunn_test(S ~ campagne, p.adjust.method = "bonferroni") 
pwc21 <- pwc21 %>% add_xy_position(x = "campagne")
gg21 <- ggboxplot(data.Met.ch[-c(ind),], x = "campagne", y = "S", fill="comp_site") +
  stat_pvalue_manual(pwc21, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Soufre") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc22 <- data.Met.ch[-c(ind),] %>% dunn_test(Sb ~ campagne, p.adjust.method = "bonferroni") 
pwc22 <- pwc22 %>% add_xy_position(x = "campagne")
gg22 <- ggboxplot(data.Met.ch[-c(ind),], x = "campagne", y = "Sb", fill="comp_site") +
  stat_pvalue_manual(pwc22, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Antimoine") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc23 <- data.Met.ch[-c(ind),] %>% dunn_test(Se ~ campagne, p.adjust.method = "bonferroni") 
pwc23 <- pwc23 %>% add_xy_position(x = "campagne")
gg23 <- ggboxplot(data.Met.ch[-c(ind),], x = "campagne", y = "Se", fill="comp_site") +
  stat_pvalue_manual(pwc23, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Sélénium") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc24 <- data.Met.ch[-c(ind),] %>% dunn_test(Si ~ campagne, p.adjust.method = "bonferroni") 
pwc24 <- pwc24 %>% add_xy_position(x = "campagne")
gg24 <- ggboxplot(data.Met.ch[-c(ind),], x = "campagne", y = "Si", fill="comp_site") +
  stat_pvalue_manual(pwc24, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Silicium") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc25 <- data.Met.ch[-c(ind),] %>% dunn_test(Sn ~ campagne, p.adjust.method = "bonferroni") 
pwc25 <- pwc25 %>% add_xy_position(x = "campagne")
gg25 <- ggboxplot(data.Met.ch[-c(ind),], x = "campagne", y = "Sn", fill="comp_site") +
  stat_pvalue_manual(pwc25, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Étain") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc26 <- data.Met.ch[-c(ind),] %>% dunn_test(Sr ~ campagne, p.adjust.method = "bonferroni") 
pwc26 <- pwc26 %>% add_xy_position(x = "campagne")
gg26 <- ggboxplot(data.Met.ch[-c(ind),], x = "campagne", y = "Sr", fill="comp_site") +
  stat_pvalue_manual(pwc26, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Strontium") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc27 <- data.Met.ch[-c(ind),] %>% dunn_test(Ti ~ campagne, p.adjust.method = "bonferroni") 
pwc27 <- pwc27 %>% add_xy_position(x = "campagne")
gg27 <- ggboxplot(data.Met.ch[-c(ind),], x = "campagne", y = "Ti", fill="comp_site") +
  stat_pvalue_manual(pwc27, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Titane") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc28 <- data.Met.ch[-c(ind),] %>% dunn_test(V ~ campagne, p.adjust.method = "bonferroni") 
pwc28 <- pwc28 %>% add_xy_position(x = "campagne")
gg28 <- ggboxplot(data.Met.ch[-c(ind),], x = "campagne", y = "V", fill="comp_site") +
  stat_pvalue_manual(pwc28, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Vanadium") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc29 <- data.Met.ch[-c(ind),] %>% dunn_test(Zn ~ campagne, p.adjust.method = "bonferroni") 
pwc29 <- pwc29 %>% add_xy_position(x = "campagne")
gg29 <- ggboxplot(data.Met.ch[-c(ind),], x = "campagne", y = "Zn", fill="comp_site") +
  stat_pvalue_manual(pwc29, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Zinc") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

plot_grid(gg1,gg2,gg3,gg4,gg5,gg6,gg7,gg8,gg9,gg10,gg11,gg12,gg13,gg14,gg15,gg16,
          gg17,gg18,gg19,gg20,gg21,gg22,gg23,gg24,gg25,gg26,gg27,gg28,gg29, 
          nrow = 6, ncol = 5, labels = "AUTO",label_size = 12, align = "v")

rm(gg1,gg2,gg3,gg4,gg5,gg6,gg7,gg8,gg9,gg10,gg11,gg12,gg13,gg14,gg15,gg16,
   gg17,gg18,gg19,gg20,gg21,gg22,gg23,gg24,gg25,gg26,gg27,gg28,gg29,
   pwc1,pwc2,pwc3,pwc4,pwc5,pwc6,pwc7,pwc8,pwc9,pwc10,pwc11,pwc12,pwc13,pwc14,pwc15,pwc16,
   pwc17,pwc18,pwc19,pwc20,pwc21,pwc22,pwc23,pwc24,pwc25,pwc26,pwc27,pwc28,pwc29,ind)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~ metaux independant coq ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~ nettoyage des donnees ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#IDENTIFIE LES OUTLIERS EXTREMES ET LES SUPPRIME
out <- data.Met.coq %>% group_by(comp_site) %>% identify_outliers(Ag)
ext <- subset(out,out$is.extreme=="TRUE")
out <- data.Met.coq %>% group_by(comp_site) %>% identify_outliers(Al)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.coq %>% group_by(comp_site) %>% identify_outliers(As)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.coq %>% group_by(comp_site) %>% identify_outliers(Ba)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.coq %>% group_by(comp_site) %>% identify_outliers(Be)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.coq %>% group_by(comp_site) %>% identify_outliers(Ca)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.coq %>% group_by(comp_site) %>% identify_outliers(Cd)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.coq %>% group_by(comp_site) %>% identify_outliers(Co)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.coq %>% group_by(comp_site) %>% identify_outliers(Cr)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.coq %>% group_by(comp_site) %>% identify_outliers(Cu)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.coq %>% group_by(comp_site) %>% identify_outliers(Fe)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.coq %>% group_by(comp_site) %>% identify_outliers(Hg)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.coq %>% group_by(comp_site) %>% identify_outliers(K)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.coq %>% group_by(comp_site) %>% identify_outliers(Mg)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.coq %>% group_by(comp_site) %>% identify_outliers(Mn)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.coq %>% group_by(comp_site) %>% identify_outliers(Mo)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.coq %>% group_by(comp_site) %>% identify_outliers(Na)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.coq %>% group_by(comp_site) %>% identify_outliers(Ni)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.coq %>% group_by(comp_site) %>% identify_outliers(P)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.coq %>% group_by(comp_site) %>% identify_outliers(Pb)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.coq %>% group_by(comp_site) %>% identify_outliers(S)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.coq %>% group_by(comp_site) %>% identify_outliers(Sb)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.coq %>% group_by(comp_site) %>% identify_outliers(Se)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.coq %>% group_by(comp_site) %>% identify_outliers(Si)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.coq %>% group_by(comp_site) %>% identify_outliers(Sn)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.coq %>% group_by(comp_site) %>% identify_outliers(Sr)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.coq %>% group_by(comp_site) %>% identify_outliers(Ti)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.coq %>% group_by(comp_site) %>% identify_outliers(V)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))
out <- data.Met.coq %>% group_by(comp_site) %>% identify_outliers(Zn)
ext <- rbind(ext,subset(out,out$is.extreme=="TRUE"))

ext <- ext[!duplicated(ext$pool), ] # supprime les doublons

ind <- matrix(ncol=1,nrow=6) # trouve les index
for( i in 1:6){
  ind[i] <- which(data.Met.coq$pool==ext$pool[i])
}

rm(ext,out,i) # nettoie l'environnement


#KRUSKAL-WILLIS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#DUNN~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

data.Met.coq <- data.Met.coq %>% # ordone les donnees d'une certaine maniere
  reorder_levels(comp_site, order = c("VD_A", "SP_A", "VD_H"))

pwc1 <- data.Met.coq[-c(ind),] %>% dunn_test(Ag ~ campagne, p.adjust.method = "bonferroni") 
pwc1 <- pwc1 %>% add_xy_position(x = "campagne")
gg1 <- ggboxplot(data.Met.coq[-c(ind),], x = "campagne", y = "Ag", fill="comp_site") +
  stat_pvalue_manual(pwc1, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Argent") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc2 <- data.Met.coq[-c(ind),] %>% dunn_test( Al ~ campagne, p.adjust.method = "bonferroni") 
pwc2 <- pwc2 %>% add_xy_position(x = "campagne")
gg2 <- ggboxplot(data.Met.coq[-c(ind),], x = "campagne", y = "Al", fill="comp_site") +
  stat_pvalue_manual(pwc2, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Aluminium") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc3 <- data.Met.coq[-c(ind),] %>% dunn_test( As ~ campagne, p.adjust.method = "bonferroni") 
pwc3 <- pwc3 %>% add_xy_position(x = "campagne")
gg3 <- ggboxplot(data.Met.coq[-c(ind),], x = "campagne", y = "As", fill="comp_site") +
  stat_pvalue_manual(pwc3, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Arsenic") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc4 <- data.Met.coq[-c(ind),] %>% dunn_test( Ba ~ campagne, p.adjust.method = "bonferroni") 
pwc4 <- pwc4 %>% add_xy_position(x = "campagne")
gg4 <- ggboxplot(data.Met.coq[-c(ind),], x = "campagne", y = "Ba", fill="comp_site") +
  stat_pvalue_manual(pwc4, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Barium") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc5 <- data.Met.coq[-c(ind),] %>% dunn_test( Be ~ campagne, p.adjust.method = "bonferroni") 
pwc5 <- pwc5 %>% add_xy_position(x = "campagne")
gg5 <- ggboxplot(data.Met.coq[-c(ind),], x = "campagne", y = "Be", fill="comp_site") +
  stat_pvalue_manual(pwc5, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Bérillyum") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc6 <- data.Met.coq[-c(ind),] %>% dunn_test( Ca ~ campagne, p.adjust.method = "bonferroni") 
pwc6 <- pwc6 %>% add_xy_position(x = "campagne")
gg6 <- ggboxplot(data.Met.coq[-c(ind),], x = "campagne", y = "Ca", fill="comp_site") +
  stat_pvalue_manual(pwc6, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Calcium") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc7 <- data.Met.coq[-c(ind),] %>% dunn_test( Cd ~ campagne, p.adjust.method = "bonferroni") 
pwc7 <- pwc7 %>% add_xy_position(x = "campagne")
gg7 <- ggboxplot(data.Met.coq[-c(ind),], x = "campagne", y = "Cd", fill="comp_site") +
  stat_pvalue_manual(pwc7, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Cadmium") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc8 <- data.Met.coq[-c(ind),] %>% dunn_test( Co ~ campagne, p.adjust.method = "bonferroni") 
pwc8 <- pwc8 %>% add_xy_position(x = "campagne")
gg8 <- ggboxplot(data.Met.coq[-c(ind),], x = "campagne", y = "Co", fill="comp_site") +
  stat_pvalue_manual(pwc8, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Cobalt") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc9 <- data.Met.coq[-c(ind),] %>% dunn_test( Cr ~ campagne, p.adjust.method = "bonferroni") 
pwc9 <- pwc9 %>% add_xy_position(x = "campagne")
gg9 <- ggboxplot(data.Met.coq[-c(ind),], x = "campagne", y = "Cr", fill="comp_site") +
  stat_pvalue_manual(pwc9, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("chrome") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc10 <- data.Met.coq[-c(ind),] %>% dunn_test( Cu ~ campagne, p.adjust.method = "bonferroni") 
pwc10 <- pwc10 %>% add_xy_position(x = "campagne")
gg10 <- ggboxplot(data.Met.coq[-c(ind),], x = "campagne", y = "Cu", fill="comp_site") +
  stat_pvalue_manual(pwc10, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Cuivre") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc11 <- data.Met.coq[-c(ind),] %>% dunn_test( Fe ~ campagne, p.adjust.method = "bonferroni") 
pwc11 <- pwc11 %>% add_xy_position(x = "campagne")
gg11 <- ggboxplot(data.Met.coq[-c(ind),], x = "campagne", y = "Fe", fill="comp_site") +
  stat_pvalue_manual(pwc11, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Fer") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc12 <- data.Met.coq[-c(ind),] %>% dunn_test( Hg ~ campagne, p.adjust.method = "bonferroni") 
pwc12 <- pwc12 %>% add_xy_position(x = "campagne")
gg12 <- ggboxplot(data.Met.coq[-c(ind),], x = "campagne", y = "Hg", fill="comp_site") +
  stat_pvalue_manual(pwc12, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Mercure") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc13 <- data.Met.coq[-c(ind),] %>% dunn_test( K ~ campagne, p.adjust.method = "bonferroni") 
pwc13 <- pwc13 %>% add_xy_position(x = "campagne")
gg13 <- ggboxplot(data.Met.coq[-c(ind),], x = "campagne", y = "K", fill="comp_site") +
  stat_pvalue_manual(pwc13, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Potassium") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc14 <- data.Met.coq[-c(ind),] %>% dunn_test( Mg ~ campagne, p.adjust.method = "bonferroni") 
pwc14 <- pwc14 %>% add_xy_position(x = "campagne")
gg14 <- ggboxplot(data.Met.coq[-c(ind),], x = "campagne", y = "Mg", fill="comp_site") +
  stat_pvalue_manual(pwc14, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Magnésium") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc15 <- data.Met.coq[-c(ind),] %>% dunn_test( Mn ~ campagne, p.adjust.method = "bonferroni") 
pwc15 <- pwc15 %>% add_xy_position(x = "campagne")
gg15 <- ggboxplot(data.Met.coq[-c(ind),], x = "campagne", y = "Mn", fill="comp_site") +
  stat_pvalue_manual(pwc15, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Manganèse") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc16 <- data.Met.coq[-c(ind),] %>% dunn_test( Mo ~ campagne, p.adjust.method = "bonferroni") 
pwc16 <- pwc16 %>% add_xy_position(x = "campagne")
gg16 <- ggboxplot(data.Met.coq[-c(ind),], x = "campagne", y = "Mo", fill="comp_site") +
  stat_pvalue_manual(pwc16, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Molybdène") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc17 <- data.Met.coq[-c(ind),] %>% dunn_test(Na ~ campagne, p.adjust.method = "bonferroni") 
pwc17 <- pwc17 %>% add_xy_position(x = "campagne")
gg17 <- ggboxplot(data.Met.coq[-c(ind),], x = "campagne", y = "Na", fill="comp_site") +
  stat_pvalue_manual(pwc17, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Sodium") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc18 <- data.Met.coq[-c(ind),] %>% dunn_test(Ni ~ campagne, p.adjust.method = "bonferroni") 
pwc18 <- pwc18 %>% add_xy_position(x = "campagne")
gg18 <- ggboxplot(data.Met.coq[-c(ind),], x = "campagne", y = "Ni", fill="comp_site") +
  stat_pvalue_manual(pwc18, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Nickel") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc19 <- data.Met.coq[-c(ind),] %>% dunn_test(P ~ campagne, p.adjust.method = "bonferroni") 
pwc19 <- pwc19 %>% add_xy_position(x = "campagne")
gg19 <- ggboxplot(data.Met.coq[-c(ind),], x = "campagne", y = "P", fill="comp_site") +
  stat_pvalue_manual(pwc19, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Phosphore") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc20 <- data.Met.coq[-c(ind),] %>% dunn_test(Pb ~ campagne, p.adjust.method = "bonferroni") 
pwc20 <- pwc20 %>% add_xy_position(x = "campagne")
gg20 <- ggboxplot(data.Met.coq[-c(ind),], x = "campagne", y = "Pb", fill="comp_site") +
  stat_pvalue_manual(pwc20, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Plomb") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc21 <- data.Met.coq[-c(ind),] %>% dunn_test(S ~ campagne, p.adjust.method = "bonferroni") 
pwc21 <- pwc21 %>% add_xy_position(x = "campagne")
gg21 <- ggboxplot(data.Met.coq[-c(ind),], x = "campagne", y = "S", fill="comp_site") +
  stat_pvalue_manual(pwc21, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Soufre") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc22 <- data.Met.coq[-c(ind),] %>% dunn_test(Sb ~ campagne, p.adjust.method = "bonferroni") 
pwc22 <- pwc22 %>% add_xy_position(x = "campagne")
gg22 <- ggboxplot(data.Met.coq[-c(ind),], x = "campagne", y = "Sb", fill="comp_site") +
  stat_pvalue_manual(pwc22, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Antimoine") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc23 <- data.Met.coq[-c(ind),] %>% dunn_test(Se ~ campagne, p.adjust.method = "bonferroni") 
pwc23 <- pwc23 %>% add_xy_position(x = "campagne")
gg23 <- ggboxplot(data.Met.coq[-c(ind),], x = "campagne", y = "Se", fill="comp_site") +
  stat_pvalue_manual(pwc23, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Sélénium") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc24 <- data.Met.coq[-c(ind),] %>% dunn_test(Si ~ campagne, p.adjust.method = "bonferroni") 
pwc24 <- pwc24 %>% add_xy_position(x = "campagne")
gg24 <- ggboxplot(data.Met.coq[-c(ind),], x = "campagne", y = "Si", fill="comp_site") +
  stat_pvalue_manual(pwc24, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Silicium") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc25 <- data.Met.coq[-c(ind),] %>% dunn_test(Sn ~ campagne, p.adjust.method = "bonferroni") 
pwc25 <- pwc25 %>% add_xy_position(x = "campagne")
gg25 <- ggboxplot(data.Met.coq[-c(ind),], x = "campagne", y = "Sn", fill="comp_site") +
  stat_pvalue_manual(pwc25, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Étain") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc26 <- data.Met.coq[-c(ind),] %>% dunn_test(Sr ~ campagne, p.adjust.method = "bonferroni") 
pwc26 <- pwc26 %>% add_xy_position(x = "campagne")
gg26 <- ggboxplot(data.Met.coq[-c(ind),], x = "campagne", y = "Sr", fill="comp_site") +
  stat_pvalue_manual(pwc26, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Strontium") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc27 <- data.Met.coq[-c(ind),] %>% dunn_test(Ti ~ campagne, p.adjust.method = "bonferroni") 
pwc27 <- pwc27 %>% add_xy_position(x = "campagne")
gg27 <- ggboxplot(data.Met.coq[-c(ind),], x = "campagne", y = "Ti", fill="comp_site") +
  stat_pvalue_manual(pwc27, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Titane") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc28 <- data.Met.coq[-c(ind),] %>% dunn_test(V ~ campagne, p.adjust.method = "bonferroni") 
pwc28 <- pwc28 %>% add_xy_position(x = "campagne")
gg28 <- ggboxplot(data.Met.coq[-c(ind),], x = "campagne", y = "V", fill="comp_site") +
  stat_pvalue_manual(pwc28, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5)) +
  ylab("Concentration [ppm]") + ggtitle("Vanadium") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc29 <- data.Met.coq[-c(ind),] %>% dunn_test(Zn ~ campagne, p.adjust.method = "bonferroni") 
pwc29 <- pwc29 %>% add_xy_position(x = "campagne")
gg29 <- ggboxplot(data.Met.coq[-c(ind),], x = "campagne", y = "Zn", fill="comp_site") +
  stat_pvalue_manual(pwc29, hide.ns = TRUE)  +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5))+
  ylab("Concentration [ppm]") + ggtitle("Zinc") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

plot_grid(gg1,gg2,gg3,gg4,gg5,gg6,gg7,gg8,gg9,gg10,gg11,gg12,gg13,gg14,gg15,gg16,
          gg17,gg18,gg19,gg20,gg21,gg22,gg23,gg24,gg25,gg26,gg27,gg28,gg29, 
          nrow = 6, ncol = 5, labels = "AUTO",label_size = 12, align = "v")

rm(gg1,gg2,gg3,gg4,gg5,gg6,gg7,gg8,gg9,gg10,gg11,gg12,gg13,gg14,gg15,gg16,
   gg17,gg18,gg19,gg20,gg21,gg22,gg23,gg24,gg25,gg26,gg27,gg28,gg29,
   pwc1,pwc2,pwc3,pwc4,pwc5,pwc6,pwc7,pwc8,pwc9,pwc10,pwc11,pwc12,pwc13,pwc14,pwc15,pwc16,
   pwc17,pwc18,pwc19,pwc20,pwc21,pwc22,pwc23,pwc24,pwc25,pwc26,pwc27,pwc28,pwc29,ind)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# on fait une comparaison pour tous les metaux
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~ nettoyage des donnees ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#IDENTIFIE LES OUTLIERS EXTREMES ET LES GARDES DANS UN DF _ Chair
as.data.frame(data.Met.ch %>% group_by(comp_site) %>% identify_outliers(tot))

#IDENTIFIE LES OUTLIERS EXTREMES ET LES GARDES DANS UN DF _ Coquille
as.data.frame(data.Met.coq %>% group_by(comp_site) %>% identify_outliers(tot))

# pas de outliers extreme pour la chair et la coquille !!!

#TEST DE NORMALITE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
#la p-value doit etre supérieure à 0.05
data.Met.ch %>% group_by(comp_site) %>% shapiro_test(tot)
data.Met.coq %>% group_by(comp_site) %>% shapiro_test(tot)

# les données sont normalement distribuées (p-value > 0.05)

#TEST D'HOMOGENEITE DES VARIANCES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
#la p-value doit etre supérieure à 0.05
data.Met.ch %>% levene_test(tot ~ comp_site)
data.Met.coq%>% levene_test(tot ~ comp_site)

# Il n'y a pas de differences significatives entre les variances d’un groupe a l’autre
# Par conséquent, nous pouvons supposer l’homogénéité des variances dans les différents groupes de traitement.

#TEST ANOVA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

data.Met.ch <- data.Met.ch %>% # ordone les donnees d'une certaine maniere
  reorder_levels(comp_site, order = c("VD_A", "SP_A", "VD_H"))
data.Met.coq <- data.Met.coq %>% # ordone les donnees d'une certaine maniere
  reorder_levels(comp_site, order = c("VD_A", "SP_A", "VD_H"))

# TEST ANOVA
res.aov.ch <- data.Met.ch %>% anova_test(tot ~ comp_site)
res.aov.ch
res.aov.coq <- data.Met.coq %>% anova_test(tot ~ comp_site)
res.aov.coq

# Test post-hoc
pwc.ch <- data.Met.ch %>% tukey_hsd(tot ~ comp_site)
pwc.ch
pwc.coq <- data.Met.coq %>% tukey_hsd(tot ~ comp_site)
pwc.coq

# plot
pwc.ch <- pwc.ch %>% add_xy_position(x = "comp_site")
ggplot.ch.all <- ggboxplot(data.Met.ch, x = "comp_site", y = "tot", fill="comp_site") +
  stat_pvalue_manual(pwc.ch, hide.ns = TRUE) +
  labs(subtitle = get_test_label(res.aov.ch, detailed = TRUE),caption = get_pwc_label(pwc.ch)) +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold")) +
  ylab("Concentration [ppm]") + ggtitle("Métaux totaux - chair") + xlab("Campagne")

pwc.coq <- pwc.coq %>% add_xy_position(x = "comp_site")
ggplot.coq.all <- ggboxplot(data.Met.coq, x = "comp_site", y = "tot", fill="comp_site") +
  stat_pvalue_manual(pwc.coq, hide.ns = TRUE) +
  labs(subtitle = get_test_label(res.aov.coq, detailed = TRUE),caption = get_pwc_label(pwc.coq)) +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold")) +
  ylab("Concentration [ppm]") + ggtitle("Métaux totaux - coquille") + xlab("Campagne")

rm(pwc.ch,pwc.coq,res.aov.ch,res.aov.coq)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# on refait mais cette fois-ci sans le Ca (forte conc.)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

data.Met.ch$totWCa <- data.Met.ch$tot - data.Met.ch$Ca
data.Met.coq$totWCa <- data.Met.coq$tot - data.Met.coq$Ca

#IDENTIFIE LES OUTLIERS EXTREMES ET LES GARDES DANS UN DF _ Chair
as.data.frame(data.Met.ch %>% group_by(comp_site) %>% identify_outliers(totWCa))

#IDENTIFIE LES OUTLIERS EXTREMES ET LES GARDES DANS UN DF _ Coquille
out <- as.data.frame(data.Met.coq %>% group_by(comp_site) %>% identify_outliers(totWCa))
ext <- subset(out,out$is.extreme=="TRUE")

ind <- matrix(ncol=1,nrow=1) # trouve les index
for( i in 1:1){
  ind[i] <- which(data.Met.coq$pool==ext$pool[i])
}

rm(ext,out,i) # nettoie l'environnement

# la normalite et l'homogeneite des variances ont etes testes est c'est OK !

#TEST ANOVA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

data.Met.ch <- data.Met.ch %>% # ordone les donnees d'une certaine maniere
  reorder_levels(comp_site, order = c("VD_A", "SP_A", "VD_H"))
data.Met.coq <- data.Met.coq %>% # ordone les donnees d'une certaine maniere
  reorder_levels(comp_site, order = c("VD_A", "SP_A", "VD_H"))

# TEST ANOVA
res.aov.ch.WA <- data.Met.ch %>% anova_test(totWCa ~ comp_site)
res.aov.ch.WA
res.aov.coq.WA <- data.Met.coq[-c(ind),] %>% anova_test(totWCa ~ comp_site)
res.aov.coq.WA

# Test post-hoc
pwc.ch.WA  <- data.Met.ch %>% tukey_hsd(totWCa ~ comp_site)
pwc.ch.WA 
pwc.coq.WA  <- data.Met.coq[-c(ind),] %>% tukey_hsd(totWCa ~ comp_site)
pwc.coq.WA 

# plot
pwc.ch.WA  <- pwc.ch.WA  %>% add_xy_position(x = "comp_site")
ggplot.ch.W.Ca  <- ggboxplot(data.Met.ch, x = "comp_site", y = "totWCa", fill="comp_site") +
  stat_pvalue_manual(pwc.ch.WA, hide.ns = TRUE) +
  labs(subtitle = get_test_label(res.aov.ch.WA, detailed = TRUE),caption = get_pwc_label(pwc.ch.WA)) +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold")) +
  ylab("Concentration [ppm]") + ggtitle("Métaux totaux - chair (sans Ca)") + xlab("Campagne")

pwc.coq.WA <- pwc.coq.WA %>% add_xy_position(x = "comp_site")
ggplot.coq.W.Ca <- ggboxplot(data.Met.coq[-c(ind),], x = "comp_site", y = "totWCa", fill="comp_site") +
  stat_pvalue_manual(pwc.coq.WA, hide.ns = TRUE) +
  labs(subtitle = get_test_label(res.aov.coq.WA, detailed = TRUE),caption = get_pwc_label(pwc.coq.WA)) +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold")) +
  ylab("Concentration [ppm]") + ggtitle("Métaux totaux - coquille (sans Ca)") + xlab("Campagne")

rm(pwc.ch.WA,pwc.coq.WA,res.aov.ch.WA,res.aov.coq.WA,ind)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# on refait mais cette fois-ci qu'avec le plot des metaux cibles
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# calcul des metaux cibles
Met.tox.ch <- as.data.frame(data.Met.ch$Cd+data.Met.ch$Cu+data.Met.ch$Hg+
                            data.Met.ch$Pb+data.Met.ch$Zn)
colnames(Met.tox.ch) <- c("Met.tox.ch")
data.Met.ch <- cbind(data.Met.ch,Met.tox.ch)
rm(Met.tox.ch)

Met.tox.coq <- as.data.frame(data.Met.coq$Cd+data.Met.coq$Cu+data.Met.coq$Hg+
                               data.Met.coq$Pb+data.Met.coq$Zn)
colnames(Met.tox.coq) <- c("Met.tox.coq")
data.Met.coq <- cbind(data.Met.coq,Met.tox.coq)
rm(Met.tox.coq)

#IDENTIFIE LES OUTLIERS EXTREMES ET LES GARDES DANS UN DF _ Chair
as.data.frame(data.Met.ch %>% group_by(comp_site) %>% identify_outliers(Met.tox.ch))
as.data.frame(data.Met.coq %>% group_by(comp_site) %>% identify_outliers(Met.tox.coq))
# pas de outliers extreme pour la chair et la coquille !!!

# l'homogeneite des varainces n'est pas assuree pour la chair

data.Met.ch <- data.Met.ch %>% # ordone les donnees d'une certaine maniere
  reorder_levels(comp_site, order = c("VD_A", "SP_A", "VD_H"))
data.Met.coq <- data.Met.coq %>% # ordone les donnees d'une certaine maniere
  reorder_levels(comp_site, order = c("VD_A", "SP_A", "VD_H"))


res.kruskal1 <- data.Met.ch %>% kruskal_test(Met.tox.ch ~ comp_site)

# TEST ANOVA
res.aov.coq <- data.Met.coq  %>% anova_test(Met.tox.coq ~ comp_site)
res.aov.coq
# Test post-hoc
pwc.coq <- data.Met.coq  %>% tukey_hsd(Met.tox.coq ~ comp_site)
pwc.coq

# plot
pwc.ch <- data.Met.ch %>% dunn_test(Met.tox.ch ~ comp_site, p.adjust.method = "bonferroni") 
pwc.ch <- pwc.ch %>% add_xy_position(x = "comp_site")
ggplot.ch.tox <- ggboxplot(data.Met.ch, x = "comp_site", y = "Met.tox.ch", fill="comp_site") +
  stat_pvalue_manual(pwc.ch, hide.ns = TRUE) +
  labs(subtitle = get_test_label(res.kruskal1, detailed = TRUE),caption = get_pwc_label(pwc.ch)) +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold")) +
  ylab("Concentration [ppm]") + ggtitle("Métaux cibles - chair") + xlab("Campagne")

pwc.coq <- pwc.coq %>% add_xy_position(x = "comp_site")
ggplot.coq.tox <- ggboxplot(data.Met.coq , x = "comp_site", y = "Met.tox.coq", fill="comp_site") +
  stat_pvalue_manual(pwc.coq, hide.ns = TRUE) +
  labs(subtitle = get_test_label(res.aov.coq, detailed = TRUE),caption = get_pwc_label(pwc.coq)) +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold")) +
  ylab("Concentration [ppm]") + ggtitle("Métaux cibles - coquille") + xlab("Campagne")

rm(pwc.ch,pwc.coq,res.aov.coq,res.kruskal1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# on refait mais cette fois-ci qu'avec le plot du Ca
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

out <- as.data.frame(data.Met.ch %>% group_by(comp_site) %>% identify_outliers(Ca))
ext <- subset(out,out$is.extreme=="TRUE")
ind <- matrix(ncol=1,nrow=1) # trouve les index
for( i in 1:1){
  ind[i] <- which(data.Met.ch$pool==ext$pool[i])
}
rm(ext,out,i) # nettoie l'environnement

as.data.frame(data.Met.coq %>% group_by(comp_site) %>% identify_outliers(Ca))

# la normalite n'est pas assure pour au moins une serie de donnees --> Kruskal-Willis

#KRUSKAL-WILLIS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#H₀ : toutes les medianes sont egales
#H₁ : au moins une mediane est differente
#seuil de significtion : 0.05

res.kruskal1 <- data.Met.ch[-c(ind),] %>% kruskal_test(Ca ~ comp_site)
res.kruskal2 <- data.Met.coq %>% kruskal_test(Ca ~ comp_site)

#le test de Kruskal-Willis nous indique q'au moins un groupe est significativement
# differente des autres. Procedons a un test de Dunn 

#DUNN~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data.Met.ch <- data.Met.ch %>% # ordone les donnees d'une certaine maniere
  reorder_levels(comp_site, order = c("VD_A", "SP_A", "VD_H"))

pwc1 <- data.Met.ch[-c(ind),] %>% dunn_test(Ca ~ campagne, p.adjust.method = "bonferroni") 
pwc1 <- pwc1 %>% add_xy_position(x = "campagne")
ggplot.ch.Ca <- ggboxplot(data.Met.ch[-c(ind),], x = "campagne", y = "Ca", fill="comp_site") +
  stat_pvalue_manual(pwc1, hide.ns = TRUE)  +
  labs(subtitle = get_test_label(res.kruskal1, detailed = TRUE),caption = get_pwc_label(pwc1)) +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold")) +
  ylab("Concentration [ppm]") + ggtitle("Ca total - chair") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

pwc2 <- data.Met.coq %>% dunn_test(Ca ~ campagne, p.adjust.method = "bonferroni") 
pwc2 <- pwc2 %>% add_xy_position(x = "campagne")
ggplot.coq.Ca <- ggboxplot(data.Met.coq, x = "campagne", y = "Ca", fill="comp_site") +
  stat_pvalue_manual(pwc2, hide.ns = TRUE)  +
  labs(subtitle = get_test_label(res.kruskal2, detailed = TRUE),caption = get_pwc_label(pwc2)) +
  theme(legend.position="none", plot.title = element_text(size=20,hjust = 0.5),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold")) +
  ylab("Concentration [ppm]") + ggtitle("Ca total - coquille") + xlab("Campagne") +
  scale_x_discrete(labels=c("1" = "VD_A", "2" = "SP_A","3" = "VD_H"))

rm(ind,pwc1,pwc2,res.kruskal1,res.kruskal2)

# on plot tous les resultats (annexes)

tit <- ggdraw() + draw_label("Analyse de variances des métaux", fontface='bold', size = 25)
tot <- plot_grid(ggplot.ch.all,ggplot.coq.all,ggplot.ch.W.Ca,ggplot.coq.W.Ca,
                 ggplot.ch.Ca,ggplot.coq.Ca,ggplot.ch.tox,ggplot.coq.tox,nrow = 4, ncol = 2, labels = "AUTO", 
                 label_size = 12, align = "v")
plot_grid(tit,tot, ncol=1, rel_heights=c(0.1, 1))

# on plot tous les resultats (resultat et discussion)

tot <- plot_grid(ggplot.ch.W.Ca,ggplot.coq.W.Ca, 
                 ggplot.ch.Ca,ggplot.coq.Ca,ggplot.ch.tox,ggplot.coq.tox,nrow = 3, ncol = 2, labels = "AUTO", 
                 label_size = 12, align = "v")
plot_grid(tit,tot, ncol=1, rel_heights=c(0.1, 1))

rm(ggplot.ch.all,ggplot.ch.Ca,ggplot.ch.tox,ggplot.ch.W.Ca,
   ggplot.coq.all,ggplot.coq.Ca,ggplot.coq.tox,ggplot.coq.W.Ca,tit,tot)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~ fonctions diverses ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor*r)
}

"coldiss" <- function(D,
                      nc = 4,
                      byrank = TRUE,
                      diag = FALSE) {
  require(gclus)
  
  D <- as.dist(as.matrix(D))
  
  if (max(D) > 1)
    D <- D / max(D)
  
  if (byrank) {
    spe.color <- dmat.color(1 - D, cm.colors(nc))
  }
  else {
    spe.color <- dmat.color(1 - D, byrank = FALSE, cm.colors(nc))
  }
  
  spe.o <- order.single(1 - D)
  speo.color <- spe.color[spe.o, spe.o]
  
  op <- par(mfrow = c(1, 2), pty = "s")
  
  if (diag) {
    plotcolors(
      spe.color,
      rlabels = attributes(D)$Labels,
      main = "Dissimilarity Matrix",
      dlabels = attributes(D)$Labels
    )
    plotcolors(
      speo.color,
      rlabels = attributes(D)$Labels[spe.o],
      main = "Ordered Dissimilarity Matrix",
      dlabels = attributes(D)$Labels[spe.o]
    )
  }
  else {
    plotcolors(spe.color, rlabels = attributes(D)$Labels,
               main = "Dissimilarity Matrix")
    plotcolors(speo.color,
               rlabels = attributes(D)$Labels[spe.o],
               main = "Ordered Dissimilarity Matrix")
  }
  
  par(op)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~ analyse de correlation IC-IC ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pairs(
  data.morpho[which(data.morpho$campagne==1),6:12],
  lower.panel = panel.smooth,
  upper.panel = panel.cor,
  no.col = TRUE,
  diag.panel = panel.hist,
  method = "kendall",
  main = "Corrélations des indices de condition - Vidy Automne (Kendall)", cex.cor=2
)

pairs(
  data.morpho[which(data.morpho$campagne==2),6:12],
  lower.panel = panel.smooth,
  upper.panel = panel.cor,
  no.col = TRUE,
  diag.panel = panel.hist,
  method = "kendall",
  main = "Corrélations des indices de condition - St-Prex Automne (Kendall)", cex.cor=2
)

pairs(
  data.morpho[which(data.morpho$campagne==3),6:12],
  lower.panel = panel.smooth,
  upper.panel = panel.cor,
  no.col = TRUE,
  diag.panel = panel.hist,
  method = "kendall",
  main = "Corrélations des indices de condition - Vidy Hiver (Kendall)", cex.cor=2
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~ analyse de correlation MTs-IC ~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# calculs et sauvegarde des valeurs de tau
tau <- as.data.frame(matrix(ncol = 7, nrow = 3))
p.value <- as.data.frame(matrix(ncol = 7, nrow = 3))

for (j in 1:3){
  for( i in 6:12){
  
    tau[j,i-5] <- cor(data.MTs[which(data.MTs$campagne==j),i],
                          data.MTs[which(data.MTs$campagne==j),13], method = "kendall")
    
    resum <- cor.test(data.MTs[which(data.MTs$campagne==j),i],
                      data.MTs[which(data.MTs$campagne==j),13], method = "kendall")
    p.value[j,i-5] <- resum$p.value
  
  }
}

colnames(tau) <- c("Longeur [cm]", "Largeur [cm]","Épaisseur [cm]","Poids total [mg]",
                   "Poids total sans eau [mg]","Poids tissu [mg]","Poids coquille [mg]")
rownames(tau) <- c("VD_A","SP_A","VD_H")
colnames(p.value) <- c("Longeur [cm]", "Largeur [cm]","Épaisseur [cm]","Poids total [mg]",
                   "Poids total sans eau [mg]","Poids tissu [mg]","Poids coquille [mg]")
rownames(p.value) <- c("VD_A","SP_A","VD_H")
kendall.MTs <- rbind(tau,p.value)
write.csv(kendall.MTs, file = "Kendall_MTs.csv")

rm(i,j,resum,kendall.MTs)

for(k in 1:3){
  par(mfrow=c(2,4))
  for(u in 6:12){
    model <- lm(data.MTs[which(data.MTs$campagne==k),13] ~
                  data.MTs[which(data.MTs$campagne==k),u], data = data.MTs)
    plot(data.MTs[which(data.MTs$campagne==k),u],
         data.MTs[which(data.MTs$campagne==k),13],
         pch=19, ylab = "Métallothionéines [μg/g]", xlab = colnames(data.MTs[u]))
    mtext(paste("tau-value = ",round(tau[k,u-5],3), "; p-value = ",round(p.value[k,u-5],3)))
    abline(model, col = "red",lwd=3)
    if(u==8 && k == 1){mtext(side=3, line=2, at=-0.07, adj=0.35, cex=1.2, 
                             "Correlation : indices de condition - la métallothionéine (VD_A)")}
    if(u==8 && k == 2){mtext(side=3, line=2, at=-0.07, adj=-0.35, cex=1.2, 
                             "Correlation : indices de condition - la métallothionéine (SP_A)")}
    if(u==8 && k == 3){mtext(side=3, line=2, at=-0.07, adj=0.35, cex=1.2, 
                             "Correlation : indices de condition - la métallothionéine (VD_H)")}
  }
}

rm(model,tau,u,k,p.value)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~ analyse de correlation MTs-met ~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# tentative de correlation entre metaux et metallothioneines
test_corr <- cbind(data.Met.ch[,c(2,14:43,45:46)],data.MTs[c(1:10,20:28,40:54),13])

# calculs et sauvegarde des valeurs de tau
tau <- as.data.frame(matrix(ncol = 32, nrow = 3))
p.value <- as.data.frame(matrix(ncol = 32, nrow = 3))

for (j in 1:3){
  for( i in 2:33){
    
    tau[j,i-1] <- cor(test_corr[which(test_corr$campagne==j),i],
                      test_corr[which(test_corr$campagne==j),34], method = "kendall")
    
    resum <- cor.test(test_corr[which(test_corr$campagne==j),i],
                      test_corr[which(test_corr$campagne==j),34], method = "kendall")
    p.value[j,i-1] <- resum$p.value
    
  }
}

rm(tau,p.value,test_corr,i,j,resum)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~ analyse de correlation IC-met ~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculs et sauvegarde des valeurs de tau et des p.value chair
tau.ch <- as.data.frame(matrix(ncol = 1, nrow = 3*4))
p.value <- as.data.frame(matrix(ncol = 1, nrow = 3*4))
c <- 1
l <- 43 # met tot

for (j in 1:12){
    
    tau.ch[j,] <- cor(data.Met.ch[which(data.Met.ch$campagne==c),6],
                         data.Met.ch[which(data.Met.ch$campagne==c),l], method = "kendall")
    p_val <- cor.test(data.Met.ch[which(data.Met.ch$campagne==c),6],
                      data.Met.ch[which(data.Met.ch$campagne==c),l], method = "kendall")
    p.value[j,] <- p_val$p.value
  
  c <- c + 1
  if(c==4){c<-1}
  if(j==3){l<-45} # met total sans ca
  if(j==6){l<-19} # met  Ca
  if(j==9){l<-46} # met cib
  
}

colnames(tau.ch) <- c("Poids tissus total [mg]")
rownames(tau.ch) <- c("tot_1","tot_2","tot_3","Wca_1","Wca_2","Wca_3",
                      "Ca_1","Ca_2","Ca_3","cib_1","cib_2","cib_3")
colnames(p.value) <- c("Poids tissus total [mg]")
rownames(p.value) <- c("tot_1","tot_2","tot_3","Wca_1","Wca_2","Wca_3",
                       "Ca_1","Ca_2","Ca_3","cib_1","cib_2","cib_3")

rm(j,c,l,p_val)

mu <- 1
met <- c("Métaux totaux", "Métaux totaux (sans Ca)", "Calcium", "Métaux cible")

for(j in 1:4){
  
  par(mfrow=c(1,3))
  l <- 43
  
  for(k in 1:3){
    if(j==2){l<-45} # met total sans Ca
    if(j==3){l<-19} # Ca
    if(j==4){l<-46} # met cib
      model <- lm(data.Met.ch[which(data.Met.ch$campagne==k),l] ~
                  data.Met.ch[which(data.Met.ch$campagne==k),6], data = data.Met.ch)
      plot(data.Met.ch[which(data.Met.ch$campagne==k),6],
           data.Met.ch[which(data.Met.ch$campagne==k),l],
           pch=19, ylab = paste("Concentrations",met[j],"[ppm]"), 
           xlab = colnames(data.Met.ch[6]))
      mtext(paste("tau-value = ",round(tau.ch[mu,],3), "; p-value = ", round(p.value[mu,],3)))
      abline(model, col = "red",lwd=3)
      
      if(k == 1){title(main="IC - métaux chair (VD_A)", cex.main=1.2)}
      
      if(k == 2){title(main="IC - métaux chair (SP_A)", cex.main=1.2)}
      
      if(k == 3){title(main="IC - métaux chair (VD_H)", cex.main=1.2)}
      
      mu <- mu + 1
      
  }
}

rm(model,j,k,l,mu)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculs et sauvegarde des valeurs de tau et p.value - coquille
tau.coq <- as.data.frame(matrix(ncol = 1, nrow = 3*4))
p.value.coq <- as.data.frame(matrix(ncol = 1, nrow = 3*4))
c <- 1
l <- 43 # met tot

for (j in 1:12){
  
  tau.coq[j,] <- cor(data.Met.coq[which(data.Met.coq$campagne==c),7],
                    data.Met.coq[which(data.Met.coq$campagne==c),l], method = "kendall")
  p_val <- cor.test(data.Met.coq[which(data.Met.coq$campagne==c),7],
                    data.Met.coq[which(data.Met.coq$campagne==c),l], method = "kendall")
  p.value.coq[j,] <- p_val$p.value
  
  c <- c + 1
  if(c==4){c<-1}
  if(j==3){l<-45} # met total sans ca
  if(j==6){l<-19} # met  Ca
  if(j==9){l<-46} # met cib
  
}

colnames(tau.coq) <- c("Poids coquille total [mg]")
rownames(tau.coq) <- c("tot_1","tot_2","tot_3","Wca_1","Wca_2","Wca_3",
                       "Ca_1","Ca_2","Ca_3","cib_1","cib_2","cib_3")
colnames(p.value.coq) <- c("Poids coquille total [mg]")
rownames(p.value.coq) <- c("tot_1","tot_2","tot_3","Wca_1","Wca_2","Wca_3",
                       "Ca_1","Ca_2","Ca_3","cib_1","cib_2","cib_3")

rm(j,c,l,p_val)

mu <- 1

for(j in 1:4){
  
  par(mfrow=c(1,3))
  l <- 43
  
  for(k in 1:3){
    if(j==2){l<-45} # met total sans Ca
    if(j==3){l<-19} # Ca
    if(j==4){l<-46} # met cib
    model <- lm(data.Met.coq[which(data.Met.coq$campagne==k),l] ~
                  data.Met.coq[which(data.Met.coq$campagne==k),7], data = data.Met.coq)
    plot(data.Met.coq[which(data.Met.coq$campagne==k),7],
         data.Met.coq[which(data.Met.coq$campagne==k),l],
         pch=19, ylab =  paste("Concentrations",met[j],"[ppm]"), 
         xlab = colnames(data.Met.coq[7]))
    mtext(paste("tau-value = ",round(tau.coq[mu,],3), "; p-value = ", round(p.value.coq[mu,],3)))
    abline(model, col = "red",lwd=3)
    
    if(k == 1){title(main="IC - métaux coquille (VD_A)", cex.main=1.2)}
    
    if(k == 2){title(main="IC - métaux coquille  (SP_A)", cex.main=1.2)}
    
    if(k == 3){title(main="IC - métaux coquille  (VD_H)", cex.main=1.2)}
    
    mu <- mu + 1
    
  }
}

# enregistrement du tableau de données

col.tau <- cbind(tau.ch,tau.coq)
col.pvalue <- cbind(p.value,p.value.coq)
kendall <- rbind(col.tau,col.pvalue)

write.csv(kendall,file = "kendallCorrMet.csv")

rm(model,met,j,k,l,mu,tau.coq,coldiss,panel.cor,panel.hist,
   p.value,col.tau,col.pvalue,tau.ch,p.value.coq,kendall)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~ statistiques multivariees ~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Analyse discriminante linéaire ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
# lda met tox ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# creation tableau avec que les metaux cibles
Met.tox.ch <- cbind(data.Met.ch[15:17],data.Met.ch[20],data.Met.ch[22:25],data.Met.ch[28],data.Met.ch[31],
                    data.Met.ch[33],data.Met.ch[35],data.Met.ch[42])

Met.tox.coq <- cbind(data.Met.coq[15:17],data.Met.coq[20],data.Met.coq[22:25],data.Met.coq[28],
                     data.Met.coq[31],data.Met.coq[33],data.Met.coq[35],data.Met.coq[42])

Met.tox.ch <- cbind(paste(data.Met.ch$comp_site,"_ch"),Met.tox.ch)
Met.tox.coq <- cbind(paste(data.Met.coq$comp_site,"_coq"),Met.tox.coq)

rownames(Met.tox.ch) <- data.Met.ch[,1]
rownames(Met.tox.coq) <- data.Met.ch[,1]
colnames(Met.tox.ch) <- c("campagne","Al","As","Ba","Cd","Cr","Cu","Fe","Hg","Mn","Ni","Pb","Sb","Zn")
colnames(Met.tox.coq) <- c("campagne","Al","As","Ba","Cd","Cr","Cu","Fe","Hg","Mn","Ni","Pb","Sb","Zn")
Met.tox <- rbind(Met.tox.ch,Met.tox.coq)
rm(Met.tox.ch,Met.tox.coq)

# 70% pour le training et 30 pour le test
sample <- sample(c(TRUE, FALSE), nrow(Met.tox), replace=TRUE, prob=c(0.6,0.4))
train1 <- Met.tox[sample, ]
test <- Met.tox[!sample, ] 

#fit le modele LDA
model <- lda(campagne~., data=train1)
#Output
model

# utilise la model pour la prediction sur les donnees test 
predicted <- predict(model, test)
names(predicted)

# la precision du model
prec <- round(mean(predicted$class==test$campagne),2)

# defini les donnees a ploter
lda_plot <- cbind(train1, predict(model)$x)

#creer le plot
lda.tox <- ggplot(lda_plot, aes(LD1, LD2)) +
           geom_point(aes(color = train1$campagne)) +
           theme(plot.title = element_text(size=15,hjust = 0.5),
                 plot.subtitle = element_text(hjust = 0.5),
                 plot.caption = element_text(size=12),
                 legend.position = "none") +
           labs(title = "ADL - métaux cibles", 
                subtitle = paste("Précision du modèle : ",prec),
                caption = "training : 60%, test 40 %",
                color = "Campagne") +
           geom_hline(yintercept=0) + geom_vline(xintercept=0)

rm(model,predicted,prec,test,sample,lda_plot)
              

# lda met tot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
# creation tableau avec tous les metaux

ch <- cbind(paste(data.Met.ch$comp_site,"_ch"),data.Met.ch)
coq <- cbind(paste(data.Met.coq$comp_site,"_coq"),data.Met.coq)
names.ch <- colnames(ch)
names.ch[1] <- c("campagne")
colnames(ch) <- names.ch
names.coq <- colnames(coq)
names.coq[1] <- c("campagne")
colnames(coq) <- names.coq
Met.tot <- rbind(ch[,c(1:13,15:45)],coq[,c(1:13,15:45)])
rm(ch,coq,names.ch,names.coq)

Met.tot.CA <- Met.tot[,c(1,14:42)]

# 70% pour le training et 30 pour le test
sample <- sample(c(TRUE, FALSE), nrow(Met.tot.CA), replace=TRUE, prob=c(0.6,0.4))
train2 <- Met.tot.CA[sample, ]
test <- Met.tot.CA[!sample, ] 

#fit le modele LDA
model <- lda(campagne~., data=train2)
#Output
model

# utilise la model pour la prediction sur les donnees test 
predicted <- predict(model, test)
names(predicted)

# la precision du model
prec <- round(mean(predicted$class==test$campagne),2)

# defini les donnees a ploter
lda_plot <- cbind(train2, predict(model)$x)

#creer le plot
lda.tot <-  ggplot(lda_plot, aes(LD1, LD2)) +
            geom_point(aes(color = train2$campagne)) +
            theme(plot.title = element_text(size=15,hjust = 0.5),
                  plot.subtitle = element_text(hjust = 0.5),
                  plot.caption = element_text(size=12),
                  legend.position = "none") +
            labs(title = "ADL - métaux totaux", 
                 subtitle = paste("Précision du modèle : ",prec),
                 caption = "training : 60%, test 40 %", 
                 color = "Campagne") +
            geom_hline(yintercept=0) + geom_vline(xintercept=0)

rm(model,predicted,test,prec,sample,lda_plot)

# lda met sans Ca ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# creation tableau avec tous les metaux sauf Ca

Met.tot.WCA <- Met.tot[,c(1,14:18,20:42)]

# 70% pour le training et 30 pour le test
sample <- sample(c(TRUE, FALSE), nrow(Met.tot.WCA), replace=TRUE, prob=c(0.6,0.4))
train3 <- Met.tot.WCA[sample, ]
test <- Met.tot.WCA[!sample, ] 

#fit le modele LDA
model <- lda(campagne~., data=train3)
#Output
model

# utilise la model pour la prediction sur les donnees test 
predicted <- predict(model, test)
names(predicted)

# la precision du model
prec <- round(mean(predicted$class==test$campagne),2)

# defini les donnees a ploter
lda_plot <- cbind(train3, predict(model)$x)

#creer le plot
lda.totWCa <-  ggplot(lda_plot, aes(LD1, LD2)) +
  geom_point(aes(color = train3$campagne)) +
  theme(plot.title = element_text(size=15,hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(size=12),
        legend.position = c(0.83, 0.12)) +
  labs(title = "ADL - métaux totaux sans Ca", 
       subtitle = paste("Précision du modèle : ",prec),
       caption = "training : 60%, test 40 %", 
       color = "Campagne") +
  geom_hline(yintercept=0) + geom_vline(xintercept=0)

rm(model,predicted,test,prec,sample,lda_plot)


plot_grid(lda.tot,lda.totWCa,lda.tox, nrow = 1, ncol = 3,labels = "AUTO", 
          label_size = 12, align = "v")

rm(train1,train2,train3,lda.tot,lda.tox,lda.totWCa)

# lda MTs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# creation tableau 
MTs.lda <- data.MTs[,c(3,6:13)]

# 70% pour le training et 30 pour le test
sample <- sample(c(TRUE, FALSE), nrow(MTs.lda), replace=TRUE, prob=c(0.7,0.3))
train <- MTs.lda[sample, ]
test <- MTs.lda[!sample, ] 

#fit le modele LDA
model <- lda(comp_site~., data=train)
#Output
model

# utilise la model pour la prediction sur les donnees test 
predicted <- predict(model, test)
names(predicted)

# la precision du model
prec <- round(mean(predicted$class==test$comp_site),2)

# defini les donnees a ploter
lda_plot <- cbind(train, predict(model)$x)

#creer le plot
ggplot(lda_plot, aes(LD1, LD2)) +
  geom_point(aes(color = train$comp_site)) +
  theme(plot.title = element_text(size=15,hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = c(0.83, 0.12)) +
  labs(title = "ADL - métaux totaux et métallothionéines", 
       subtitle = paste("Précision du modèle : ",prec),
       caption = "training : 60%, test 40 %", 
       color = "Campagne") +
  geom_hline(yintercept=0) + geom_vline(xintercept=0)

rm(model,predicted,test,prec,sample,lda_plot,MTs.lda,train)

# PCA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# analyse PCA
met.tot.CA.pca <- PCA(Met.tot.CA[,2:30], scale.unit = TRUE, ncp = 5, graph = TRUE)
met.tot.WCA.pca <- PCA(Met.tot.WCA[,c(2:29)], scale.unit = TRUE, ncp = 5, graph = TRUE)
met.tox.pca <- PCA(Met.tox[,c(5,7,9,12,14)], scale.unit = TRUE, ncp = 5, graph = TRUE)

# Contributions of variables to PC1
cont1 <- fviz_contrib(met.tot.CA.pca, choice = "var", axes = 1, title = "Contribution dim 1 - métaux totaux")
cont2 <- fviz_contrib(met.tox.pca, choice = "var", top = 15, axes = 1, title = "Contribution dim 1 - métaux cibles totaux")
# Contributions of variables to PC2
cont3 <- fviz_contrib(met.tot.CA.pca, choice = "var", top = 15, axes = 2,title = "Contribution dim 2 - métaux totaux")
cont4 <- fviz_contrib(met.tox.pca, choice = "var", top = 15, axes = 2, title = "Contribution dim 2 - métaux cibles totaux")

plot_grid(cont1,cont2,cont3,cont4,ncol = 2, nrow = 2,labels = "AUTO", 
          label_size = 12, align = "v")

rm(cont1,cont2,cont3,cont4)

# contribution de chaque variables
var.CA <- get_pca_var(met.tot.CA.pca)
var.WCA <- get_pca_var(met.tot.WCA.pca)
var.tox <- get_pca_var(met.tox.pca)

pca_1 <- fviz_pca_var(met.tot.CA.pca, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE, 
             title = "ACP (contribution) - métaux totaux"
)
pca_2 <- fviz_pca_var(met.tot.WCA.pca, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE, 
             title = "ACP (contribution) - métaux totaux sans Ca"
)
pca_3 <- fviz_pca_var(met.tox.pca, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE, 
             title = "ACP (contribution) - métaux cibles totaux"
)

# plot nnexes 
plot_grid(pca_1,pca_2,pca_3,ncol = 2, nrow = 2,labels = "AUTO", 
          label_size = 12, align = "v")
#plot rapport
plot_grid(pca_1,ncol=1,nrow=1,label_size = 12)


# visualisation des groupes
pca_1 <- fviz_pca_biplot(met.tot.CA.pca,
                         geom.ind = "point",
                         col.ind = Met.tot.CA$campagne, 
                         addEllipses = TRUE, 
                         legend.title = "Groupes",
                         title = "ACP groupes - Métaux totaux"
)
pca_2 <- fviz_pca_biplot(met.tot.WCA.pca,
                         geom.ind = "point", 
                         col.ind = Met.tot.WCA$campagne, 
                         addEllipses = TRUE, 
                         legend.title = "Groupes",
                         title = "ACP groupes - Métaux totaux sans Ca"
)
pca_3 <- fviz_pca_biplot(met.tox.pca,
                         geom.ind = "point", 
                         col.ind = Met.tox$campagne, 
                         addEllipses = TRUE, 
                         legend.title = "Groupes",
                         title = "ACP groupes - Métaux cibles totaux"
)

# plot annexes
plot_grid(pca_1,pca_2,pca_3,ncol = 2, nrow = 2,labels = "AUTO", 
          label_size = 12, align = "v")
# plot rapport
plot_grid(pca_1,ncol = 1, nrow = 1, 
          label_size = 12, align = "v")

# analyse PCA 2eme fois en supprimant des parametre environnement redondant
met.tot.CA.pca <- PCA(Met.tot.CA[,c(6:8,13:15,19,27,30)], scale.unit = TRUE, ncp = 5, graph = TRUE)
met.tot.WCA.pca <- PCA(Met.tot.WCA[,c(6,8,13:15,19,27,29)], scale.unit = TRUE, ncp = 5, graph = TRUE)
met.tox.pca <- PCA(Met.tox[,c(9,11,13,14)], scale.unit = TRUE, ncp = 5, graph = TRUE)

# visualisation des groupes
pca_1 <- fviz_pca_biplot(met.tot.CA.pca,
             geom.ind = "point",
             col.ind = Met.tot.CA$campagne, 
             addEllipses = TRUE, 
             legend.title = "Groupes",
             title = "ACP groupes - Métaux totaux"
)
pca_2 <- fviz_pca_biplot(met.tot.WCA.pca,
                      geom.ind = "point", 
                      col.ind = Met.tot.WCA$campagne, 
                      addEllipses = TRUE, 
                      legend.title = "Groupes",
                      title = "ACP groupes - Métaux totaux sans Ca"
)
pca_3 <- fviz_pca_biplot(met.tox.pca,
                      geom.ind = "point", 
                      col.ind = Met.tox$campagne, 
                      addEllipses = TRUE, 
                      legend.title = "Groupes",
                      title = "ACP groupes - Métaux cibles totaux"
)

#plotannexes
plot_grid(pca_1,pca_2,pca_3,ncol = 2, nrow = 2,labels = "AUTO", 
          label_size = 12, align = "v")
#plot rapport
plot_grid(pca_1,ncol = 1, nrow = 1,
          label_size = 12, align = "v")

rm(met.tot.CA.pca,met.tot.WCA.pca,met.tox.pca,pca_1,pca_2,pca_3,var.CA,
   var.tox,var.WCA)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~ analyses complémentaires ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Comparaison des métaux dans les sediments, la chair et la coquille Vidy ~~~~~~

sedi.vd <- apply(sedi[18,2:10],2,mean)
met.ch.vd.a <- apply(Met.tot[1:10,c(20,22,23,24,25,31,42,33,28)],2,mean)
#met.ch.vd.h <- apply(Met.tot[20:34,c(20,22,23,24,25,31,42,33,28)],2,mean)
met.coq.vd.a <- apply(Met.tot[35:44,c(20,22,23,24,25,31,42,33,28)],2,mean)
#met.coq.vd.h <- apply(Met.tot[54:68,c(20,22,23,24,25,31,42,33,28)],2,mean)

sd.sedi <- apply(sedi[18,2:10],2,sd)
sd.met.ch <- apply(Met.tot[1:10,c(20,22,23,24,25,31,42,33,28)],2,sd)
sd.met.coq <- apply(Met.tot[35:44,c(20,22,23,24,25,31,42,33,28)],2,sd)

sedi.moules <- rbind(sedi.vd,met.ch.vd.a,met.coq.vd.a)
sedi.moules.sd <- rbind(sd.sedi,sd.met.ch,sd.met.coq)
  
rm(sedi.vd,met.ch.vd.a,met.coq.vd.a,sd.sedi,sd.met.ch,sd.met.coq)

base_r_barplot <- barplot(abs(log(sedi.moules)),
                          main = "Comparaison des concentrations des métaux - Vidy",
                          xlab = "Métal",
                          ylab = ("log(moyenne de la concentration)"),
                          col = wes_palette(n=3, name="Moonrise2"),
                          beside = TRUE,
                          ylim=c(0,15))

arrows(x0 = base_r_barplot,                           # Add error bars
       y0 = abs(log(sedi.moules)) + abs(log(sedi.moules.sd)),
       y1 = abs(log(sedi.moules)) - 0,
       angle = 90,
       code = 3,
       length = 0.05)

legend(20,15,
       c("Sédiments","Chair","Coquille"),
       fill = wes_palette(n=3, name="Moonrise2")
)

rm(sedi.moules,sedi.moules.sd,base_r_barplot)

# Comparaison des métaux dans les sediments, la chair et la coquille St-Prex ~~~

sedi.sp <- apply(sedi[19,2:10],2,mean)
met.ch.sp.a <- apply(Met.tot[11:19,c(20,22,23,24,25,31,42,33,28)],2,mean)
#met.ch.vd.h <- apply(Met.tot[20:34,c(20,22,23,24,25,31,42,33,28)],2,mean)
met.coq.sp.a <- apply(Met.tot[45:64,c(20,22,23,24,25,31,42,33,28)],2,mean)
#met.coq.vd.h <- apply(Met.tot[54:68,c(20,22,23,24,25,31,42,33,28)],2,mean)


sd.sedi <- apply(sedi[19,2:10],2,sd)
sd.met.ch <- apply(Met.tot[11:19,c(20,22,23,24,25,31,42,33,28)],2,sd)
sd.met.coq <- apply(Met.tot[45:64,c(20,22,23,24,25,31,42,33,28)],2,sd)

sedi.moules <- rbind(sedi.sp,met.ch.sp.a,met.coq.sp.a)
sedi.moules.sd <- rbind(sd.sedi,sd.met.ch,sd.met.coq)

rm(sedi.sp,met.ch.sp.a,met.coq.sp.a,sd.sedi,sd.met.ch,sd.met.coq)

base_r_barplot <- barplot(abs(log(sedi.moules)),
                          main = "Comparaison des concentrations des métaux - St-Prex",
                          xlab = "Métal",
                          ylab = ("log(moyenne de la concentration)"),
                          col = wes_palette(n=3, name="Moonrise2"),
                          beside = TRUE,
                          ylim=c(0,15))

arrows(x0 = base_r_barplot,                           # Add error bars
       y0 = abs(log(sedi.moules)) + abs(log(sedi.moules.sd)),
       y1 = abs(log(sedi.moules)) - 0,
       angle = 90,
       code = 3,
       length = 0.05)

legend(20,15,
       c("Sédiments","Chair","Coquille"),
       fill = wes_palette(n=3, name="Moonrise2")
)

rm(sedi.moules,sedi.moules.sd,base_r_barplot)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~ tentatives supplementaires ~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# RDA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

inter1 <- data.MTs
rownames(inter1) <- c(1:length(inter1[,1]))
mts <- as.data.frame(inter1[c(1:10,20:28,37:51),13])
colnames(mts) <- c("Mts [ug/g]")
inter2 <- cbind(Met.tot[1:34,1:43],mts)
Data.RDA <- cbind(inter2)

rm(inter1,inter2,mts)

reponse <- Data.RDA[1:10,c(20,22,24,25,31,42,33,28,44)]
explica <- sedi[9:18,c(2,3,5:10)]

rep.hel <- decostand(reponse, "hellinger")
rep.rda <- rda(rep.hel ~ ., explica)
summary(rep.rda)

(R2adj <- RsquareAdj(rep.rda)$adj.r.squared)

par(mfrow=c(1,2))
plot(rep.rda, disp=c("sp", "lc", "cn"), scaling=1, 
     main = "RDA - scale 1")
plot(rep.rda, disp=c("sp", "lc", "cn"), scaling=2, 
     main = "RDA - scale 2")

anova.cca(rep.rda, step = 1000)

rm(rep.hel, R2adj, rep.rda, reponse, explica)

# regression linéaire

model1 <- lm(Data.RDA$`Mts [ug/g]` ~ Ag + Al + As + Ba +  Be + Ca + Cd +Co 
            + Cr  + Cu + Fe + Hg + K + Mg + Mn + Mo + Na + Ni + P + Pb + S + Sb 
            + Se + Si + Sn + Sr + Ti + V + Zn , data = Data.RDA) 
summary(model1)

model2 <- lm(Data.RDA$`Mts [ug/g]` ~ As +  Be + Ca + Cd + Co 
             + Cr  + Fe + Hg + K + Mg + Mo + Ni + P + S  
             + Se + Si  + Ti + Zn , data = Data.RDA) 
summary(model2)

model3 <- lm(Data.RDA$`Mts [ug/g]` ~ Ni + Zn + As + Cr + Al + Pb + Mn + Fe +
               Hg + Sb + Cd + Fe, data = Data.RDA) 
summary(model3)

rm(model1,model2,model3)
