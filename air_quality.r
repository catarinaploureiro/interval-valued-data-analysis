###
##
# Apply the method proposed by Lin et al. (2022) to an air quality dataset from Entrecampos, Lisbon.
##
###

source("./setup.r")

set.seed(13) #for reproducibility

library("readxl")
library("data.table")
library("ggplot2")
library("imputeTS")
library("zoo")
library("car")
library("geigen")
library("factoextra")
library("extRemes")
library("rrcov")
library("isotree")
library("gplots")

###
# Preparing the dataset
###

#Import data & convert 1st column to date
entrecampos_2019<-as.data.frame(read_excel('./Data/Entrecampos_2019-01-01_2019-12-31.xlsx', 
                    col_types = c("text","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric")))
entrecampos_2019[,1]<-as.POSIXct(entrecampos_2019[,1], format="%Y-%m-%d %H:%M:%S")

entrecampos_2020<-as.data.frame(read_excel('./Data/Entrecampos_2020-01-01_2020-12-31.xlsx', 
                    col_types = c("text","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric")))
entrecampos_2020[,1]<-as.POSIXct(entrecampos_2020[,1], format="%Y-%m-%d %H:%M:%S")

entrecampos_2021<-as.data.frame(read_excel('./Data/Entrecampos_2021-01-01_2021-12-31.xlsx', 
                    col_types = c("text","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric")))
entrecampos_2021[,1]<-as.POSIXct(entrecampos_2021[,1], format="%Y-%m-%d %H:%M:%S")

#Due to change to "summer time", clocks went forward from 01h to 02h, so as.POSIXct returns NA
entrecampos_2019[2138,1] <- as.POSIXct("2019-03-31 00:59:59", format="%Y-%m-%d %H:%M:%S")
entrecampos_2020[2114,1] <- as.POSIXct("2020-03-29 00:59:59", format="%Y-%m-%d %H:%M:%S")
entrecampos_2021[2066,1] <- as.POSIXct("2021-03-28 00:59:59", format="%Y-%m-%d %H:%M:%S")

#Create new column with Day
entrecampos_2019[,11]<-as.Date(entrecampos_2019[,1])
colnames(entrecampos_2019)[11]<-'Day'
entrecampos_2020[,11]<-as.Date(entrecampos_2020[,1])
colnames(entrecampos_2020)[11]<-'Day'
entrecampos_2021[,11]<-as.Date(entrecampos_2021[,1])
colnames(entrecampos_2021)[11]<-'Day'

#Append all 3 years into 1 data frame
entrecampos<-rbind(entrecampos_2019,entrecampos_2020,entrecampos_2021)

#Check missing values
sort(colSums(is.na(entrecampos))/26304) #Benzeno 20% NAs, Ozono 10% NAs
statsNA(entrecampos[,7])

#Remove Benzene
entrecampos<-entrecampos[, -7]

colnames(entrecampos)<-c("Timestamp", "Sulphur Dioxide (µg/m3)", "Particles < 10 µm (µg/m3)",
                    "Ozone (µg/m3)", "Nitrogen Dioxide (µg/m3)", "Carbon Monoxide (µg/m3)",
                    "Particles < 2.5 µm (µg/m3)", "Nitrogen Oxides (µg/m3)", "Nitrogen Monoxide (µg/m3)","Day")


#Plot histograms
par(mfrow=c(2,4))
for (i in 2:9){
    hist(entrecampos[,i], xlab=NULL, main=paste("Histogram of", colnames(entrecampos)[i]))
}

#Data transformation
entrecampos_log10<-entrecampos
entrecampos_log10[,6]<-1000*entrecampos_log10[,6] #reduce to same unit
entrecampos_log10[,2:9]<-log10(entrecampos_log10[,2:9]+1)

#Aggregate the data with min, max and mean
entrecampos_log10_max<-aggregate(entrecampos_log10[,2:9], by=list(entrecampos[,10]), FUN="max")
entrecampos_log10_min<-aggregate(entrecampos_log10[,2:9], by=list(entrecampos[,10]), FUN="min")
entrecampos_log10_mean<-aggregate(entrecampos_log10[,2:9], by=list(entrecampos[,10]), FUN="mean")

#Interpolation
entrecampos_log10_max_ma<-na_ma(entrecampos_log10_max, k=20, weighting = "linear")
entrecampos_log10_min_intp<-na_interpolation(entrecampos_log10_min, option="linear")
entrecampos_log10_mean_ma<-na_ma(entrecampos_log10_mean, k=20, weighting = "linear")

#Plot daily average mean
entrecampos_mean_1col_log10<-melt(as.data.table(entrecampos_log10_mean_ma), 
                                    id.vars=c('Group.1'), variable.name = 'Pollutant')
entrecampos_mean_1col_log10$Pollutant<-factor(entrecampos_mean_1col_log10$Pollutant, 
                                        levels = c("Carbon Monoxide (µg/m3)", "Nitrogen Dioxide (µg/m3)",
                                        "Nitrogen Monoxide (µg/m3)", "Nitrogen Oxides (µg/m3)",
                                        "Ozone (µg/m3)",  "Particles < 10 µm (µg/m3)",
                                        "Particles < 2.5 µm (µg/m3)",  "Sulphur Dioxide (µg/m3)"))

pdf(file="./Figures/Log/Mean_all_lines_log.pdf", width = 12, height = 4)
ggplot(entrecampos_mean_1col_log10, aes(Group.1, value)) + geom_line(aes(colour = Pollutant)) + 
    theme(text = element_text(size=20), axis.title.x = element_blank(), axis.title.y = element_blank())
dev.off()

#Plot symbolic data with interpolated values
varnames<-c("Group.1", "Sulphur Dioxide (µg/m3)",  "Particles < 10 µm (µg/m3)",  "Ozone (µg/m3)", 
            "Nitrogen Dioxide (µg/m3)",  "Carbon Monoxide (µg/m3)", "Particles < 2.5 µm (µg/m3)", 
            "Nitrogen Oxides (µg/m3)", "Nitrogen Monoxide (µg/m3)")
pdf(file="./Figures/Log/intervals_all_with_imputation_log.pdf", width=12, height=5, pointsize = 14)
par(mfrow=c(2,4), mar = c(3, 4, 1, 1))
for (i in 2:9){
    min_dif<-entrecampos_log10_min_intp[entrecampos_log10_min_intp[,1] %in% entrecampos_log10_min[is.na(entrecampos_log10_min[,i]),1], c(1,i)]
    max_dif<-entrecampos_log10_max_ma[entrecampos_log10_max_ma[,1] %in% entrecampos_log10_max[is.na(entrecampos_log10_max[,i]),1], c(1,i)]
    plot(c(entrecampos_log10_min[, 1],entrecampos_log10_max[,1]), 
        c(entrecampos_log10_min[, i],entrecampos_log10_max[, i]), 
        pch=20,xlab="",ylab=varnames[i], col=c("blue"))
    rect(
        xleft = entrecampos_log10_min[,1],
        xright = entrecampos_log10_max[,1],
        ybottom = entrecampos_log10_min[, i],
        ytop = entrecampos_log10_max[, i],
        density = 50,
        col = "black"
    )
    points(c(min_dif[, 1],max_dif[,1]), c(min_dif[, 2],max_dif[, 2]), pch=20, col=c("red"))
    rect(
        xleft = min_dif[,1],
        xright = max_dif[,1],
        ybottom = min_dif[, 2],
        ytop = max_dif[, 2],
        density = 50,
        col = "darkred"
    )
}
dev.off()

###
# Apply the method from Lin et al. (2022)
###

#Obtain symbolic covariance matrix
alpha<-0.375
n<-24
co_matrix_log10<-kendall_corr(entrecampos_log10_min_intp[,-1],
                                entrecampos_log10_max_ma[,-1],n,alpha,"clayton")

#PCA symbolic
eng_sym_log10<-geigen(co_matrix_log10$cov_est, diag(8))

eng_sym_log10$values #eigenvalues
sum(eng_sym_log10$values[6:8])/sum(eng_sym_log10$values) #80% explained variability with 3 PCs

#Compute pc scores
pca1_sym_log10<-(-eng_sym_log10$vectors[,8])
pca1_sym_score_log10<-pca_scores_symbolic(entrecampos_log10_min_intp[,-1],
                                        entrecampos_log10_max_ma[,-1],pca1_sym_log10)
pca2_sym_log10<-(-eng_sym_log10$vectors[,7])
pca2_sym_score_log10<-pca_scores_symbolic(entrecampos_log10_min_intp[,-1],
                                        entrecampos_log10_max_ma[,-1],pca2_sym_log10)
pca3_sym_log10<-(-eng_sym_log10$vectors[,6])
pca3_sym_score_log10<-pca_scores_symbolic(entrecampos_log10_min_intp[,-1],
                                        entrecampos_log10_max_ma[,-1],pca3_sym_log10)

#Fit a GEV distribution to the pc scores
Max_PC1_log10<-as.numeric(pca1_sym_score_log10$MaximumPC)
Min_PC1_log10<-as.numeric(pca1_sym_score_log10$MinimumPC)
Max_extr_PC1_log10<-fevd(Max_PC1_log10, type="GEV", method = "GMLE")
Min_extr_PC1_log10<-fevd(Min_PC1_log10, type="GEV", method = "MLE")

Max_PC2_log10<-as.numeric(pca2_sym_score_log10$MaximumPC)
Min_PC2_log10<-as.numeric(pca2_sym_score_log10$MinimumPC)
Max_extr_PC2_log10<-fevd(Max_PC2_log10, type="GEV", method = "Lmoments")
Min_extr_PC2_log10<-fevd(Min_PC2_log10, type="GEV", method = "MLE")

Max_PC3_log10<-as.numeric(pca3_sym_score_log10$MaximumPC)
Min_PC3_log10<-as.numeric(pca3_sym_score_log10$MinimumPC)
Max_extr_PC3_log10<-fevd(Max_PC3_log10, type="GEV", method = "Lmoments")
Min_extr_PC3_log10<-fevd(Min_PC3_log10, type="GEV", method = "MLE")

#Plot PC scores histogram and fitted GEV distribution
pdf(file="./Figures/Log/extr_hist_pc_log.pdf", width=12, height=8)
par(mfrow=c(2,3), mar=c(2,2.5,2,1))
plot(Min_extr_PC1_log10, type="hist", col="lightskyblue", main="Minimum PC1s", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, cex.sub=1.8)
plot(Min_extr_PC2_log10, type="hist", ylim=c(0,2), col="lightskyblue", main="Minimum PC2s", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, cex.sub=1.8)
plot(Min_extr_PC3_log10, type="hist", col="lightskyblue", main="Minimum PC3s", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, cex.sub=1.8)
plot(Max_extr_PC1_log10, type="hist", ylim=c(0,0.8), col="lightskyblue", main="Maximum PC1s", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, cex.sub=1.8)
plot(Max_extr_PC2_log10, type="hist", xlim=c(1,3), col="lightskyblue", main="Maximum PC2s", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, cex.sub=1.8)
plot(Max_extr_PC3_log10, type="hist", ylim=c(0,4), col="lightskyblue", main="Maximum PC3s", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, cex.sub=1.8)
dev.off()

#Obtain quantiles
p_max1_log10 <- Max_extr_PC1_log10$results$par
qtl99_max1_log10<-extRemes::qevd(0.99, loc = p_max1_log10[1], 
                                scale = p_max1_log10[2], shape = p_max1_log10[3])
p_min1_log10 <- Min_extr_PC1_log10$results$par
qtl99_min1_log10<-extRemes::qevd(0.01, loc = p_min1_log10[1], 
                                scale = p_min1_log10[2], shape = p_min1_log10[3])

p_max2_log10 <- Max_extr_PC2_log10$results
qtl99_max2_log10<-extRemes::qevd(0.99, loc = p_max2_log10[1], 
                                scale = p_max2_log10[2], shape = p_max2_log10[3])
p_min2_log10 <- Min_extr_PC2_log10$results$par
qtl99_min2_log10<-extRemes::qevd(0.01, loc = p_min2_log10[1], 
                                scale = p_min2_log10[2], shape = p_min2_log10[3])

p_max3_log10 <- Max_extr_PC3_log10$results
qtl99_max3_log10<-extRemes::qevd(0.99, loc = p_max3_log10[1], 
                                scale = p_max3_log10[2], shape = p_max3_log10[3])
p_min3_log10 <- Min_extr_PC3_log10$results$par
qtl99_min3_log10<-extRemes::qevd(0.01, loc = p_min3_log10[1], 
                                scale = p_min3_log10[2], shape = p_min3_log10[3])

#Outliers
outliers_gev_pc1<-Min_PC1_log10<qtl99_min1_log10|Max_PC1_log10>qtl99_max1_log10
outliers_gev_pc2<-Min_PC2_log10<qtl99_min2_log10|Max_PC2_log10>qtl99_max2_log10
outliers_gev_pc3<-Min_PC3_log10<qtl99_min3_log10|Max_PC3_log10>qtl99_max3_log10

outliers_gev<-outliers_gev_pc1|outliers_gev_pc2|outliers_gev_pc3

#Plot the charts
pdf(file="./Figures/Log/PCs_scores_chart.pdf", width = 9, height = 9)
par(mfrow=c(3,1),mar = c(2, 2.5, 1, 1))
plot(c(entrecampos_log10_min[, 1],entrecampos_log10_max[,1]), 
        c(Min_PC1_log10,Max_PC1_log10), pch=20, xlab="",ylab="", 
        col=c("#00BE67"), cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
rect(
    xleft = entrecampos_log10_min[,1],
    xright = entrecampos_log10_max[,1],
    ybottom = Min_PC1_log10,
    ytop = Max_PC1_log10,
    density = 50
)
points(entrecampos_log10_min[which(Min_PC1_log10<qtl99_min1_log10),1],
        Min_PC1_log10[Min_PC1_log10<qtl99_min1_log10], pch=20, col="#ff0000")
points(entrecampos_log10_min[which(Max_PC1_log10>qtl99_max1_log10),1],
        Max_PC1_log10[Max_PC1_log10>qtl99_max1_log10], pch=20, col="#ff0000")
abline(h=c(qtl99_min1_log10,qtl99_max1_log10))
plot(c(entrecampos_log10_min[, 1],entrecampos_log10_max[,1]), 
        c(Min_PC2_log10,Max_PC2_log10), pch=20, xlab="",ylab="", 
        col=c("#00BE67"), cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
rect(
    xleft = entrecampos_log10_min[,1],
    xright = entrecampos_log10_max[,1],
    ybottom = Min_PC2_log10,
    ytop = Max_PC2_log10,
    density = 50
)
points(entrecampos_log10_min[which(Max_PC2_log10>qtl99_max2_log10),1],
        Max_PC2_log10[Max_PC2_log10>qtl99_max2_log10], pch=20, col="#ff0000")
points(entrecampos_log10_min[which(Min_PC2_log10<qtl99_min2_log10),1],
        Min_PC2_log10[Min_PC2_log10<qtl99_min2_log10], pch=20, col="#ff0000")
abline(h=c(qtl99_min2_log10,qtl99_max2_log10))
plot(c(entrecampos_log10_min[, 1],entrecampos_log10_max[,1]), 
        c(Min_PC3_log10,Max_PC3_log10), pch=20, xlab="",ylab="", 
        col=c("#00BE67"), cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
rect(
    xleft = entrecampos_log10_min[,1],
    xright = entrecampos_log10_max[,1],
    ybottom = Min_PC3_log10,
    ytop = Max_PC3_log10,
    density = 50
)
points(entrecampos_log10_min[which(Max_PC3_log10>qtl99_max3_log10),1],
        Max_PC3_log10[Max_PC3_log10>qtl99_max3_log10], pch=20, col="#ff0000")
points(entrecampos_log10_min[which(Min_PC3_log10<qtl99_min3_log10),1],
        Min_PC3_log10[Min_PC3_log10<qtl99_min3_log10], pch=20, col="#ff0000")
abline(h=c(qtl99_min3_log10,qtl99_max3_log10))
dev.off()

###
#Classic PCA on the daily average
pca_mean_log10<-prcomp(entrecampos_log10_mean_ma[,2:9])$x
pca_mean_pcs_log10<-prcomp(entrecampos_log10_mean_ma[,2:9])$rotation
fviz_eig(prcomp(entrecampos_log10_mean_ma[,2:9])) #choose 2 PCs
get_eig(prcomp(entrecampos_log10_mean_ma[,2:9])) #87.02% of explained variability with 2 PCs

#Outliers
outliers_pca_pc1<-pca_mean_log10[,1]>(3*sd(pca_mean_log10[,1]))|pca_mean_log10[,1]<(-3*sd(pca_mean_log10[,1]))
outliers_pca_pc2<-pca_mean_log10[,2]>(3*sd(pca_mean_log10[,2]))|pca_mean_log10[,2]<(-3*sd(pca_mean_log10[,2]))
outliers_pca<-outliers_pca_pc1|outliers_pca_pc2

#Plot the charts
pdf(file="./Figures/Log/PCm_chart_log.pdf", width = 12, height = 8)
par(mfrow=c(2,1), mar=c(2,2.5,1,1))
plot(entrecampos_log10_mean_ma[,1],pca_mean_log10[,1],type = "l",ylab="PC1m scores", 
        cex.lab=1.3, cex.axis=1.3, cex.main=1.3, cex.sub=1.3)
points(entrecampos_log10_mean_ma[,1],pca_mean_log10[,1], pch = 20, col=c("#00BE67"))
abline(h=c(mean(pca_mean_log10[,1]),3*sd(pca_mean_log10[,1]),-3*sd(pca_mean_log10[,1])))
points(entrecampos_log10_mean_ma[which(pca_mean_log10[,1]>(3*sd(pca_mean_log10[,1]))),1],
        pca_mean_log10[pca_mean_log10[,1]>(3*sd(pca_mean_log10[,1])),1], pch=20, col="#ff0000")
plot(entrecampos_log10_mean_ma[,1],pca_mean_log10[,2],type = "l",ylab="PC2m scores", 
        cex.lab=1.3, cex.axis=1.3, cex.main=1.3, cex.sub=1.3)
points(entrecampos_log10_mean_ma[,1],pca_mean_log10[,2], pch = 20, col=c("#00BE67"))
abline(h=c(mean(pca_mean_log10[,2]),3*sd(pca_mean_log10[,2]),-3*sd(pca_mean_log10[,2])))
points(entrecampos_log10_mean_ma[which(pca_mean_log10[,2]>(3*sd(pca_mean_log10[,2]))),1],
        pca_mean_log10[pca_mean_log10[,2]>(3*sd(pca_mean_log10[,2])),2], pch=20, col="#ff0000")
dev.off()

#Plot PCs
PCs_log10 <- data.frame(PC1s=pca1_sym_log10, PC1c=pca_mean_pcs_log10[,1],
                    PC2s=pca2_sym_log10, PC2c=pca_mean_pcs_log10[,2],
                    PC3s=pca3_sym_log10)
PCs_log10$id <- c("Sulphur Dioxide (µg/m3)",  "Particles < 10 µm (µg/m3)",  
                    "Ozone (µg/m3)", "Nitrogen Dioxide (µg/m3)",  
                    "Carbon Monoxide (µg/m3)", "Particles < 2.5 µm (µg/m3)", 
                    "Nitrogen Oxides (µg/m3)", "Nitrogen Monoxide (µg/m3)")
PCs_1col_log10 <- reshape2::melt(PCs_log10,id.vars = "id")
PCs_1col_log10$id<-factor(PCs_1col$id, levels = c("Carbon Monoxide (µg/m3)", 
                    "Nitrogen Dioxide (µg/m3)", "Nitrogen Monoxide (µg/m3)", 
                    "Nitrogen Oxides (µg/m3)",  "Ozone (µg/m3)",  
                    "Particles < 10 µm (µg/m3)", "Particles < 2.5 µm (µg/m3)",  
                    "Sulphur Dioxide (µg/m3)"))

pdf(file="./Figures/Log/PCs_barplot_log.pdf", width = 12, height = 4)
ggplot(PCs_1col_log10,aes(x=value, y = id, fill=id)) + 
    facet_wrap(~variable,nrow=1) + geom_bar(stat="identity") +
    theme(legend.position = "none",axis.title.x=element_blank(),
        axis.title.y=element_blank(), text = element_text(size=20))
dev.off()

###
# ROBPCA
###

#Obtain ROBPCA
robpca_mean_log10 <- PcaHubert(entrecampos_log10_mean_ma[,2:9], crit.pca.distances=0.99)

#Outliers
outliers_robpca<-!robpca_mean_log10$flag

#Outlier map
pdf(file="./Figures/Log/ROBPCA_log.pdf", width = 6, height = 6)
col = rep('black',1096); col[robpca_mean_log10@flag==FALSE] = 'red'
plot(robpca_mean_log10, pch=20, col=col)
dev.off()

###
# Isolation Forest
###

iso_forest_mean_log10 <- isolation.forest(
    entrecampos_log10_mean_ma[,2:9],
    ndim=1, sample_size=256,
    ntrees=100,
    missing_action="fail")

#Obtain outliers scores
pred_iso_forest_mean_log10 <- predict(iso_forest_mean_log10, entrecampos_log10_mean_ma[,2:9])

#Outliers
outliers_iso_forest<-pred_iso_forest_mean_log10>0.55

#Plot outlier scores of each observation
pdf(file="./Figures/Log/iso_forest_log.pdf", width = 6, height = 6)
col = rep('black',1096); col[outliers_iso_forest] = 'red'
plot(entrecampos_log10_mean_ma[,1],pred_iso_forest_mean_log10,
    pch=20,col=col,xlab="Observations",ylab="Outlier Score", main="Isolation Forest", cex.lab=1.3, cex.axis=1.3, cex.main=1.3, cex.sub=1.3)
abline(h=0.55)
dev.off()

###
# Plot outliers from all methods
###

outliers_all<-outliers_gev|outliers_pca|outliers_robpca|outliers_iso_forest
outliers_total<-as.data.frame(cbind(outliers_gev,outliers_pca,outliers_robpca,outliers_iso_forest))
outliers_total<-cbind("date"=entrecampos_log10_max_ma[,1],outliers_total)

outliers_total_true<-outliers_total[outliers_all,]

outliers_heatmap <- 1*(outliers_total[,-1]==TRUE)
rownames(outliers_heatmap)<-outliers_total[,1]
outliers_heatmap[outliers_gev_pc3,1]<-4
outliers_heatmap[outliers_gev_pc2,1]<-3
outliers_heatmap[outliers_gev_pc1,1]<-2
outliers_heatmap[outliers_gev_pc1 & outliers_gev_pc2,1]<-2.5
outliers_heatmap[outliers_pca_pc1,2]<-2
outliers_heatmap[outliers_pca_pc2,2]<-3
outliers_heatmap<-outliers_heatmap[outliers_all,]

pdf(file="./Figures/Log/Outliers_Log_Total.pdf", width = 6.5, height = 25)
heatmap.2(outliers_heatmap, 
        scale = "none", Rowv = NA, Colv = NA, dendrogram = "none",
        col=c("#eceeef", "#00B0F6","#F8766D", "#A3A500", "#E76BF3", "#00BF7D"),
        margins=c(2,19),srtCol=0,labRow = outliers_total_true[,1],
        labCol = c("Symbolic\nPCA", "Conventional\nPCA", "ROBPCA", "Isolation\nForest"),
        cexRow = 1,cexCol = 1,adjCol=c(0.5,0.5),
        key=FALSE, trace = "none", density.info = "none",
        lwid = c(1,50), lhei=c(1, 60), colsep=1:4, rowsep=1:156)
legend(x="topright",pch=15,pt.cex=2,
        col=c("#00B0F6","#F8766D", "#A3A500", "#E76BF3", "#00BF7D"),
        legend = c("Outliers","Outliers - PC1", "Outliers - PC1&2", 
                    "Outliers - PC2", "Outliers - PC3"))
dev.off()