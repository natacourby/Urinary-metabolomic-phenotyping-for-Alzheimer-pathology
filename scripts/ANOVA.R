#########################################
# Load metabolomic matrices 
#########################################
library(dplyr)
library(ggpubr)
source("./scripts/feature_selection.R")

# Load data
annotated_metabolites <- read.csv("./results/annotated_metabolites.txt",sep="\t", header=T)
covariates <- c("Gender","Centre","Age","Diagnosis")
data <- feature_selection_data_preparation("D","I")
data_annotated<- data[,colnames(data) %in% c(as.character(annotated_metabolites$Metabolic_features),covariates)]
data_annotated$Diagnosis <- as.character(data_annotated$Diagnosis)
data_annotated$Diagnosis[data_annotated$Diagnosis=="ADC"] <- "AD"
data_annotated$Diagnosis <- as.factor(data_annotated$Diagnosis)

result <- data.frame("Metabolite"="","ANOVA_p_value"=0,"ANOVA_F_value"=0,
                     "TukeyHSD_cMCI_AD"=0,"TukeyHSD_CTL_AD"=0,"TukeyHSD_sMCI_AD"=0,
                     "TukeyHSD_CTL_cMCI"=0,"TukeyHSD_sMCI_cMCI"=0,"TukeyHSD_sMCI_CTL"=0,
                     "MANOVA_p_value"=0,"MANOVA_F_value"=0,"MANOVA_Pillai_trace"=0)
# Select metabolite
for (m in unique(annotated_metabolites$Metabolite)){
  selected_metabolic_features<- c(as.character(annotated_metabolites[annotated_metabolites$Metabolite==m,]$Metabolic_features))
  expr <- data_annotated[,colnames(data_annotated) %in% c(selected_metabolic_features,covariates)]
  
  data_anova <- data.frame(value=rowMeans(expr[colnames(expr) %in% selected_metabolic_features]),group=expr$Diagnosis)
  grouped <- group_by(data_anova,group)
  
  # Mean, SD and #N of samples per diagnosis
  summarise(grouped,mean=mean(value),sd=sd(value), nr=n())

  # Compute the analysis of variance (ANOVA)
  res.aov <- aov(value ~ group, data = data_anova)
  
  # Summary of the analysis
  summary(res.aov)

  #Tukey multiple comparisons of means
  TukeyHSD(res.aov)
  
  TurkeyHSD_res <- TukeyHSD(res.aov)$group[,4]
  
  # MANOVA
  if (length(selected_metabolic_features)>1) {
    Y <- cbind(expr[colnames(expr) %in% selected_metabolic_features])
    Y <- sapply( Y, as.numeric )
    group <- expr$Diagnosis
  
    fit <- manova(Y ~ group)
    summary(fit, test="Pillai")
    metabolite_manova_result <- c(summary(fit)$stats[, "Pr(>F)"][[1]],summary(fit)$stats[, "approx F"][[1]],summary(fit)$stats[, "Pillai"][[1]])
  }
  else{
    metabolite_manova_result <- c(0,0,0)
  }
  metabolite_result <- data.frame("Metabolite"=m,
                                  "ANOVA_p_value"=summary(res.aov)[[1]][["Pr(>F)"]][[1]],
                                  "ANOVA_F_value"=  summary(res.aov)[[1]][["F value"]][[1]],
                                  "TukeyHSD_cMCI_AD"=TurkeyHSD_res[1],
                                  "TukeyHSD_CTL_AD"=TurkeyHSD_res[2],
                                  "TukeyHSD_sMCI_AD"=TurkeyHSD_res[3],
                                  "TukeyHSD_CTL_cMCI"=TurkeyHSD_res[4],
                                  "TukeyHSD_sMCI_cMCI"=TurkeyHSD_res[5],
                                  "TukeyHSD_sMCI_CTL"=TurkeyHSD_res[6],
                                  "MANOVA_p_value"=metabolite_manova_result[1],
                                  "MANOVA_F_value"=metabolite_manova_result[2],
                                  "MANOVA_Pillai_trace"=metabolite_manova_result[3])
  
  result <- rbind(result, metabolite_result)
  # Plot averages per diagnosis
  png(paste("./images/annotated_metabolites/",m,".png",sep=""), width = 5, height = 5, units = 'in', res = 300)
  ggboxplot(data_anova, x = "group", y = "value", 
            color = "group", palette = c("#00AFBB", "#E7B800", "#FC4E07","#4F7942"),
            order = c("AD", "cMCI", "CTL","sMCI"),
            ylab = m, xlab = "Diagnosis")
  dev.off()
  
  #Another version of the plot
  #require(ggplot2)
  # ggplot(data = data_anova, aes(x=group, y=value)) + 
  #   geom_boxplot(aes(fill=group)) + 
  #   labs(y=m, x = "Diagnosis")
}

result <- result[-1,]
rownames(result) <- unique(annotated_metabolites$Metabolite)
write.csv(result,"./results/annotated_metabolites_ANOVA.csv",row.names = F)





