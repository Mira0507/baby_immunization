library(ggplot2)
library(broom)
library(tidyverse)
library(data.table)
library(gridExtra)
library(pheatmap)
library(factoextra)

##################################### Importing and Cleaning Data ##################################### 
# importing a data file
nut <- fread("HNP_StatsData.csv", header = TRUE) 

names(nut) <- str_replace_all(colnames(nut), " ", "_")

nut <- nut[, -c("Indicator_Code",
                "V65")]        

# subsetting data table w indicator of interest 
indicators <- unique(nut$Indicator_Name)

IndicatorOfInterest <- c(
        "Mortality rate, neonatal (per 1,000 live births)",
        "Mortality rate, infant (per 1,000 live births)",
        "Mortality rate, under-5 (per 1,000)",
        "Incidence of tuberculosis (per 100,000 people)",
        "Immunization, BCG (% of one-year-old children)",
        "Immunization, DPT (% of children ages 12-23 months)",
        "Immunization, measles (% of children ages 12-23 months)",
        "Immunization, Pol3 (% of one-year-old children)")

IndicatorFilter <- nut$Indicator_Name %in% IndicatorOfInterest

nut1 <- nut[IndicatorFilter, -as.character(1960:1979)] %>%
        gather(Year, Value, -c("Country_Name", 
                               "Country_Code", 
                               "Indicator_Name")) %>% 
        spread(Indicator_Name, Value) 

NewNames <- c("Country_name",
              "Country_code",
              "Year",
              "Immunization_BCG",
              "Immunization_DPT",
              "Immunization_measles",
              "Immunization_Pol3",
              "Incidence_of_tuberculosis",
              "Mortality_rate_infant",
              "Mortality_rate_neonatal",
              "Mortality_rate_under5")

names(nut1) <- NewNames



# Removing countries that are NOT countries
# e.g. "World", "OECD", etc
region <- c("World",
            "East Asia",
            "income",
            "Europe",
            "Middle",
            "East",
            "North", 
            "OECD",
            "Pacific",
            "Post",
            "Sub-Saharan",
            "small",
            "countries",
            "Fragile",
            "Islands",
            "demographic",
            "Latin",
            "area",
            "Asia")

for (x in region) {
        nut1 <- nut1[!str_detect(nut1$Country_name, x), ]
}


nut2 <- nut1[complete.cases(nut1), ]

# nesting 
nut2_nst <- nut2 %>% 
        group_by(Country_name) %>%
        nest() 



################################ Supervised Learning: Linear Regression ################################ 

# Setting formulae 
FormulaInfant <- as.formula(Mortality_rate_infant ~ 
                                    Immunization_BCG + 
                                    Immunization_DPT + 
                                    Immunization_measles +
                                    Immunization_Pol3)
FormulaNeonatal <- as.formula(Mortality_rate_neonatal ~ 
                                      Immunization_BCG + 
                                      Immunization_DPT + 
                                      Immunization_measles +
                                      Immunization_Pol3) 
FormulaUnder5 <- as.formula(Mortality_rate_under5 ~ 
                                    Immunization_BCG + 
                                    Immunization_DPT + 
                                    Immunization_measles +
                                    Immunization_Pol3) 


# constructing linear regression models, calculating fit (Fit_*), 
# predicting (Pred_*), and calculating coefficients (Coef_*) by country
nut3_nst <- nut2_nst %>%
        mutate(Model_neonatal = map(data, ~ lm(FormulaNeonatal, .x)),
               Model_infant = map(data, ~ lm(FormulaInfant, .x)),
               Model_under5 = map(data, ~ lm(FormulaUnder5, .x)),
               Fit_neonatal = map(Model_neonatal, ~ glance(.x)),
               Fit_infant = map(Model_infant, ~ glance(.x)),
               Fit_under5 = map(Model_under5, ~ glance(.x)),
               Pred_neonatal = map(Model_neonatal, ~ augment(.x)),
               Pred_infant = map(Model_infant, ~ augment(.x)),
               Pred_under5 = map(Model_under5, ~ augment(.x)),
               Coef_neonatal = map(Model_neonatal, ~ tidy(.x)),
               Coef_infant = map(Model_infant, ~ tidy(.x)),
               Coef_under5 = map(Model_under5, ~ tidy(.x))) %>%
        as.data.table()


# data cleaning for squared distribution across the countries 
nut3_lmR2 <- nut3_nst[, c("Country_name", 
                          "Fit_neonatal",
                          "Fit_infant",
                          "Fit_under5")] %>%
        unnest(Fit_neonatal, Fit_infant, Fit_under5) %>% 
        select(Country_name, r.squared, r.squared1, r.squared2) %>%
        rename(Neonatal = r.squared,
               Infant = r.squared1, 
               Under5 = r.squared2) %>%
        gather(key = "Mortality_Group", 
               value = R2, -Country_name) %>%
        mutate(Mortality_Group = factor(Mortality_Group,
                                        levels = c("Neonatal",
                                                   "Infant",
                                                   "Under5"))) %>%
        as.data.table()


# Plotting distribution of Rsquared 
nut3_lmR2_plot <- ggplot(nut3_lmR2, 
                         aes(x = Mortality_Group, 
                             y = R2, 
                             color = Mortality_Group, 
                             fill = Mortality_Group)) +
        geom_boxplot(alpha = 0.5) +
        theme_bw() + 
        ggtitle("Distribution of Rsquared in Linear Regression Models") + 
        xlab("Mortality Group") + 
        ylab("Rsquared")

# Top and bottom Rsquared countries 
R2ToptoBottom <- nut3_lmR2 %>% 
        arrange(desc(R2)) %>%
        as.data.table()

# Filtering countries whose R2 from Under
Neonatal_R2Over0.9 <- R2ToptoBottom[R2 >= 0.9 & Mortality_Group == "Neonatal"]
Infant_R2Over0.9 <- R2ToptoBottom[R2 >= 0.9 & Mortality_Group == "Infant"]
Under5_R2Over0.9 <- R2ToptoBottom[R2 >= 0.9 & Mortality_Group == "Under5"]


# Data cleaning for regression coefficients 
nut3_Coef <- nut3_nst[, .(Country_name, 
                          Coef_neonatal,
                          Coef_infant,
                          Coef_under5)] %>%
        unnest(Coef_neonatal,
               Coef_infant,
               Coef_under5) %>%
        select(Country_name, term, starts_with("estimate")) %>%
        rename(Est_neonatal = estimate, 
               Est_infant = estimate1, 
               Est_under5 = estimate2) 




################################ Unsupervised Learning: Clustering ################################ 


# Imputation

library(bnstruct)


CleanData_fn <- function(est, group) {
        
        # data cleaning
        dt <- nut3_Coef[, c("Country_name", "term", est)] %>%
                spread(term, est) %>%
                inner_join(nut3_lmR2[Mortality_Group == group,
                                     c("Country_name", "R2")],
                           by = "Country_name") %>%
                rename(Intercept = "(Intercept)") 
        
        # converting numerical columns to matrix
        dt_mtx <- as.matrix(dt[, -1])
        
        # knn imputation 
        dt_imp <- knn.impute(dt_mtx)
        
        # combining imputed matrix with country names 
        dt1 <- cbind(dt[, 1], as.data.table(dt_imp)) %>%
                mutate(Country_name = ifelse(Country_name == "Korea, Dem. Peopleâ€™s Rep.",
                                             "Korea, Dem. Peoples Rep.", 
                                             Country_name),
                       Group = group) %>%
                as.data.table()
        
        
        return(dt1)
}

# Imputated Coefficients Table        
Coef_neonatal <- CleanData_fn("Est_neonatal", "Neonatal")
Coef_infant <- CleanData_fn("Est_infant", "Infant")
Coef_under5 <- CleanData_fn("Est_under5", "Under5")

CoefTable <- rbind(Coef_neonatal, Coef_infant, Coef_under5) %>%
        nest(1:7)


####################### tSNE


# Running tSNE
library(Rtsne)
tSNE_fn <- function(dt, pp) {
        
        dt1 <- dt[, -c("Country_name")]
        
        # row labeling 
        rownames(dt1) <- dt$Country_name
        dt1 <- as.matrix(dt1)
        
        # tSNE
        set.seed(113)
        Rtsne(dt1,
              PCA = T, 
              perplexity = pp, 
              max_iter = 2000)
}


iterPlot_fn <- function(dt) {
        
        plot(dt$itercosts, 
             type = "l",
             ylab = "Total K-L Divergence Cost",
             xlab = "Gradient Descent (50 Steps Each)",
             main = "Optimal Number of Iterations")
        
}



CoefTable0 <- CoefTable %>%
        
        # Running tSNE with perplexity 5, 7, 10
        mutate(tsne5 = map(data, ~ tSNE_fn(.x, 5)),
               tsne7 = map(data, ~ tSNE_fn(.x, 7)),
               tsne10 = map(data, ~ tSNE_fn(.x, 10))) 

CoefTable1 <- CoefTable0 %>%
        
        # Extracting coordinates when perplexity = 5, 7, 10
        mutate(data = map2(data, tsne5, ~ mutate(.x, tsne5X = .y$Y[, 1])),
               data = map2(data, tsne5, ~ mutate(.x, tsne5Y = .y$Y[, 2])),
               data = map2(data, tsne7, ~ mutate(.x, tsne7X = .y$Y[, 1])),
               data = map2(data, tsne7, ~ mutate(.x, tsne7Y = .y$Y[, 2])),
               data = map2(data, tsne10, ~ mutate(.x, tsne10X = .y$Y[, 1])),
               data = map2(data, tsne10, ~ mutate(.x, tsne10Y = .y$Y[, 2]))) 

# Determining perplexity 
CoefTable_ComparePP <- CoefTable1 %>%
        
        # data cleaning
        unnest(data) %>%
        select(Country_name, Group, ends_with(c("X", "Y")))


ComparePP_Plot_fn <- function(dt, xvar, yvar, group, tit) {
        
        ggplot(dt, aes(x = xvar, y = yvar, color = group)) + 
                geom_point(alpha = 0.5, size = 1.5) +
                facet_grid(Group ~ .) + 
                theme_bw() + 
                theme(legend.position = "none") + 
                labs(title = tit,
                     x = "X", 
                     y = "Y") + 
                xlim(-50, 50) + 
                ylim(-50, 50)
}

grid.arrange(ComparePP_Plot_fn(CoefTable_ComparePP, 
                               CoefTable_ComparePP$tsne5X,
                               CoefTable_ComparePP$tsne5Y,
                               CoefTable_ComparePP$Group, 
                               "tSNE (Perplexity = 5)"),
             ComparePP_Plot_fn(CoefTable_ComparePP, 
                               CoefTable_ComparePP$tsne7X,
                               CoefTable_ComparePP$tsne7Y,
                               CoefTable_ComparePP$Group, 
                               "tSNE (Perplexity = 7)"),
             ComparePP_Plot_fn(CoefTable_ComparePP, 
                               CoefTable_ComparePP$tsne10X,
                               CoefTable_ComparePP$tsne10Y,
                               CoefTable_ComparePP$Group, 
                               "tSNE (Perplexity = 10)"),
             nrow = 1) # Perplexity = 7 is the best



# Determining k: scree plot
Scree_fn <- function(dt, tit, intercept) {
        
        # measuring total within ss 
        set.seed(1987)
        ttWithinss <- map_dbl(1:10, 
                              function(k) {
                                      km <- kmeans(x = dt, 
                                                   centers = k,
                                                   nstart = 25)
                                      km$tot.withinss})
        screeDT <- data.table(
                k = 1:10, 
                Total_Within_SS = ttWithinss
        )
        
        # creating a scree plot
        ggplot(screeDT,
               aes(x = k,
                   y = Total_Within_SS)) + 
                geom_line(size = 1, color = "blue") +
                geom_point(size = 2, color = "blue") +
                ggtitle(tit) +
                theme_bw() + 
                ylab("Total Within Cluster Sum of Squares") + 
                scale_x_continuous(breaks = 1:10, 
                                   minor_breaks = NULL) +
                geom_vline(xintercept = intercept, 
                           size = 1,
                           color = "red")
}

CoefTable2 <- CoefTable1 %>%
        
        # data cleaning
        mutate(data = map(data, ~ column_to_rownames(.x, 
                                                     var = "Country_name")),
               data = map(data, ~ select(.x, 
                                         starts_with("tsne7"))),
               # Heatmap
               HeatMap = map2(data, Group,
                              ~ pheatmap(.x, main = .y)),
               
               # Scree plot to determine k
               KMScree = map(data, 
                             ~ Scree_fn(.x, paste("Scree Plot:", Group), 4)))

CoefTable3 <- CoefTable2 %>%
               
               mutate(# clustering k = 4, 5, 6
                       data = map(data, 
                                  ~ mutate(.x, kmCluster4 = factor(kmeans(.x, centers = 4, 
                                                                          nstart = 25)$cluster))),
                       data = map(data, 
                                  ~ mutate(.x, kmCluster5 = factor(kmeans(.x, centers = 5, 
                                                                          nstart = 25)$cluster))),
                       
                       data = map(data, 
                                  ~ mutate(.x, kmCluster6 = factor(kmeans(.x, centers = 6, 
                                                                          nstart = 25)$cluster))),
                       data = map(data, 
                                  ~ mutate(.x, hCluster4 = factor(cutree(hclust(dist(.x), 
                                                                                method = "average"),
                                                                         k = 4)))),
                       data = map(data, 
                                  ~ mutate(.x, hCluster5 = factor(cutree(hclust(dist(.x), 
                                                                                method = "average"),
                                                                         k = 5)))),
                       data = map(data, 
                                  ~ mutate(.x, hCluster6 = factor(cutree(hclust(dist(.x), 
                                                                                method = "average"),
                                                                         k = 6)))),
                       # adding country names back 
                       data = map(data,
                                  ~ mutate(.x, Country_name = Coef_neonatal$Country_name))) %>%
        unnest(data)

# plotting tSNE-based clustering
tSNEClustering_Viz_fn <- function(df, hc, kmc, tit) {
        
        # data cleaning
        gather(df, Clustered_by, Cluster, c(hc, kmc)) %>%
                mutate(Clustered_by = ifelse(Clustered_by == hc, 
                                             "Hierarchical", 
                                             "K-means")) %>%
                
                # plotting
                ggplot(aes(x = tsne7X, y = tsne7Y, color = Cluster)) + 
                geom_point(alpha = 0.5, size = 1.5) +
                facet_grid(Group ~ Clustered_by) + 
                theme_bw() +
                ggtitle(tit) + 
                xlab("X") + 
                ylab("Y")
        
}

# k = 4
tSNE_k4_plot <- tSNEClustering_Viz_fn(CoefTable3, 
                                      "hCluster4", 
                                      "kmCluster4", 
                                      "Clustering and Visualization by tSNE (k = 4)")

# k = 5 (looks like the best)
tSNE_k5_plot <- tSNEClustering_Viz_fn(CoefTable3, 
                                      "hCluster5", 
                                      "kmCluster5", 
                                      "Clustering and Visualization by tSNE (k = 5)")

# k = 6
tSNE_k6_plot <- tSNEClustering_Viz_fn(CoefTable3, 
                                      "hCluster6", 
                                      "kmCluster6", 
                                      "Clustering and Visualization by tSNE (k = 6)")


####################### PCA

# running PCA



PCA_fn <- function(dt) {
        
        dt1 <- dt[, -1]
        
        # row labeling 
        rownames(dt1) <- dt$Country_name
        
        # PCA 
        prcomp(dt1, 
               scale = TRUE, 
               center = TRUE)
}

# Scree plots
PCAScree_fn <- function(pca, pc, tit) {
        
        # Computing cumulative proportion of var explained
        v <- pca$sdev^2
        pve <- v / sum(v)
        
        
        dt <- data.table(CumSum = cumsum(pve) * 100,
                         PC = pc)
        
        # plotting
        ggplot(dt, aes(x = PC, 
                       y = CumSum)) + 
                geom_line(size = 1, color = "blue") + 
                geom_point(size = 1.5, color = "blue") + 
                theme_bw() +
                scale_x_continuous(n.breaks = length(pc)) + 
                ylab("Cumulative Proportion of Variance Explained (%)") + 
                ggtitle(tit)
}

CoefTable4 <- CoefTable %>%
        
        # PCA
        mutate(PCA = map(data, ~ PCA_fn(.x)),
               
               # biplots
               biplot1 = map2(PCA, Group,~ fviz_pca_var(.x, col.var="contrib",
                                                gradient.cols = c("blue", "red"),
                                                repel = TRUE,
                                                title = .y)),
               biplot2 = map2(PCA, Group, ~ fviz_pca_biplot(.x, 
                                                    repel = TRUE,
                                                    geom.ind = "point",
                                                    title = .y)),
               
               # scree plots
               PCA_scree = map2(PCA, Group, ~ PCAScree_fn(.x, 1:6, .y)),
               
               # extracting PC1-4 coordinates
               data = map2(data, PCA, ~ mutate(.x, PC1 = .y$x[, 1])),
               data = map2(data, PCA, ~ mutate(.x, PC2 = .y$x[, 2])),
               data = map2(data, PCA, ~ mutate(.x, PC3 = .y$x[, 3])),
               data = map2(data, PCA, ~ mutate(.x, PC4 = .y$x[, 4])))

CoefTable5 <- CoefTable4 %>%
        
        # data cleaning
        mutate(data = map(data, ~ column_to_rownames(.x, var = "Country_name")),
               data = map(data, ~ select(.x, starts_with("PC"))),
               
               # heatmap
               HeatMap = map2(data, Group, ~ pheatmap(.x, main = .y)),
               
               # Scree plot to determine k
               KMScree = map2(data, Group,
                             ~ Scree_fn(.x, paste("Scree Plot:", .y), 5)))


CoefTable6 <- CoefTable5 %>%
        
        # clustering k = 5, 6, 7
        mutate(data = map(data, 
                          ~ mutate(.x, kmCluster5 = factor(kmeans(.x, centers = 5, 
                                                                  nstart = 25)$cluster))),
               
               data = map(data, 
                          ~ mutate(.x, kmCluster6 = factor(kmeans(.x, centers = 6, 
                                                                  nstart = 25)$cluster))),
               
               data = map(data, 
                          ~ mutate(.x, kmCluster7 = factor(kmeans(.x, centers = 7, 
                                                                  nstart = 25)$cluster))),
               
               data = map(data, 
                          ~ mutate(.x, hCluster5 = factor(cutree(hclust(dist(.x), 
                                                                        method = "average"), 
                                                                 k = 5)))),
               
               data = map(data, 
                          ~ mutate(.x, hCluster6 = factor(cutree(hclust(dist(.x), 
                                                                        method = "average"), 
                                                                 k = 6)))),
               
               data = map(data, 
                          ~ mutate(.x, hCluster7 = factor(cutree(hclust(dist(.x), 
                                                                        method = "average"), 
                                                                 k = 7)))),
               
               # adding country names back 
               data = map(data,
                          ~ mutate(.x, Country_name = Coef_neonatal$Country_name))) %>%
        unnest(data)
                                   
                                                          
# plotting PCA-based clustering
PCAClustering_Viz_fn <- function(df, hc, kmc, tit) {
        
        # data cleaning
        gather(df, Clustered_by, Cluster, c(hc, kmc)) %>%
                mutate(Clustered_by = ifelse(Clustered_by == hc, 
                                             "Hierarchical", 
                                             "K-means")) %>%
                
                # plotting
                ggplot(aes(x = PC1, y = PC2, color = Cluster)) + 
                geom_point(alpha = 0.5, size = 1.5) +
                facet_grid(Group ~ Clustered_by) + 
                theme_bw() +
                ggtitle(tit) + 
                xlab("PC1") + 
                ylab("PC2")
        
}

# k = 5
PCA_k5_plot <- PCAClustering_Viz_fn(CoefTable6,
                                    "hCluster5",
                                    "kmCluster5",
                                    "Clustering and Visualization by PCA (k = 5)")
# k = 6
PCA_k6_plot <- PCAClustering_Viz_fn(CoefTable6,
                                    "hCluster6",
                                    "kmCluster6",
                                    "Clustering and Visualization by PCA (k = 6)")
# k = 7
PCA_k7_plot <- PCAClustering_Viz_fn(CoefTable6,
                                    "hCluster7",
                                    "kmCluster7",
                                    "Clustering and Visualization by PCA (k = 7)")

                             
################################ Data Interpretation Based on tSNE-clustering ################################


# Data cleaning: Combining coefficient values with cluster results 
CoefTable7 <- CoefTable1[, c("Group", "data")] %>%
        
        mutate(data = map(data, ~ select(.x, !starts_with("tsne")))) %>%
        
        unnest(data) %>%
        
        inner_join(CoefTable3[, c("Group", 
                                  "kmCluster5",
                                  "hCluster5",
                                  "Country_name")],
                   by = c("Group", "Country_name")) %>%
        
        as.data.table()


####################### Distribution of coefficient by cluster


CoefDist_fn <- function(dt, column_compare, method){
        
        # data cleaning
        dt %>% gather(Coefficient, 
                      Value, 
                      -c("Country_name", column_compare, "Group")) %>%
                
                mutate(Group = factor(Group,
                                      levels = c("Neonatal",
                                                 "Infant",
                                                 "Under5"))) %>%
                
                rename(Cluster = column_compare) %>%
                
                # plotting
                ggplot(aes(x = Coefficient, 
                           y = Value,
                           fill = Coefficient)) + 
                geom_boxplot() + 
                theme_bw() + 
                theme(axis.text.x = element_blank()) + 
                facet_grid(Cluster ~ Group, scales = "free_y") + 
                ggtitle(paste("Distribution of Coefficients:", method))}
        
        

# Plotting distribution of coefficient by cluster
DistCoef_HC_plot1 <- CoefDist_fn(CoefTable7[, !c("kmCluster5", "R2")], 
                                          "hCluster5",
                                          "Hierarchical Clustering")     

DistCoef_HC_plot2 <- DistCoef_HC_plot1 +
        ylab("Value (Log Scale)") +
        scale_y_log10() 
        

DistCoef_KMC_plot1 <- CoefDist_fn(CoefTable7[, !c("hCluster5", "R2")], 
                                  "kmCluster5",
                                  "K-means Clustering") 

DistCoef_KMC_plot2 <- DistCoef_KMC_plot1 +
        ylab("Value (Log Scale)") +
        scale_y_log10() 


####################### Distribution of R2 by cluster

                      
DistR2_plot <- CoefTable7 %>%
        
        # data cleaning
        gather(Clustering, Cluster, c("kmCluster5", "hCluster5")) %>%
        mutate(Clustering = ifelse(Clustering == "kmCluster5", "K-means", "Hierarchical")) %>%
        
        # plotting
        ggplot(aes(x = Cluster, y = R2, fill = Cluster)) +
        geom_boxplot() + 
        facet_grid(Clustering ~ Group) + 
        theme_bw()