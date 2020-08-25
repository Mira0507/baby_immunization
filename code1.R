library(ggplot2)
library(broom)
library(tidyverse)
library(data.table)
library(gridExtra)

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
library(factoextra)

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
                                             Country_name))
        
        
        return(dt1)
}

# Imputated Coefficients Table        
Coef_neonatal <- CleanData_fn("Est_neonatal", "Neonatal")
Coef_infant <- CleanData_fn("Est_infant", "Infant")
Coef_under5 <- CleanData_fn("Est_under5", "Under5")



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



pca_neonatal <- PCA_fn(Coef_neonatal)
pca_infant <- PCA_fn(Coef_infant)
pca_under5 <- PCA_fn(Coef_under5)

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


grid.arrange(PCAScree_fn(pca_neonatal, 1:6, "PCA: Neonatal"),
             PCAScree_fn(pca_infant, 1:6, "PCA: Infant"),
             PCAScree_fn(pca_under5, 1:6, "PCA: Under5"),
             nrow = 1)

# bioplots 
fviz_pca_biplot(pca_neonatal)
fviz_pca_biplot(pca_infant)
fviz_pca_biplot(pca_under5)


# Find a PC where cumulative proportion gets over 90%: PC4  
summary(pca_neonatal)
summary(pca_infant)
summary(pca_under5)

# Extracting PC1-4 coordinates 
pcaX_neonatal <- pca_neonatal$x[, 1:4]
pcaX_infant <- pca_infant$x[, 1:4]
pcaX_under5 <- pca_under5$x[, 1:4]


# hierarchical clustering & heatmap 
library(pheatmap)
PCA_Heatmap_Neonatal <- pheatmap(pcaX_neonatal,
                                 main = "Immunization and Mortality Rate: Neonatal") 
PCA_Heatmap_Infant <- pheatmap(pcaX_infant,
                               main = "Immunization and Mortality Rate: Infant")
PCA_Heatmap_Under5 <- pheatmap(pcaX_under5,
                               main = "Immunization and Mortality Rate: Under5")



# Distance & clustering (manually) 
hc_neonatal <- hclust(dist(pcaX_neonatal), 
                        method = "average")

hc_infant <- hclust(dist(pcaX_infant), 
                    method = "average")

hc_under5 <- hclust(dist(pcaX_under5), 
                    method = "average")

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
    


# Scree plots 
grid.arrange(Scree_fn(pcaX_neonatal, "K-means: Neonatal", 6),
             Scree_fn(pcaX_infant, "K-means: Infant", 5), 
             Scree_fn(pcaX_under5, "K-means: Under5", 5),
             nrow = 1)

# Extracting clustering results             
Extract_KMCluster_fn <- function(dt, k) {
        
        set.seed(1987)
        factor(kmeans(x = dt, 
                      centers = k, 
                      nstart = 25)$cluster)
}


# Combining PC1/2 coordinates with clustering results 
PCA_Clustering <- data.table(Country_name = Coef_neonatal$Country_name,
                             HCluster_Neonatal = factor(cutree(hc_neonatal, k = 6)),
                             HCluster_Infant = factor(cutree(hc_infant, k = 5)),
                             HCluster_Under5 = factor(cutree(hc_infant, k = 5)),
                             KMCluster_Neonatal = Extract_KMCluster_fn(pcaX_neonatal, 6),
                             KMCluster_Infant = Extract_KMCluster_fn(pcaX_infant, 5),
                             KMCluster_Under5 = Extract_KMCluster_fn(pcaX_under5, 5)) %>%
        
        mutate(Neonatal_R2 = ifelse(Country_name %chin% Neonatal_R2Over0.9$Country_name, 
                                    ">= 0.9", "< 0.9"),
               Infant_R2 = ifelse(Country_name %chin% Infant_R2Over0.9$Country_name, 
                                  ">= 0.9", "< 0.9"),
               Under5_R2 = ifelse(Country_name %chin% Under5_R2Over0.9$Country_name, 
                                  ">= 0.9", "< 0.9"),
               PC1_Neonatal = pca_neonatal$x[, 1],
               PC2_Neonatal = pca_neonatal$x[, 2],
               PC1_Infant = pca_infant$x[, 1],
               PC2_Infant = pca_infant$x[, 2],
               PC1_Under5 = pca_under5$x[, 1],
               PC2_Under5 = pca_under5$x[, 2])

PCA_Clustering <- PCA_Clustering[, c("Neonatal_R2",
                                     "Infant_R2",
                                     "Under5_R2") := 
                                         .(factor(Neonatal_R2, levels = c(">= 0.9", "< 0.9")),
                                           factor(Infant_R2, levels = c(">= 0.9", "< 0.9")),
                                           factor(Under5_R2, levels = c(">= 0.9", "< 0.9")))]


# Clustering visualization
ggplot(PCA_Clustering, 
       aes(x = PC1_Neonatal,
           y = PC2_Neonatal,
           shape = Neonatal_R2,
           color = KMCluster_Neonatal)) + 
        geom_point(alpha = 0.5, size = 2) + 
        facet_grid(KMCluster_Neonatal ~ Neonatal_R2) +
        theme_bw() + 
        ggtitle("K-means Clustering: Neonatal")


ggplot(PCA_Clustering, 
       aes(x = PC1_Neonatal,
           y = PC2_Neonatal,
           shape = Neonatal_R2,
           color = HCluster_Neonatal)) + 
        geom_point(alpha = 0.5, size = 2) + 
        facet_grid(HCluster_Neonatal ~ Neonatal_R2) + 
        theme_bw() + 
        ggtitle("Hierarchical Clustering: Neonatal")


ggplot(PCA_Clustering, 
       aes(x = PC1_Infant,
           y = PC2_Infant,
           shape = Infant_R2,
           color = KMCluster_Infant)) + 
        geom_point(alpha = 0.5, size = 2) + 
        facet_grid(KMCluster_Infant ~ Infant_R2) +
        theme_bw() + 
        ggtitle("K-means Clustering: Infant")

ggplot(PCA_Clustering, 
       aes(x = PC1_Infant,
           y = PC2_Infant,
           shape = Infant_R2,
           color = HCluster_Infant)) + 
        geom_point(alpha = 0.5, size = 2) + 
        facet_grid(HCluster_Infant ~ Infant_R2) +
        theme_bw() + 
        ggtitle("Hierarchical Clustering: Infant")

ggplot(PCA_Clustering, 
       aes(x = PC1_Under5,
           y = PC2_Under5,
           shape = Under5_R2,
           color = KMCluster_Under5)) + 
        geom_point(alpha = 0.5, size = 2) + 
        facet_grid(KMCluster_Under5 ~ Under5_R2) +
        theme_bw() + 
        ggtitle("K-means Clustering: Under5")

ggplot(PCA_Clustering, 
       aes(x = PC1_Under5,
           y = PC2_Under5,
           shape = Under5_R2,
           color = HCluster_Under5)) + 
        geom_point(alpha = 0.5, size = 2) + 
        facet_grid(HCluster_Under5 ~ Under5_R2) +
        theme_bw() + 
        ggtitle("Hierarchical Clustering: Under5")

nut3_nst_PCACluster <- nut3_nst %>%
        inner_join(PCA_Clustering, by = "Country_name")

nut3_nst_PCACluster[Neonatal_R2 == ">= 0.9" & KMCluster_Neonatal == 2]$Coef_neonatal
nut3_nst_PCACluster[Neonatal_R2 == ">= 0.9" & KMCluster_Neonatal == 2]$data                                                 

nut3_nst_PCACluster[KMCluster_Neonatal == 6]$Coef_neonatal
nut3_nst_PCACluster[KMCluster_Neonatal == 6]$data                                                 




####################### tSNE


# Running tSNE
library(Rtsne)
tSNE_fn <- function(dt, pp) {
        
        dt1 <- dt[, -1]
        
        # row labeling 
        rownames(dt1) <- dt$Country_name
        dt1 <- as.matrix(dt1)
        
        # tSNE
        set.seed(24)
        Rtsne(dt1,
              PCA = T, 
              perplexity = pp, 
              max_iter = 2000)
}

# Running tSNE: neonatal
tsne_neonatal5 <- tSNE_fn(Coef_neonatal, 5)
tsne_neonatal10 <- tSNE_fn(Coef_neonatal, 10)
tsne_neonatal25 <- tSNE_fn(Coef_neonatal, 25)

# Running tSNE: infant
tsne_infant5 <- tSNE_fn(Coef_infant, 5)
tsne_infant10 <- tSNE_fn(Coef_infant, 10)
tsne_infant25 <- tSNE_fn(Coef_infant, 25)

# Running tSNE: under5
tsne_under5 <- tSNE_fn(Coef_under5, 5)
tsne_under10 <- tSNE_fn(Coef_under5, 10)
tsne_under25 <- tSNE_fn(Coef_under5, 25)

# Iteration inspection
plot(tsne_neonatal25$itercosts, 
     type = "l",
     ylab = "Total K-L Divergence Cost",
     xlab = "Gradient Descent (50 Steps Each)",
     main = "Optimal Number of Iterations")

plot(tsne_neonatal25$costs, type = "l")


# Extracting tSNE coordinates 
tSNECoordTable_fn <- function(tsneObject, pp, group) {
        
        set.seed(113)
        data.frame(X = tsneObject$Y[, 1], 
                   Y = tsneObject$Y[, 2],
                   Perplexity = pp,
                   Group = group) %>%
                mutate(Perplexity = factor(Perplexity)) 
        
}

tSNE_Coordinates <- rbind(tSNECoordTable_fn(tsne_neonatal5, 5, "Neonatal"),
                          tSNECoordTable_fn(tsne_neonatal10, 10, "Neonatal"),
                          tSNECoordTable_fn(tsne_neonatal25, 25, "Neonatal"),
                          tSNECoordTable_fn(tsne_infant5, 5, "Infant"),
                          tSNECoordTable_fn(tsne_infant10, 10, "Infant"),
                          tSNECoordTable_fn(tsne_infant25, 25, "Infant"),
                          tSNECoordTable_fn(tsne_under5, 5, "Under5"),
                          tSNECoordTable_fn(tsne_under10, 10, "Under5"),
                          tSNECoordTable_fn(tsne_under25, 25, "Under5")) %>%
        mutate(Group = factor(Group, 
                              levels = c("Neonatal", "Infant", "Under5")))

# Visualizing tSNE results with perplexity 5, 10, and 25
tSNE_Perplexity_plot <- ggplot(tSNE_Coordinates,
                               aes(x = X, 
                                   y = Y,
                                   color = Perplexity,
                                   shape = Group)) + 
        geom_point(alpha = 0.5, size = 1.5) + 
        facet_grid(Perplexity ~ Group) + 
        theme_bw() + 
        ggtitle("tSNE with Perplexity 5, 10, and 25")


# clustering with tSNE coordinates where perplexity = 10
set.seed(1987)
tSNE_Coordinates_pp10 <- subset(tSNE_Coordinates, 
                                Perplexity == 10) %>% 
        
        mutate(Country_name = nut3_lmR2$Country_name) %>%
        
        group_by(Group, Perplexity) %>%
        
        # nesting
        nest() %>%
        
        # Labeling row names w country names
        mutate(data = map(data, 
                          ~ column_to_rownames(.x, var = "Country_name")),
               
               # Heatmap
               HeatMap = map(data, 
                             ~ pheatmap(.x,
                                        main = "Immunization and Mortality Rate")),
               
               # Scree plot to determine k
               KMScree = map(data, 
                             ~ Scree_fn(.x, "Scree Plot", 3))) %>%
        
        # clustering k = 3 or 4
        mutate(data = map(data, 
                          ~ mutate(.x, kmCluster3 = factor(kmeans(.x, centers = 3, 
                                                                  nstart = 25)$cluster))),
               data = map(data, 
                          ~ mutate(.x, kmCluster4 = factor(kmeans(.x, centers = 4, 
                                                                  nstart = 25)$cluster))),
               data = map(data, 
                          ~ mutate(.x, hCluster3 = factor(cutree(hclust(dist(.x), 
                                                                        method = "average"),
                                                                 k = 3)))),
               data = map(data, 
                          ~ mutate(.x, hCluster4 = factor(cutree(hclust(dist(.x), 
                                                                        method = "average"),
                                                                 k = 4))))) 