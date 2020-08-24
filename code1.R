library(ggplot2)
library(broom)
library(tidyverse)
library(data.table)

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



##################################### Constructing Models ##################################### 

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

        
Coef_neonatal <- CleanData_fn("Est_neonatal", "Neonatal")
Coef_infant <- CleanData_fn("Est_infant", "Infant")
Coef_under5 <- CleanData_fn("Est_under5", "Under5")



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


# scree plots
fviz_eig(pca_neonatal)
fviz_eig(pca_infant)
fviz_eig(pca_under5)

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

# K = 4? 

plot(hc_neonatal)
plot(hc_infant)
plot(hc_under5)
