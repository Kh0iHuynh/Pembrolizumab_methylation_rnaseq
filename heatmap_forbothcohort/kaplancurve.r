# Install and load required packages
library(RColorBrewer)
library(survival)
library(survminer)
library(tidyverse)


pdf("kaplan_meier.pdf", width = 7, height = 7)
data <- read.table("kaplandata.txt", header = TRUE, sep = " ")  # Adjust path and separator if needed

#------------
# combined data curve
#-----------

# Step 3: Create the survival object
surv_object <- Surv(data$Days)

# Step 4: Kaplan-Meier survival model by Cohort
km_model <- survfit(surv_object ~ Cohort, data = data)

# Step 5: Plot the Kaplan-Meier curve
ggsurvplot(km_model, 
           data = data, 
           pval = TRUE,        # Add p-value for log-rank test
           risk.table = TRUE,  # Show the risk table
           legend.title = "Cohort",  # Title for the legend
           legend.labs = unique(data$Cohort),  # Labels for the legend (based on your Cohort column)
           palette = c("#E69F00", "#56B4E9", "#009E73"),  # Customize colors if you have multiple cohorts
           title = "Kaplan-Meier Survival Curve by Cohort",  # Title of the plot
           xlab = "Time (Days)",  # Label for the x-axis
           ylab = "Survival Probability")  # Label for the y-axis


#----------
# separate by cohort
#-----------
data <- read.table("kaplandata.txt", header = TRUE, sep = " ")  # Adjust path and separator if needed
data1 = data %>% filter(Cohort == "A")
# Step 3: Create the survival object
surv_object <- Surv(data1$Days)

# Step 4: Kaplan-Meier survival model by Cohort
km_model <- survfit(surv_object ~ Overall_Survival, data = data1)

# Step 5: Plot the Kaplan-Meier curve
ggsurvplot(km_model,
           data = data1,
           pval = TRUE,        # Add p-value for log-rank test
           risk.table = TRUE,  # Show the risk table
           legend.title = "Overall Survival",  # Title for the legend
           legend.labs = unique(data1$Overall_Survival),  # Labels for the legend (based on your Cohort column)
           palette = c("#E69F00", "#56B4E9", "#009E73"),  # Customize colors if you have multiple cohorts
           title = "Kaplan-Meier Survival Curve by Over Survival for Cohort A",  # Title of the plot
           xlab = "Time (Days)",  # Label for the x-axis
           ylab = "Survival Probability")  # Label for the y-axis



data1 = data %>% filter(Cohort == "B")
# Step 3: Create the survival object
surv_object <- Surv(data1$Days)

# Step 4: Kaplan-Meier survival model by Cohort
km_model <- survfit(surv_object ~ Overall_Survival, data = data1)

# Step 5: Plot the Kaplan-Meier curve
ggsurvplot(km_model,
           data = data1,
           pval = TRUE,        # Add p-value for log-rank test
           risk.table = TRUE,  # Show the risk table
           legend.title = "Overall Survival",  # Title for the legend
           legend.labs = unique(data1$Overall_Survival),  # Labels for the legend (based on your Cohort column)
           palette = c("#E69F00", "#56B4E9", "#009E73"),  # Customize colors if you have multiple cohorts
           title = "Kaplan-Meier Survival Curve by Overall Survival for Cohort B",  # Title of the plot
           xlab = "Time (Days)",  # Label for the x-axis
           ylab = "Survival Probability")  # Label for the y-axis


#####################
data <- read.table("kaplandata.txt", header = TRUE, sep = " ")  # Adjust path and separator if needed
data$Cohort_Overall_Survival <- paste0(data$Cohort,"_", data$Overall_Survival)

# Step 3: Create the survival object
surv_object <- Surv(data$Days)

# Step 4: Kaplan-Meier survival model by the new combined factor 'Cohort_Overall_Survival'
km_model <- survfit(surv_object ~ data$Cohort_Overall_Survival, data = data)

print("done make model;")
# Step 5: Plot the Kaplan-Meier curve
print("plotting")
ggsurvplot(km_model,
           data = data,
           pval = TRUE,        # Add p-value for log-rank test
           risk.table = TRUE,  # Show the risk table
           legend.title = "Cohort & Overall Survival",  # Title for the legend
           legend.labs = unique(data$Cohort_Overall_Survival),  # Labels for the legend (based on the new combined factor)
           palette = "Set1",  # Use a predefined color palette (you can customize if needed)
           title = "Kaplan-Meier Survival Curve by Cohort and by Overall Survival",  # Title of the plot
           xlab = "Time (Days)",  # Label for the x-axis
           ylab = "Survival Probability")  # Label for the y-axis

dev.off()
