## Title:  Survival Analysis and Treatment Outcome in Leukaemia Patients
## Prepared by:  Patrick Mensah 
## Dataset : Gehan dataset. 


### Setting working directory 
setwd("C:/Users/patme/OneDrive/Desktop/Data_Science/R_programming/data")

## Installing pacman and loading all the packages involved########## 

install.packages(c("survival", "tidyverse", "MASS"))
library(survival)
library(tidyverse)
library(MASS)
library(ggplot2)
?gehan





# Exploring the dataset#########################

# Summary statistics for remission times (overall and by treatment group)
summary(gehan)
by(gehan$time, gehan$treat, summary)

# Count the number of patients in each treatment group
table(gehan$treat)

## Boxplot 
### to compare the distribution of remission times between the treatment groups visually.
# Boxplot comparing remission times by treatment group

ggplot(gehan, aes(x = treat, y = time, fill = treat)) +
  geom_boxplot() +
  labs(x = "Treatment", y = "Remission Time (Weeks)",
       title = "Remission Times of Leukaemia Patients",
       caption = paste("Prepared by Patrick Mensah -", format(Sys.Date(), "%Y-%m-%d"))) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 1, margin = margin(t = 20, unit = "pt")))


# comparing the remission times by the censoring status. 
# The censoring status indicates whether the patient's remission time is observed (event)
# or censored (withdrawn from the trial without an event)
# Create a boxplot comparing remission times by censoring status
ggplot(gehan, aes(x = cens, y = time, fill = factor(cens))) +
  geom_boxplot() +
  labs(x = "Censoring Status", y = "Remission Time (Weeks)",
       title = "Remission Times of Leukaemia Patients",
       caption = paste("Prepared by Patrick Mensah -", format(Sys.Date(), "%Y-%m-%d"))) +
  scale_fill_discrete(name = "Censoring", labels = c("Event", "Censored")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 1, margin = margin(t = 20, unit = "pt")))


## boxplot to compare the remission times of leukemia patients based on the pair labels.
## This will help us understand how the remission times vary within each pair. 
# Create a boxplot comparing remission times by pair labels
ggplot(gehan, aes(x = factor(pair), y = time, fill = factor(pair))) +
  geom_boxplot() +
  labs(x = "Pair Labels", y = "Remission Time (Weeks)",
       title = "Remission Times of Leukaemia Patients",
       caption = paste("Prepared by Patrick Mensah -", format(Sys.Date(), "%Y-%m-%d"))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for better readability
        plot.title = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 1, margin = margin(t = 20, unit = "pt")))




## Censoring plots################### 
## censoring plot can be helpful in visualizing the censoring status of the data.
# Censoring plot
ggplot(gehan, aes(x = time, y = cens, color = factor(cens))) +
  geom_point(size = 3) +
  scale_color_manual(name = "Censoring", values = c("blue", "red"),
                     labels = c("Event", "Censored")) +
  labs(x = "Time (Weeks)", y = "Censoring Status",
       title = "Censoring Plot (Blue: Event, Red: Censoring)",
       caption = paste("Prepared by Patrick Mensah -", format(Sys.Date(), "%Y-%m-%d"))) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 1, margin = margin(t = 20, unit = "pt")),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_continuous(breaks = seq(0, max(gehan$time), by = 1))



## Survival Plot ######
## Survival Plot (Kaplan-Meier Plot)
## A Kaplan-Meier plot is commonly used in survival analysis to visualize the survival 
## probabilities over time. It shows the estimated survival function for each treatment group.
# Kaplan-Meier survival plot by treatment group
install.packages("survminer")
library(survival)
library(survminer)

# Create the Kaplan-Meier survival plot using ggsurvplot with a title
ggsurvplot(
  survfit(Surv(time, cens) ~ treat, data = gehan),
  data = gehan,
  pval = TRUE,   # Display p-value comparing the survival curves
  conf.int = TRUE,  # Display confidence intervals
  risk.table = TRUE, # Display a table with the number at risk
  xlim = c(0, 100),  # Set the x-axis limit for better visualization
  xlab = "Time (Weeks)", ylab = "Survival Probability",
  surv.plot.height = 0.7,  # Adjust the height of the survival plot
  surv.plot.width = 0.8,   # Adjust the width of the survival plot
  tables.theme = theme_cleantable(),  # Use a cleaner theme for the tables
  tables.title = "Number at Risk",  # Set a title for the risk table
  title = "Kaplan-Meier Survival Plot for Leukemia Patients",  # Set the title for the plot
  ggtheme = theme_minimal() + 
    theme(plot.title = element_text(size = 18, face = "bold", color = "darkblue", hjust = 0.5))
)


## Hazard plots ######
## Hazard plots, also known as hazard rate plots or hazard functions, 
## are commonly used in survival analysis to visualize the instantaneous failure rate 
## (hazard) over time. The hazard represents the probability of an event occurring at a given time,
## given that the subject has survived up to that time.



# Fit the survival model and calculate the Nelson-Aalen cumulative hazard
surv_fit <- survfit(Surv(time, cens) ~ treat, data = gehan)
hazard_est <- cumsum(surv_fit$n.event / surv_fit$n.risk)

# Create a visually appealing hazard plot
plot(surv_fit$time, hazard_est, type = "s", lwd = 2, col = c("blue", "red"),
     xlab = "Time (Weeks)", ylab = "Cumulative Hazard",
     main = "Nelson-Aalen Cumulative Hazard Plot for Leukemia Patients",
     cex.main = 1.5, cex.lab = 1.2, cex.axis = 1.2)

legend("topright", legend = c("Control", "6-MP"), col = c("blue", "red"),
       lty = 1, cex = 1, bty = "n")

mtext("Prepared by Patrick Mensah - July 2023",
      side = 1, line = 5, cex = 0.8, col = "gray")



#### ANALYSIS #####

## The log rank test####
## The log-rank test is a non-parametric test used to compare the survival curves of two or more groups. 
## In this case, we will use it to compare the survival curves of the control and 6-MP treatment groups in the Gehan dataset.

# Fit the survival model and calculate the Nelson-Aalen cumulative hazard
surv_fit <- survfit(Surv(time, cens) ~ treat, data = gehan)

# Perform the log-rank test
logrank_test <- survdiff(Surv(time, cens) ~ treat, data = gehan)
# Print the results
print(logrank_test)

## Since the p-value is less than the typical significance level of 0.05,
## we reject the null hypothesis and conclude that there is a statistically significant
## difference in survival between the control and 6-MP treatment groups in the Gehan dataset



## Cox Proportional Hazards Model 
## The Cox proportional hazards model is a semi-parametric survival regression model that 
## allows us to assess the effect of one or more predictor variables on the hazard of an event (e.g., relapse)
##while assuming that the hazard ratios are constant over time


# Fit the Cox Proportional Hazards Model
cox_model <- coxph(Surv(time, cens) ~ treat, data = gehan)
# Print the summary of the model
summary(cox_model)

## The significant p-values in all three tests (LR test, Wald test, and Score test) indicate 
##that the treatment group variable has a strong and statistically significant effect on the hazard of relapse.
##The hazard ratio of 4.817 means that the hazard of relapse is significantly higher in the 6-MP treatment
##group compared to the control group.

##Overall, the Cox proportional hazards model confirms that the treatment group has a significant impact
## on the survival outcome in the Gehan dataset, and patients in the 6-MP treatment group are at a higher 
## risk of relapse compared to those in the control group


## Censoring Analysis #####
## Censoring analysis is an important step in survival analysis to understand the pattern 
##of censoring and its potential impact on the study results. In the Gehan dataset, the "cens" variable indicates
## whether an event (relapse) was observed or if the data is right-censored (0 = event observed, 1 = right-censored). 



# Censoring analysis
censoring_summary <- table(gehan$cens)
censoring_percentage <- prop.table(censoring_summary) * 100

# Print the summary
cat("Censoring Summary:\n")
cat("===================\n")
cat("Observed Events (Censored = 0):", censoring_summary[1], "\n")
cat("Right-Censored Events (Censored = 1):", censoring_summary[2], "\n\n")

cat("Censoring Percentage:\n")
cat("=====================\n")
cat("Observed Events (%):", sprintf("%.2f%%", censoring_percentage[1]), "\n")
cat("Right-Censored Events (%):", sprintf("%.2f%%", censoring_percentage[2]), "\n")

## The high percentage of right-censored events (71.43%) suggests that a substantial portion of 
## patients in the Gehan dataset did not experience the event (relapse) during the study and were 
## censored at the end of the follow-up period. This information is crucial to consider when 
## interpreting the results of survival analyses and drawing conclusions about the impact of treatment
## or other covariates on survival outcomes.


## Time to event predictions
## To perform time-to-event prediction using a survival model, we can use the fitted Cox proportional 
## hazards model to estimate the survival probabilities for new leukemia patients based on their characteristics.
## This allows us to predict the time it would take for these new patients to experience an event (e.g., relapse)



# Fit the Cox Proportional Hazards Model and specify variable names explicitly
cox_model <- coxph(Surv(time, cens) ~ treat, data = gehan)

# Create a new data frame with predictor values for the new patients
# Ensure that variable names match those used in the gehan dataset
new_patients <- data.frame(time = c(10, 20, 30), cens = c(0, 1, 0), treat = c("control", "6-MP", "control"))

# Predict the survival probabilities for the new patients
surv_probs <- predict(cox_model, newdata = new_patients, type = "survival")

# Print the predicted survival probabilities
print(surv_probs)

## Output ##### 
##The output you provided from the time-to-event prediction represents the estimated survival probabilities for the new patients at different time points (10, 20, and 30 weeks) based on their treatment group and the fitted Cox proportional hazards model.

##Let's interpret the results:

##For patient 1 (treatment group: control), the predicted survival probability at 10 weeks is approximately 0.347, at 20 weeks is approximately 0.633, and at 30 weeks is approximately 0.020.

##For patient 2 (treatment group: 6-MP), the predicted survival probability at 10 weeks is approximately 0.633, at 20 weeks is approximately 0.633, and at 30 weeks is approximately 0.020.

##For patient 3 (treatment group: control), the predicted survival probability at 10 weeks is approximately 0.347, at 20 weeks is approximately 0.633, and at 30 weeks is approximately 0.020.

##Observations:

###The predicted survival probabilities for patients 1 and 3 are the same at all time points. This is because they both belong to the "control" treatment group, and the Cox proportional hazards model estimated a constant hazard ratio for this group.

##Patient 2 (treatment group: 6-MP) has a higher predicted survival probability at all time points compared to patients 1 and 3. This suggests that, based on the model, patients in the "6-MP" treatment group have a higher chance of survival over time.

##The survival probabilities decrease over time, which is expected in survival analysis, as the likelihood of experiencing the event (relapse) increases with time.





