# load packages
library(JM)
library(nlme)
library(lme4)
library(tidyverse)
library(gridExtra)
library(survival)
library(survminer)
library(gtsummary)
library(forecast)
library(gt)
library(labelled)
library(here)

data(pbc2.id, package = "JM")
pbc2.id$status2 <- as.numeric(pbc2.id$status != "alive" || pbc2.id$status != "transplanted")

# histologic to factor
pbc2$histologic <- factor(pbc2$histologic)
pbc2.id$histologic <- factor(pbc2.id$histologic)

### Summary Statistics
pbc2_summary <- pbc2.id %>%
  select(drug, sex, status2, years, age, ascites, hepatomegaly, spiders, edema, serBilir, 
         albumin, alkaline, SGOT, platelets, prothrombin, histologic)

uv_cox_fit %>%
  add_global_p() %>%
  bold_labels()

pbc2_summary %>% tbl_summary() %>%
  bold_labels() %>% 
  italicize_levels()

for (i in 1:length(pbc2_summary$status2)) {
  if (pbc2_summary$status2[i] == 0) {
    pbc2_summary$status2[i] <- "alive"
  } else{
    pbc2_summary$status2[i] <- "dead"
  }
}

pbc2_summary %>% tbl_summary() %>%
  bold_labels() %>% 
  italicize_levels()

ggplot(pbc2_summary, aes(x=serBilir)) + 
  geom_histogram(fill = "orange", aes(y=..density..), color = 'orange') + 
  xlab("Baseline Biilrubin") + 
  ylab("Count") + 
  ggtitle("Histogram of Baseline Bilirubin")

### Data Exploration
# Histogram for continuous variables
ggplot(data = pbc2.id, mapping = aes(x = age)) +
  geom_histogram()

g1 <- ggplot(data = pbc2, mapping = aes(x = serBilir)) +
  geom_histogram()

serBilir_lambda <- BoxCox.lambda( pbc2$serBilir )
transform_serBilir <- BoxCox(pbc2$serBilir, serBilir_lambda)

g2 <- ggplot(data = pbc2, mapping = aes(x = transform_serBilir)) +
  geom_histogram()
grid.arrange(g1, g2, ncol = 2)

ggplot(data = pbc2, mapping = aes(x = albumin)) +
  geom_histogram()

g5 <- ggplot(data = pbc2, mapping = aes(x = alkaline)) +
  geom_histogram()
g6 <- ggplot(data = pbc2, mapping = aes(x = log(alkaline))) +
  geom_histogram()
grid.arrange(g5, g6, ncol = 2)

g7 <- ggplot(data = pbc2, mapping = aes(x = SGOT)) +
  geom_histogram()
g8 <- ggplot(data = pbc2, mapping = aes(x = log(SGOT))) +
  geom_histogram()
grid.arrange(g7, g8, ncol = 2)

g9 <- ggplot(data = pbc2, mapping = aes(x = platelets)) +
  geom_histogram()
g10 <- ggplot(data = pbc2, mapping = aes(x = log(platelets))) +
  geom_histogram()
grid.arrange(g9, g10, ncol = 2)

g11 <- ggplot(data = pbc2, mapping = aes(x = prothrombin)) +
  geom_histogram()

g12 <- ggplot(mapping = aes(x = log(platelets))) +
  geom_histogram()

grid.arrange(g11, g12, ncol = 2)

### Survival Analysis
# Kaplan-Meier method
surv_object <- Surv(time = pbc2.id$years, event = pbc2.id$status2)

fit0 <- survfit(formula = surv_object ~ 1, data = pbc2.id)
summary(fit0)
ggsurvplot(fit0, data = pbc2.id, conf.int = TRUE, risk.table = TRUE, 
           xlab = "Time in years", surv.median.line = "hv", xlim = c(0, 14), break.x.by = 2, 
           title = "Overall survival probability")

#sex
fit1 <- survfit(formula = surv_object ~ sex, data = pbc2.id)
summary(fit1)
ggsurvplot(fit1, data = pbc2.id, conf.int = TRUE, risk.table = TRUE, pval = TRUE, 
           xlab = "Time in years", surv.median.line = "hv", risk.table.height = 0.3, xlim = c(0, 14), break.x.by = 2, 
           title = "Sex")

fit2 <- survfit(formula = surv_object ~ drug, data = pbc2.id)
summary(fit2)
ggsurvplot(fit2, data = pbc2.id, conf.int = TRUE, risk.table = TRUE, pval = TRUE, 
           xlab = "Time in years", surv.median.line = "hv", risk.table.height = 0.3, xlim = c(0, 14), break.x.by = 2, 
           title = "Treatment effect")

fit3 <- survfit(formula = surv_object ~ ascites, data = pbc2.id)
summary(fit3)
ggsurvplot(fit3, data = pbc2.id, conf.int = TRUE, risk.table = TRUE, pval = TRUE, 
           xlab = "Time in years", surv.median.line = "hv", risk.table.height = 0.3, xlim = c(0, 14), break.x.by = 2, title = "Ascites")

fit4 <- survfit(formula = surv_object ~ hepatomegaly, data = pbc2.id)
summary(fit4)
ggsurvplot(fit4, data = pbc2.id, conf.int = TRUE, risk.table = TRUE, pval = TRUE, 
           xlab = "Time in years", surv.median.line = "hv", risk.table.height = 0.3, xlim = c(0, 14), break.x.by = 2, title = "Hepatomegaly")

fit5 <- survfit(formula = surv_object ~ spiders, data = pbc2.id)
summary(fit5)
ggsurvplot(fit5, data = pbc2.id, conf.int = TRUE, risk.table = TRUE, pval = TRUE, 
           xlab = "Time in years", surv.median.line = "hv", risk.table.height = 0.3, xlim = c(0, 14), break.x.by = 2, title = "Spiders")

fit6 <- survfit(formula = surv_object ~ edema, data = pbc2.id)
summary(fit6)
ggsurvplot(fit6, data = pbc2.id, conf.int = TRUE, risk.table = TRUE, pval = TRUE, 
           xlab = "Time in years", surv.median.line = "hv", risk.table.height = 0.33, xlim = c(0, 14), break.x.by = 2, title = "Edema")

fit7 <- survfit(formula = surv_object ~ histologic, data = pbc2.id)
summary(fit7)
ggsurvplot(fit7, data = pbc2.id, conf.int = TRUE, risk.table = TRUE, pval = TRUE, 
           xlab = "Time in years", surv.median.line = "hv", risk.table.height = 0.38, xlim = c(0, 14), break.x.by = 2, 
           title = "Histologic")

pbc2.id <- pbc2.id %>% 
  mutate(age_group = ifelse(age >= 60, "old", 
                            ifelse(age <= 40, "young", "middle-aged")))
pbc2.id$age_group <- factor(pbc2.id$age_group)

fit8 <- survfit(formula = surv_object ~ age_group, data = pbc2.id)
summary(fit8)
ggsurvplot(fit8, data = pbc2.id, conf.int = TRUE, risk.table = TRUE, pval = TRUE, 
           xlab = "Time in years", surv.median.line = "hv", risk.table.height = 0.33, xlim = c(0, 14), break.x.by = 2)

pbc2.id <- pbc2.id %>% 
  mutate(bilirubin = ifelse(serBilir >= 1.2, "high > 1.2 ", "normal <= 1.2"))

fit9 <- survfit(formula = surv_object ~ bilirubin, data = pbc2.id)
summary(fit8)
ggsurvplot(fit9, data = pbc2.id, conf.int = TRUE, risk.table = TRUE, pval = TRUE, 
           xlab = "Time in years", surv.median.line = "hv", risk.table.height = 0.33, xlim = c(0, 14), break.x.by = 2, 
           title = "Baseline Bilirubin")

### Variable Selection for Cox Regression Model
#summary(pbc2.id) # platelets and serChol missing

# complete case analysis - remove missingness for platelets
pbc2.id_complete_case <- pbc2.id[!is.na(pbc2.id$platelets),]

# backward selection without drug 
cox_fit1 <- coxph(Surv(years, status2) ~ platelets + serBilir + albumin + alkaline + SGOT + prothrombin +
                    age + sex + ascites + hepatomegaly + spiders + edema + histologic, data = pbc2.id_complete_case)
summary(cox_fit1)
drop1(cox_fit1, test = "Chisq") # drop alkaline with p-value 0.9245082 

cox_fit2 <- coxph(Surv(years, status2) ~ platelets + serBilir + albumin + SGOT + prothrombin +
                    age + sex + ascites + hepatomegaly + spiders + edema + histologic, data = pbc2.id_complete_case)
drop1(cox_fit2, test = "Chisq") # drop spiders with p-value 0.7176748   

cox_fit3 <- coxph(Surv(years, status2) ~ platelets + serBilir + albumin + SGOT + prothrombin +
                    age + sex + ascites + hepatomegaly + edema + histologic, data = pbc2.id_complete_case)
drop1(cox_fit3, test = "Chisq") # drop platelets with p-value 0.4968202 

cox_fit4 <- coxph(Surv(years, status2) ~ serBilir + albumin + SGOT + prothrombin +
                    age + sex + ascites + hepatomegaly + edema + histologic, data = pbc2.id_complete_case)
drop1(cox_fit4, test = "Chisq") # drop ascites with p-value of 0.5098511  

cox_fit5 <- coxph(Surv(years, status2) ~ serBilir + albumin + SGOT + prothrombin +
                    age + sex + hepatomegaly + edema + histologic, data = pbc2.id_complete_case)
drop1(cox_fit5, test = "Chisq") # drop hepatomegaly with p-value of 0.1757678

cox_fit6 <- coxph(Surv(years, status2) ~ serBilir + albumin + SGOT + prothrombin +
                    age + sex + histologic + edema, data = pbc2.id_complete_case)
drop1(cox_fit6, test = "Chisq") # drop sex with p-value of 0.0618493

# forward selection
cox_fit1 <- coxph(Surv(years, status2) ~ serBilir, data = pbc2.id_complete_case)
add1(cox_fit1, ~ drug + serBilir + albumin + alkaline + SGOT + prothrombin + age + sex + ascites +
       hepatomegaly + spiders + edema + histologic, test = "Chisq") # add histologic with p-value of 4.102e-10

cox_fit2 <- coxph(Surv(years, status2) ~ serBilir + histologic, data = pbc2.id_complete_case)
add1(cox_fit2, ~ drug + serBilir + albumin + alkaline + SGOT + prothrombin + age + sex + ascites +
       hepatomegaly + spiders + edema + histologic, test = "Chisq") # add age with p-value of 4.633e-06

cox_fit3 <- coxph(Surv(years, status2) ~ serBilir + histologic + age, data = pbc2.id_complete_case)
add1(cox_fit3, ~ drug + serBilir + albumin + alkaline + SGOT + prothrombin + age + sex + ascites +
       hepatomegaly + spiders + edema + histologic, test = "Chisq") # add prothrombin with p-value of 6.480e-05

cox_fit4 <- coxph(Surv(years, status2) ~ serBilir + histologic + age + prothrombin, data = pbc2.id_complete_case)
add1(cox_fit4, ~ drug + serBilir + albumin + alkaline + SGOT + prothrombin + age + sex + ascites +
       hepatomegaly + spiders + edema + histologic, test = "Chisq") # add albumin with p-value of 0.0002209

cox_fit5 <- coxph(Surv(years, status2) ~ serBilir + histologic + age + prothrombin + albumin, data = pbc2.id_complete_case)
add1(cox_fit5, ~ drug + serBilir + albumin + alkaline + SGOT + prothrombin + age + sex + ascites +
       hepatomegaly + spiders + edema + histologic, test = "Chisq") # add edema with p-value of 0.01472

cox_fit6 <- coxph(Surv(years, status2) ~ serBilir + histologic + age + prothrombin + albumin + edema, data = pbc2.id_complete_case)
add1(cox_fit6, ~ drug + serBilir + albumin + alkaline + SGOT + prothrombin + age + sex + ascites +
       hepatomegaly + spiders + edema + histologic, test = "Chisq") # add SGOT with p-value of 0.02093

cox_fit7 <- coxph(Surv(years, status2) ~ serBilir + histologic + age + prothrombin + albumin + edema + SGOT, data = pbc2.id_complete_case)
add1(cox_fit7, ~ drug + serBilir + albumin + alkaline + SGOT + prothrombin + age + sex + ascites +
       hepatomegaly + spiders + edema + histologic, test = "Chisq")

### Multivariable Cox model
# Cox proportional hazards model
coxFit <- coxph(Surv(years, status2) ~ serBilir + albumin + age + edema + histologic + SGOT + prothrombin, data = pbc2.id)

# summary(coxFit)
ggforest(coxFit, data = pbc2.id, fontsize = 0.55)

coxFit %>%
  tbl_regression(exponentiate = TRUE) %>%
  add_global_p() %>%
  bold_labels() %>%
  italicize_levels()  %>%
  as_kable(caption = "Multivariable Cox Model")

### Test the Proportional Hazards Assumption of a Cox Regression
cox.zph(coxFit)

### Data for time-dependent Cox
result <- create_time(pbc2$year, pbc2$years, pbc2$status2)
start <- as.vector(unlist(result[[1]]))
stop <- as.vector(unlist(result[[2]])) 
event <- as.vector(unlist(result[[3]])) 

pbc2_visit_time <- pbc2
pbc2_visit_time$start <- start
pbc2_visit_time$end <- stop
pbc2_visit_time$event <- event

# data modification
pbc2_visit_time_merge <- merge(pbc2.id, pbc2_visit_time, by = "id", sort = FALSE)

# obstain baseline values
pbc2_visit_time <- pbc2_visit_time %>%
  mutate(
    albumin_baseline = pbc2_visit_time_merge[,"albumin.x"],
    serBilir_baseline = pbc2_visit_time_merge[,"serBilir.x"],
    serChol_baseline = pbc2_visit_time_merge[,"serChol.x"],
    alkaline_baseline = pbc2_visit_time_merge[,"alkaline.x"],
    SGOT_baseline = pbc2_visit_time_merge[,"SGOT.x"],
    platelets_baseline = pbc2_visit_time_merge[,"platelets.x"],
    prothrombin_baseline = pbc2_visit_time_merge[,"prothrombin.x"],
    ascites_baseline = pbc2_visit_time_merge[,"ascites.x"],
    hepatomegaly_baseline = pbc2_visit_time_merge[,"hepatomegaly.x"],
    spiders_baseline = pbc2_visit_time_merge[,"spiders.x"],
    edema_baseline = pbc2_visit_time_merge[,"edema.x"],
    histologic_baseline = pbc2_visit_time_merge[,"histologic.x"]
  )

### Univariable time-dependent Cox
uv_td_cox_fit <- tbl_uvregression(
  pbc2_visit_time %>% select(c(-id, -years, -status, -year, -status2, -start, -end, -event, -ascites, -hepatomegaly, -status_factor,
                -serChol_baseline, -spiders, -edema, -serBilir_baseline, -serChol, -albumin, -alkaline, -SGOT, -platelets, -prothrombin, -histologic)),
  method = coxph,
  y = Surv(pbc2_visit_time$start, pbc2_visit_time$end, pbc2_visit_time$event),
  exponentiate = TRUE
)

uv_td_cox_fit %>%
  add_global_p() %>%
  bold_labels()  %>%
  as_kable(caption = "Univariable time-dependent Cox Model")

### Complete Case Analysis for extended Cox
# summary(pbc2_visit_time) # platelets_baseline and serChol_baseline missing
# complete case analysis - remove missingness for platelets
pbc2_visit_time_complete_case <- pbc2_visit_time[!is.na(pbc2_visit_time$platelets_baseline), ]

result1 <- create_time(pbc2_visit_time_complete_case$year, pbc2_visit_time_complete_case$years, pbc2_visit_time_complete_case$status2)
start1 <- as.vector(unlist(result1[[1]]))
stop1 <- as.vector(unlist(result1[[2]])) 
event1 <- as.vector(unlist(result1[[3]])) 

pbc2_visit_time_complete_case$start1 <- start1
pbc2_visit_time_complete_case$end1 <- stop1
pbc2_visit_time_complete_case$event1 <- event1

### Variable Selection for extended Cox model
# backward selection
td_cox1 <- coxph(Surv(start1, stop1, event1) ~ platelets_baseline + serBilir + albumin_baseline +
                   alkaline_baseline + SGOT_baseline + prothrombin_baseline + ascites_baseline + age + sex + 
                   hepatomegaly_baseline + spiders_baseline + edema_baseline + histologic_baseline, data = pbc2_visit_time_complete_case)
drop1(td_cox1, test = "Chisq") # drop ascites_baseline with p-value of 0.828480 

cox_fit2 <- coxph(Surv(start1, stop1, event1) ~ platelets_baseline + serBilir + albumin_baseline +
                    alkaline_baseline + SGOT_baseline + prothrombin_baseline + age + sex + 
                    hepatomegaly_baseline + spiders_baseline + edema_baseline + histologic_baseline, data = pbc2_visit_time_complete_case)
drop1(cox_fit2, test = "Chisq") # drop spiders_baseline with p-value of 0.696068

cox_fit3 <- coxph(Surv(start1, stop1, event1) ~ platelets_baseline + serBilir + albumin_baseline +
                    alkaline_baseline + SGOT_baseline + prothrombin_baseline + age + sex + 
                    hepatomegaly_baseline + edema_baseline + histologic_baseline, data = pbc2_visit_time_complete_case)
drop1(cox_fit3, test = "Chisq") # drop alkaline_baseline with p-value of 0.604864  

cox_fit4 <- coxph(Surv(start1, stop1, event1) ~ platelets_baseline + serBilir + albumin_baseline +
                    SGOT_baseline + prothrombin_baseline + age + sex + 
                    hepatomegaly_baseline + edema_baseline + histologic_baseline, data = pbc2_visit_time_complete_case)
drop1(cox_fit4, test = "Chisq") # drop hepatomegaly_baseline with p-value of 0.577356  

cox_fit5 <- coxph(Surv(start1, stop1, event1) ~ platelets_baseline + serBilir + albumin_baseline +
                    SGOT_baseline + prothrombin_baseline + age + sex + 
                    edema_baseline + histologic_baseline, data = pbc2_visit_time_complete_case)
drop1(cox_fit5, test = "Chisq") # drop prothrombin_baseline with p-value 0.4048648 

cox_fit6 <- coxph(Surv(start1, stop1, event1) ~ platelets_baseline + serBilir + albumin_baseline +
                    SGOT_baseline + age + sex + edema_baseline + histologic_baseline, data = pbc2_visit_time_complete_case)
drop1(cox_fit6, test = "Chisq") # drop SGOT_baseline with p-value 0.2043608 

cox_fit7 <- coxph(Surv(start1, stop1, event1) ~ platelets_baseline + serBilir + albumin_baseline +
                    age + sex + edema_baseline + histologic_baseline, data = pbc2_visit_time_complete_case)
drop1(cox_fit7, test = "Chisq") # drop platelets_baseline with p-value 0.254760   

cox_fit8 <- coxph(Surv(start1, stop1, event1) ~ serBilir + albumin_baseline +
                    age + sex + edema_baseline + histologic_baseline, data = pbc2_visit_time_complete_case)
drop1(cox_fit8, test = "Chisq") # drop sex with p-value 0.1069981

# backward selection
td_cox1 <- coxph(Surv(start1, stop1, event1) ~ serBilir, data = pbc2_visit_time_complete_case)
add1(td_cox1, ~ drug + serBilir + albumin_baseline + SGOT_baseline + prothrombin_baseline + ascites_baseline + age + sex +
       platelets_baseline + edema_baseline + histologic_baseline + spiders_baseline + hepatomegaly_baseline + 
       alkaline_baseline, test = "Chisq") # add age < 2.2e-16 

td_cox2 <- coxph(Surv(start1, stop1, event1) ~ serBilir + age, data = pbc2_visit_time_complete_case)
add1(td_cox2, ~ drug + serBilir + albumin_baseline + SGOT_baseline + prothrombin_baseline + ascites_baseline + age + sex +
       platelets_baseline + edema_baseline + histologic_baseline + spiders_baseline + hepatomegaly_baseline + 
       alkaline_baseline, test = "Chisq") # add histologic_baseline 3.441e-09

td_cox3 <- coxph(Surv(start1, stop1, event1) ~ serBilir + histologic_baseline + age, data = pbc2_visit_time_complete_case)
add1(td_cox3, ~ drug + serBilir + albumin_baseline + SGOT_baseline + prothrombin_baseline + ascites_baseline + age + sex +
       platelets_baseline + edema_baseline + histologic_baseline + spiders_baseline + hepatomegaly_baseline + 
       alkaline_baseline, test = "Chisq") # add edema_baseline 9.46e-05

td_cox4 <- coxph(Surv(start1, stop1, event1) ~ serBilir + histologic_baseline + age + edema_baseline, data = pbc2_visit_time_complete_case)
add1(td_cox4, ~ drug + serBilir + albumin_baseline + SGOT_baseline + prothrombin_baseline + ascites_baseline + age + sex +
       platelets_baseline + edema_baseline + histologic_baseline + spiders_baseline + hepatomegaly_baseline + 
       alkaline_baseline, test = "Chisq") # add albumin_baseline 0.0592

td_cox5 <- coxph(Surv(start1, stop1, event1) ~ serBilir + histologic_baseline + age + edema_baseline + albumin_baseline, 
                 data = pbc2_visit_time_complete_case)
add1(td_cox5, ~ drug + serBilir + albumin_baseline + SGOT_baseline + prothrombin_baseline + ascites_baseline + age + sex +
       platelets_baseline + edema_baseline + histologic_baseline + spiders_baseline + hepatomegaly_baseline + 
       alkaline_baseline, test = "Chisq") 

### Multivariable time-dependent Cox model
td_cox_bilir <- coxph(Surv(start, stop, event) ~ serBilir + albumin_baseline + age + edema_baseline + histologic_baseline, 
                      data = pbc2_visit_time)

# summary(td_cox_bilir)
td_cox_bilir %>%
  tbl_regression(exponentiate = TRUE) %>%
  add_global_p() %>%
  bold_labels() %>%
  italicize_levels() %>%
  as_kable(caption = "Final Multivariable Time-dependent Cox Model")

### Test the Proportional Hazards Assumption of Extended Cox Model
cox.zph(td_cox_bilir)

### Longitudinal Analysis - Linear Mix-effects Model
### Description 
# plot
ggplot(pbc2, aes(x = year, y = serBilir, color = factor(id))) +
  geom_line() + geom_point() +
  theme_bw() +
  theme(legend.position = "none")

ggplot(data = pbc2, aes(x = year, y = serBilir, group = id, color = factor(status_factor))) + 
  geom_line() +
  facet_grid(. ~ status_factor) +
  ggtitle("Longitudinal trajectories of bilirubin")

ggplot(data = pbc2, aes(x = year, y = serBilir, group = id)) + 
  geom_line() +
  facet_grid(. ~ sex)

### Data for LME
# transform dataset
pbc2_merge <- merge(pbc2.id, pbc2, by = "id", sort = FALSE)

pbc2_transform <- pbc2 %>%
  mutate(
    serBilir = transform_serBilir,
    serChol = log(serChol),
    alkaline = log(alkaline), 
    SGOT = log(SGOT),
    platelets = log(platelets),
    prothrombin = log(prothrombin),
    
    serBilir_baseline = BoxCox(pbc2_merge[,"serBilir.x"], serBilir_lambda),
    serChol_baseline = log(pbc2_merge[,"serChol.x"]),
    albumin_baseline = pbc2_merge[,"albumin.x"],
    alkaline_baseline = log(pbc2_merge[,"alkaline.x"]),
    SGOT_baseline = log(pbc2_merge[,"SGOT.x"]),
    platelets_baseline = log(pbc2_merge[,"platelets.x"]),
    prothrombin_baseline = log(pbc2_merge[,"prothrombin.x"]),
    
    ascites_baseline = pbc2_merge[,"ascites.x"],
    hepatomegaly_baseline = pbc2_merge[,"hepatomegaly.x"],
    spiders_baseline = pbc2_merge[,"spiders.x"],
    edema_baseline = pbc2_merge[,"edema.x"],
    histologic_baseline = pbc2_merge[,"histologic.x"]
  )

### Data Imputation  
# remove rows with missing values for bilirubin - response variable 
for (i in 1:length(pbc2_transform$serBilir)) {
  if (is.na(pbc2_transform$serBilir[i])){
    pbc2_transform <- pbc2_transform[-c(i)]
  }
}

# impute transformed data
impute_missing_data <- function(patient, biomarker) {
  result <- c()
  result[1] <- biomarker[1]
  for (i in 2:(length(patient))) {
    if (is.na(biomarker[i])) {
      biomarker[i] <- biomarker[i-1]
    }
    else {}
    result[i] <- biomarker[i]
  }
  return(result)
}

# impute baseline data
platelets_baseline <- impute_missing_data(pbc2_transform$id, pbc2_transform$platelets_baseline)
pbc2_transform$platelets_baseline <- as.vector(unlist(platelets_baseline))

pbc2_transform_complete_case <- pbc2_transform
pbc2_transform_complete_case <- pbc2_transform_complete_case[!is.na(pbc2_transform_complete_case$platelets_baseline),]

### Univariable Analysis
pbc2_transform %>%
  select(c(-years, -status, -status2, -ascites, -hepatomegaly, -spiders, -edema, -albumin, -serChol, -serBilir_baseline, 
           serChol_baseline, -alkaline, -SGOT, -platelets, -prothrombin, -histologic, -status_factor, -serChol_baseline)) %>%
  tbl_uvregression(
    y = serBilir,
    method = lme4::lmer,
    formula = "{y} ~ {x} + (year | id)"
  ) %>%
  add_global_p() %>%
  bold_labels() %>%
  italicize_levels()

### Variable selection for Linear Mixed-Effect
# backward selection
lme_fit1 <- lmer(serBilir ~ year + ascites_baseline + sex + hepatomegaly_baseline + spiders_baseline + edema_baseline + 
                   albumin_baseline + alkaline_baseline + SGOT_baseline + prothrombin_baseline + 
                   histologic_baseline + platelets_baseline + (year | id), data = pbc2_transform) 
drop1(lme_fit1, test = "Chisq") # drop edema_baseline 0.8040889    

lme_fit2 <- lmer(serBilir ~ year + ascites_baseline + sex + hepatomegaly_baseline + spiders_baseline + 
                   albumin_baseline + alkaline_baseline + SGOT_baseline + prothrombin_baseline + 
                   histologic_baseline + platelets_baseline + (year | id), data = pbc2_transform) 
drop1(lme_fit2, test = "Chisq") # drop platelets_baseline 0.5444112  

lme_fit3 <- lmer(serBilir ~ year + ascites_baseline + sex + hepatomegaly_baseline + spiders_baseline + 
                   albumin_baseline + alkaline_baseline + SGOT_baseline + prothrombin_baseline + 
                   histologic_baseline + (year | id), data = pbc2_transform) 
drop1(lme_fit3, test = "Chisq") # drop histologic_baseline 0.3581601 

lme_fit4 <- lmer(serBilir ~ year + ascites_baseline + sex + hepatomegaly_baseline + spiders_baseline + 
                   albumin_baseline + alkaline_baseline + SGOT_baseline + prothrombin_baseline + (year | id), data = pbc2_transform) 
drop1(lme_fit4, test = "Chisq") 

### Multivariable Linear Mixed-Effects Model
final_lme_fit <- lmer(serBilir ~ year + ascites_baseline + sex + hepatomegaly_baseline + spiders_baseline + 
                        albumin_baseline + alkaline_baseline + SGOT_baseline + prothrombin_baseline + (year | id), data = pbc2_transform)

final_lme_fit %>%
  tbl_regression() %>%
  add_global_p(include = everything()) %>%
  bold_labels() %>%
  italicize_levels() %>%
  as_kable(caption = "Final Multivariable Linear Mix-Effects Model")

final_lme_fit1 <- lme(serBilir ~ year + ascites_baseline + sex + hepatomegaly_baseline + spiders_baseline + 
                        albumin_baseline + alkaline_baseline + SGOT_baseline + prothrombin_baseline, random = ~ year | id, data = pbc2_transform) 

# summary(final_lme_fit1)
summary(final_lme_fit)

### Diagnostic Plot
pbc2_transform <- pbc2_transform %>%
  mutate(
    residuals = summary(final_lme_fit)$residuals,
    fitted = summary(final_lme_fit1)$fitted[,1]
  )

ggplot(data = pbc2_transform, mapping = aes(sample = residuals)) +
  stat_qq() +
  stat_qq_line() +
  ggtitle("Residuals Q-Q")

ggplot(data = pbc2_transform, mapping = aes(x = residuals)) +
  geom_density() +
  ggtitle("Residuals")

ggplot(data = pbc2_transform, mapping = aes(y = residuals, x = fitted)) +
  geom_point() +
  ggtitle("Residuals")

### Joint Model with Serum Bilirubin as response variable
lme_fit_bilir <- lme(serBilir ~ year + ascites_baseline + sex + hepatomegaly_baseline + spiders_baseline + 
                       albumin_baseline + alkaline_baseline + SGOT_baseline + prothrombin_baseline, random = ~ year | id, 
                     data = pbc2_transform)

###  Variable Selection for Cox part of JM
# drop hepatomegaly
cox_fit_bilir1 <- coxph(Surv(years, status2) ~ albumin + alkaline + SGOT + prothrombin + age + sex + ascites +
                          hepatomegaly + spiders + edema + histologic, data = pbc2.id, x = TRUE, model = TRUE)
joint_fit_bilir1 <- jointModel(lme_fit_bilir, cox_fit_bilir1, timeVar = "year", method = "piecewise-PH-GH", iter.EM = 10)
summary(joint_fit_bilir1)

# drop SGOT
cox_fit_bilir2 <- coxph(Surv(years, status2) ~ albumin + alkaline + SGOT + prothrombin + age + sex + ascites +
                          spiders + edema + histologic, data = pbc2.id, x = TRUE, model = TRUE)
joint_fit_bilir2 <- jointModel(lme_fit_bilir, cox_fit_bilir2, timeVar = "year", method = "piecewise-PH-GH", iter.EM = 10)
summary(joint_fit_bilir2)

# drop spiders
cox_fit_bilir3 <- coxph(Surv(years, status2) ~ albumin + alkaline + prothrombin + age + sex + ascites +
                          spiders + edema + histologic, data = pbc2.id, x = TRUE, model = TRUE)
joint_fit_bilir3 <- jointModel(lme_fit_bilir, cox_fit_bilir3, timeVar = "year", method = "piecewise-PH-GH", iter.EM = 10)
summary(joint_fit_bilir3)

# drop ascites
cox_fit_bilir4 <- coxph(Surv(years, status2) ~ albumin + alkaline + prothrombin + age + sex + ascites +
                          edema + histologic, data = pbc2.id, x = TRUE, model = TRUE)
joint_fit_bilir4 <- jointModel(lme_fit_bilir, cox_fit_bilir4, timeVar = "year", method = "piecewise-PH-GH", iter.EM = 10)
summary(joint_fit_bilir4)

# drop histologic
cox_fit_bilir5 <- coxph(Surv(years, status2) ~ albumin + alkaline + prothrombin + age + sex +
                          edema + histologic, data = pbc2.id, x = TRUE, model = TRUE)
joint_fit_bilir5 <- jointModel(lme_fit_bilir, cox_fit_bilir5, timeVar = "year", method = "piecewise-PH-GH", iter.EM = 10)
summary(joint_fit_bilir5)

# drop sex
cox_fit_bilir6 <- coxph(Surv(years, status2) ~ albumin + alkaline + prothrombin + age + sex + edema, data = pbc2.id, x = TRUE, model = TRUE)
joint_fit_bilir6 <- jointModel(lme_fit_bilir, cox_fit_bilir6, timeVar = "year", method = "piecewise-PH-GH", iter.EM = 10)
summary(joint_fit_bilir6)

# drop prothrombin
cox_fit_bilir7 <- coxph(Surv(years, status2) ~ albumin + alkaline + prothrombin + age + edema, data = pbc2.id, x = TRUE, model = TRUE)
joint_fit_bilir7 <- jointModel(lme_fit_bilir, cox_fit_bilir7, timeVar = "year", method = "piecewise-PH-GH", iter.EM = 10)
summary(joint_fit_bilir7)

# drop alkaline
cox_fit_bilir8 <- coxph(Surv(years, status2) ~ albumin + alkaline + age + edema, data = pbc2.id, x = TRUE, model = TRUE)
joint_fit_bilir8 <- jointModel(lme_fit_bilir, cox_fit_bilir8, timeVar = "year", method = "piecewise-PH-GH", iter.EM = 10)
summary(joint_fit_bilir8)

# final joint model
cox_fit_bilir9 <- coxph(Surv(years, status2) ~ albumin + age + edema, data = pbc2.id, x = TRUE, model = TRUE)
joint_fit_bilir9 <- jointModel(lme_fit_bilir, cox_fit_bilir9, timeVar = "year", method = "piecewise-PH-GH", iter.EM = 100)
summary(joint_fit_bilir9)

# p-value for edema wit anova 1e-04
cox_fit_bilir10 <- coxph(Surv(years, status2) ~ albumin + age, data = pbc2.id, x = TRUE, model = TRUE)
joint_fit_bilir10 <- jointModel(lme_fit_bilir, cox_fit_bilir10, timeVar = "year", method = "piecewise-PH-GH", iter.EM = 100)
summary(joint_fit_bilir10)

exp(confint(joint_fit_bilir, parm = "Event"))

exp(1.65307287 / 3.346086)
exp(2.01092754 / 3.346086)
exp(2.36878221 / 3.346086)

### Diagnostic Tests for JM
# Wald test and Likelihood test
anova(joint_fit_bilir, process = "Longitudinal")
anova(joint_fit_bilir10, joint_fit_bilir)

# confidence intervals
confint(joint_fit, parm = "Longitudinal")
exp(confint(joint_fit, parm = "Event"))

# diagnostic for longitudinal
par(mfrow = c(2, 2))
plot(joint_fit)

# function to produce scatteplots with superimposed smooth line
plotResid <- function (x, y, col.loess = "black", ...) {
  plot(x, y, ...)
  lines(lowess(x, y), col = col.loess, lwd = 2)
  abline(h = 0, lty = 3, col = "grey", lwd = 2)
}

res_mar <- residuals(joint_fit, process = "Longitudinal", type = "Marginal")
fit_mar <- fitted(joint_fit, process = "Longitudinal", type = "Marginal")
plotResid(fit_mar, res_mar, xlab = "Fitted Values", ylab = "Marginal Residuals")

# residuals for survival part: martingale residuals 
res_mart <- residuals(joint_fit, process = "Event")
fit_mart <- fitted(joint_fit, process = "Longitudinal", type = "EventTime")

plotResid(fit_mart, res_mart, col.loess = "grey62", ylab = "Martingale Residuals",
          xlab = "Subject-Specific Fitted Values Longitudinal Outcome")

resCST <- residuals(joint_fit, process = "Event", type = "CoxSnell")
