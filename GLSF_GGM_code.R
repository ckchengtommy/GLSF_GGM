library(glmnet); library(purrr)
library(mice)
library(sas7bdat)
library(huge); library(GGally)
library(mgm); library(dplyr)
library(network); library(corrplot)
library(sna)
library(ggplot2); library(naniar)
library(SILGGM)
library(dplyr)
library(rms)
library(igraph)
library(haven)


#set.seed(538)
file_path <- file.choose() 
GLdata = read_sas(file_path)
GLdata = as_tibble(GLdata)
GLdata <- GLdata %>% rename(YearsGLF= yearglf3, HA1c = hemoglob, Age = agef4, 
                            BMI = bmif4, GLFmeals = glmealsf4, SaltfishMeals = saltfishf4, sex = gender, 
                            DiabetesMed = diabmedf4, TRIG = triglycerf4, CHOL = cholestf4,
                            DDE = ddeimpla, pcb180 = pcb180impla, pcb163_138 = pcb163138impla, 
                            pcb132_153_105 = pcb132153105impla, pcb194 = pcb194impla, 
                            pcb201 = pcb201impla, pcb170_190 = pcb170190impla, 
                            pcb187_182 = pcb187182impla, pbde47 = pbde47impla, 
                            pbde99 = pbde99impla,
                            APN = apn, CRP = crp, Diabetes = diabetes2level)

GLdata = as.data.frame(GLdata)

#check missing data
prop_miss(GLdata)
GLdata %>% is.na() %>% colSums()
miss_var_summary(GLdata)
miss_var_table(GLdata)
gg_miss_var(GLdata)
head(GLdata)

GL = subset(GLdata, select= -c(caseid))
GL[GL== "NaN"] = NA


GL_CC = GL[complete.cases(GL), ]

GL_pcb = GL_CC%>% select(pcb180, pcb163_138, pcb132_153_105, pcb194, pcb201, pcb170_190, pcb187_182)
corrplot(cor(GL_pcb),  cl.cex = 0.7, tl.cex = 0.7, method = "number")

var_names = colnames(GL)


#Rubin's Rule function 
combine_p_val <- function(p_vec){
  z = qnorm(p_vec/2)
  v = 1+var(z)
  z0 = mean(z)/sqrt(v)
  pvalue = pnorm(z0)
  pvalue
}


cont_var = var_names[-which(names(GL) %in% c("education", "sex", "Diabetes_stage", 
                                             "DiabetesMed","saltfishqt", "glmealsqt", 
                                             "Diabetes"))]

#Perform log transformation on continuous variables 
GL[, cont_var] = log(GL[, cont_var])

miss_var_names = var_names[colSums(is.na(GL))>0]
miss_cont_var = intersect(miss_var_names, cont_var)
miss_cont_var1 = c("DDE", "pcb180", "pcb163_138", "pcb132_153_105",
                   "pcb194", "pcb201", "pcb170_190", "pcb187_182")
miss_cont_var2 = miss_cont_var[-which(miss_cont_var %in% miss_cont_var1)]

GL <- GL%>% mutate(
  #  education = factor(education, order = TRUE), 
  glmealsqt = factor(glmealsqt, order = TRUE),
  saltfishqt = factor(saltfishqt, order = TRUE), 
  #Diabetes_stage = factor(Diabetes_stage, order = TRUE),
  DiabetesMed = as.factor(DiabetesMed),
  sex = as.factor(sex),
  Diabetes = as.factor(Diabetes)
)

init = mice(GL, maxit = 0)
meth = init$method
#meth["education"] = "polr"
#meth[miss_cont_var1] = "norm.boot"
#meth[miss_cont_var2] = "norm"
meth[miss_cont_var] = "norm"
#meth[miss_cont_var] = "norm.boot"
GL_imp = mice(GL, m=5, method = meth, maxit = 100, seed = 1555)
imputelist = list()
for(i in 1:5){
  imputelist[[i]] = as.data.frame(complete(GL_imp, action = i)) %>% 
    select(-c(corticosteroidf4, diagdiabf4, saltfishqt,glmealsqt )) %>% 
    mutate_if(is.factor, as.numeric) %>% mutate(Diabetes = ifelse(Diabetes == 2, 1, 0)) %>%
    mutate(DiabetesMed = ifelse(DiabetesMed == 2, 1, 0))
}


#Z-scale transformation
for(i in 1:5){
  imputelist[[i]][,1:34] = scale(imputelist[[i]][,1:34])
}


GL_list = list()
for(i in 1:5){
  GL_list[[i]] = as.data.frame(imputelist[[i]])
}


#PCA

PCB_group = list()
pbde_group = list()
pfda_group = list()

for(i in 1:5){
  PCB_group[[i]] = GL_list[[i]]%>% select(pcb180, pcb163_138, pcb132_153_105,
                                          pcb194, pcb201,
                                          pcb170_190, pcb187_182)
  pbde_group[[i]] = GL_list[[i]]%>% select(pbde47, pbde99)
  pfda_group[[i]] = GL_list[[i]]%>%select(PFDA, PFUnDA)
}

PCB_g2 = do.call("rbind", PCB_group)
pbde_g2 = do.call("rbind", pbde_group)
pfda_g2 = do.call("rbind", pfda_group)


pcb_pca <- prcomp(PCB_g2, scale = TRUE,
                  center = TRUE, retx = T)
#summary(pcb_pca)
pbde_pca <- prcomp(pbde_g2, scale = TRUE,
                   center = TRUE, retx = T)
#summary(pbde_pca)
pfda_pca <- prcomp(pfda_g2, scale = TRUE,
                   center = TRUE, retx = T)
#summary(pfda_pca)


pcb_pc1 = pcb_pca$x[, 1]
pcb_pc2 = pcb_pca$x[, 2]
pbde_pc1 = pbde_pca$x[, 1]
pbde_pc2 = pbde_pca$x[, 2]
pfda_pc1 = pfda_pca$x[, 1]
pfda_pc2 = pfda_pca$x[, 2 ]

GL_list[[1]] = as.data.frame(cbind(GL_list[[1]], pcb_pc1 = pcb_pc1[1:513], pbde_pc1 = pbde_pc1[1:513], pfda_pc1 =pfda_pc1[1:513], 
                                   pcb_pc2 = pcb_pc2[1:513], pbde_pc2 = pbde_pc2[1:513], pfda_pc2= pfda_pc2[1:513]))
GL_list[[2]] = as.data.frame(cbind(GL_list[[2]], pcb_pc1 =pcb_pc1[514:1026], pbde_pc1=pbde_pc1[514:1026], pfda_pc1= pfda_pc1[514:1026], 
                                   pcb_pc2=pcb_pc2[514:1026], pbde_pc2=pbde_pc2[514:1026], pfda_pc2=pfda_pc2[514:1026]))
GL_list[[3]] = as.data.frame(cbind(GL_list[[3]], pcb_pc1= pcb_pc1[1027:1539], pbde_pc1=pbde_pc1[1027:1539], pfda_pc1= pfda_pc1[1027:1539], 
                                   pcb_pc2=pcb_pc2[1027:1539], pbde_pc2= pbde_pc2[1027:1539], pfda_pc2= pfda_pc2[1027:1539]))
GL_list[[4]] = as.data.frame(cbind(GL_list[[4]], pcb_pc1=pcb_pc1[1540:2052], pbde_pc1=pbde_pc1[1540:2052], pfda_pc1=pfda_pc1[1540:2052], 
                                   pcb_pc2=pcb_pc2[1540:2052], pbde_pc2=pbde_pc2[1540:2052], pfda_pc2=pfda_pc2[1540:2052]))
GL_list[[5]] = as.data.frame(cbind(GL_list[[5]], pcb_pc1=pcb_pc1[2053:2565], pbde_pc1=pbde_pc1[2053:2565], pfda_pc1=pfda_pc1[2053:2565], 
                                   pcb_pc2=pcb_pc2[2053:2565], pbde_pc2=pbde_pc2[2053:2565], pfda_pc2=pfda_pc2[2053:2565]))



confounders_tag = c("Age", "sex", "BMI", "GLFmeals", "SaltfishMeals", "DiabetesMed")
confounders_tag2 = c("Age", "sex", "BMI", "GLFmeals", "SaltfishMeals")
biomarkers_tag = c("APN", "GGT", "CRP", "TRIG", "CHOL")
exposures_tag = c("pcb_pc1", "pcb_pc2", "pbde_pc1", "pbde_pc2", "pfda_pc1", 
                  "pfda_pc2", "EtFOSAA", "FOSA", "DDE","PFHpS", "PFNA", "PFHxS", "MeFOSAA", 
                  "Sm_PFOS", "n_PFOS", "n_PFOA")
outcome1_tag = "HA1c"
outcome2_tag = "Diabetes"


Con_df = list()
Con_Expo_df = list()
Con_Expo_Med_df = list()
Con_Expo_Med_out1_df = list()
Con_Expo_Med_out2_df = list()


#GL_list[[1]] %>% select(all_of(confounders_tag))

for(i in 1:5){
  Con_df[[i]]  = data.matrix(GL_list[[i]] %>% select(all_of(confounders_tag)))
  Con_Expo_df[[i]] = data.matrix(GL_list[[i]] %>% select(all_of(c(confounders_tag, exposures_tag))))
  Con_Expo_Med_df[[i]] = data.matrix(GL_list[[i]] %>% select(all_of(c(confounders_tag, exposures_tag, biomarkers_tag))))
  Con_Expo_Med_out1_df[[i]] = data.matrix(GL_list[[i]] %>% select(all_of(c(confounders_tag, exposures_tag, biomarkers_tag, outcome1_tag))))
  Con_Expo_Med_out2_df[[i]] = data.matrix(GL_list[[i]] %>% select(all_of(c(confounders_tag2, exposures_tag, biomarkers_tag, outcome2_tag))))
}

GGM_Con_label = colnames(Con_df[[1]])
GGM_Con_Expo_label = colnames(Con_Expo_df[[1]])
GGM_Con_Expo_Med_label = colnames(Con_Expo_Med_df[[1]])
GGM_Con_Expo_Med_out1_label = colnames(Con_Expo_Med_out1_df[[1]])
GGM_Con_Expo_Med_out2_label = colnames(Con_Expo_Med_out2_df[[1]])

GGM_Con_Expo_Med_out2 = GGM_Con_Expo_Med_out1 = GGM_Con_Expo_Med = GGM_Con_Expo = GGM_Con = list()
#GGM_Con_Expo = list()

for(i in 1:5){
  GGM_Con[[i]] = SILGGM(Con_df[[i]], method= "B_NW_SL", global = FALSE)
  GGM_Con_Expo[[i]] = SILGGM(Con_Expo_df[[i]], method= "B_NW_SL", global = FALSE)
  GGM_Con_Expo_Med[[i]] = SILGGM(Con_Expo_Med_df[[i]], method= "B_NW_SL", global = FALSE)
  GGM_Con_Expo_Med_out1[[i]] = SILGGM(Con_Expo_Med_out1_df[[i]], method= "B_NW_SL", global = FALSE)
  GGM_Con_Expo_Med_out2[[i]] = SILGGM(Con_Expo_Med_out2_df[[i]], method= "B_NW_SL", global = FALSE)
}

k = GGM_Con_Expo_Med_out2[[2]]$p_precision
rownames(k) = colnames(k) = GGM_Con_Expo_Med_out2_label

sort(k[,27][-27])[3] < 0.05 * (4 / (26*2))



Con_pval = Con_Expo_pval = Con_Expo_Med_pval = Con_Expo_Med_out1_pval = list()
for(i in 1:5){
  test = GGM_Con_Expo[[i]]$p_partialCor
}



#Extract p-value of the edge
Con_Expo_pval = GGM_Con_Expo %>% lapply(. %>% pluck("p_partialCor")) %>% lapply(. %>% {colnames(.) = GGM_Con_Expo_label;.}) %>% 
  lapply(. %>% {rownames(.) = GGM_Con_Expo_label;.})       

Con_Expo_Med_pval = GGM_Con_Expo_Med %>% lapply(. %>% pluck("p_partialCor")) %>% lapply(. %>% {colnames(.) = GGM_Con_Expo_Med_label;.}) %>% 
  lapply(. %>% {rownames(.) = GGM_Con_Expo_Med_label;.}) 

Con_Expo_Med_out1_pval = GGM_Con_Expo_Med_out1 %>% lapply(. %>% pluck("p_partialCor")) %>% lapply(. %>% {colnames(.) = GGM_Con_Expo_Med_out1_label;.}) %>% 
  lapply(. %>% {rownames(.) = GGM_Con_Expo_Med_out1_label;.})

Con_Expo_Med_out2_pval = GGM_Con_Expo_Med_out2 %>% lapply(. %>% pluck("p_partialCor")) %>% lapply(. %>% {colnames(.) = GGM_Con_Expo_Med_out2_label;.}) %>% 
  lapply(. %>% {rownames(.) = GGM_Con_Expo_Med_out2_label;.})


#Extract effect size of the edge
Con_Expo_partialCor = GGM_Con_Expo %>% lapply(. %>% pluck("partialCor")) %>% lapply(. %>% {colnames(.) = GGM_Con_Expo_label;.}) %>% 
  lapply(. %>% {rownames(.) = GGM_Con_Expo_label;.})       

Con_Expo_Med_partialCor = GGM_Con_Expo_Med %>% lapply(. %>% pluck("partialCor")) %>% lapply(. %>% {colnames(.) = GGM_Con_Expo_Med_label;.}) %>% 
  lapply(. %>% {rownames(.) = GGM_Con_Expo_Med_label;.}) 

Con_Expo_Med_out1_partialCor = GGM_Con_Expo_Med_out1 %>% lapply(. %>% pluck("partialCor")) %>% lapply(. %>% {colnames(.) = GGM_Con_Expo_Med_out1_label;.}) %>% 
  lapply(. %>% {rownames(.) = GGM_Con_Expo_Med_out1_label;.})

Con_Expo_Med_out2_partialCor = GGM_Con_Expo_Med_out2 %>% lapply(. %>% pluck("partialCor")) %>% lapply(. %>% {colnames(.) = GGM_Con_Expo_Med_out2_label;.}) %>% 
  lapply(. %>% {rownames(.) = GGM_Con_Expo_Med_out2_label;.})


Con_Expo_ES = Reduce("+", Con_Expo_partialCor) / length(Con_Expo_partialCor)
Con_Expo_Med

pmat1.1 = Con_Expo_pval[[1]][1:6, 7:22] 
pmat2.1 = Con_Expo_pval[[2]][1:6, 7:22]
pmat3.1 = Con_Expo_pval[[3]][1:6, 7:22] 
pmat4.1 = Con_Expo_pval[[4]][1:6, 7:22] 
pmat5.1 = Con_Expo_pval[[5]][1:6, 7:22] 


pmat1.1 = Con_Expo_pval[[1]][confounders_tag, exposures_tag] 
pmat2.1 = Con_Expo_pval[[2]][confounders_tag, exposures_tag]
pmat3.1 = Con_Expo_pval[[3]][confounders_tag, exposures_tag] 
pmat4.1 = Con_Expo_pval[[4]][confounders_tag, exposures_tag] 
pmat5.1 = Con_Expo_pval[[5]][confounders_tag, exposures_tag] 


# pcb_pc1 = rbind(pmat1.1[,1], pmat2.1[,1], pmat3.1[,1], pmat4.1[,1], pmat5.1[,1])
# pcb_pc2 = rbind(pmat1.1[,2], pmat2.1[,2], pmat3.1[,2], pmat4.1[,2], pmat5.1[,2])
# pbde_pc1= rbind(pmat1.1[,3], pmat2.1[,3], pmat3.1[,3], pmat4.1[,3], pmat5.1[,3])
# pbde_pc2= rbind(pmat1.1[,4], pmat2.1[,4], pmat3.1[,4], pmat4.1[,4], pmat5.1[,4])
# pfda_pc1= rbind(pmat1.1[,5], pmat2.1[,5], pmat3.1[,5], pmat4.1[,5], pmat5.1[,5])
# pfda_pc2= rbind(pmat1.1[,6], pmat2.1[,6], pmat3.1[,6], pmat4.1[,6], pmat5.1[,6])
# EtFOSAA = rbind(pmat1.1[,7], pmat2.1[,7], pmat3.1[,7], pmat4.1[,7], pmat5.1[,7])
# FOSA    = rbind(pmat1.1[,8], pmat2.1[,8], pmat3.1[,8], pmat4.1[,8], pmat5.1[,8])
# DDE     = rbind(pmat1.1[,9], pmat2.1[,9], pmat3.1[,9], pmat4.1[,9], pmat5.1[,9])
# PFHpS   = rbind(pmat1.1[,10], pmat2.1[,10], pmat3.1[,10], pmat4.1[,10], pmat5.1[,10])
# PFNA    = rbind(pmat1.1[,11], pmat2.1[,11], pmat3.1[,11], pmat4.1[,11], pmat5.1[,11])
# PFHxS   = rbind(pmat1.1[,12], pmat2.1[,12], pmat3.1[,12], pmat4.1[,12], pmat5.1[,12])    
# MeFOSAA = rbind(pmat1.1[,13], pmat2.1[,13], pmat3.1[,13], pmat4.1[,13], pmat5.1[,13])
# Sm_PFOS = rbind(pmat1.1[,14], pmat2.1[,14], pmat3.1[,14], pmat4.1[,14], pmat5.1[,14])      
# n_PFOS  = rbind(pmat1.1[,15], pmat2.1[,15], pmat3.1[,15], pmat4.1[,15], pmat5.1[,15])      
# n_PFOA  = rbind(pmat1.1[,16], pmat2.1[,16], pmat3.1[,16], pmat4.1[,16], pmat5.1[,16])

Graph1_pcb_pc1_pval = rbind(pmat1.1[,1], pmat2.1[,1], pmat3.1[,1], pmat4.1[,1], pmat5.1[,1])
Graph1_pcb_pc2_pval = rbind(pmat1.1[,2], pmat2.1[,2], pmat3.1[,2], pmat4.1[,2], pmat5.1[,2])
Graph1_pbde_pc1_pval = rbind(pmat1.1[,3], pmat2.1[,3], pmat3.1[,3], pmat4.1[,3], pmat5.1[,3])
Graph1_pbde_pc2_pval = rbind(pmat1.1[,4], pmat2.1[,4], pmat3.1[,4], pmat4.1[,4], pmat5.1[,4])
Graph1_pfda_pc1_pval = rbind(pmat1.1[,5], pmat2.1[,5], pmat3.1[,5], pmat4.1[,5], pmat5.1[,5])
Graph1_pfda_pc2_pval = rbind(pmat1.1[,6], pmat2.1[,6], pmat3.1[,6], pmat4.1[,6], pmat5.1[,6])
Graph1_EtFOSAA_pval = rbind(pmat1.1[,7], pmat2.1[,7], pmat3.1[,7], pmat4.1[,7], pmat5.1[,7])
Graph1_FOSA_pval = rbind(pmat1.1[,8], pmat2.1[,8], pmat3.1[,8], pmat4.1[,8], pmat5.1[,8])
Graph1_DDE_pval = rbind(pmat1.1[,9], pmat2.1[,9], pmat3.1[,9], pmat4.1[,9], pmat5.1[,9])
Graph1_PFHpS_pval = rbind(pmat1.1[,10], pmat2.1[,10], pmat3.1[,10], pmat4.1[,10], pmat5.1[,10])
Graph1_PFNA_pval = rbind(pmat1.1[,11], pmat2.1[,11], pmat3.1[,11], pmat4.1[,11], pmat5.1[,11])
Graph1_PFHxS_pval = rbind(pmat1.1[,12], pmat2.1[,12], pmat3.1[,12], pmat4.1[,12], pmat5.1[,12])    
Graph1_MeFOSAA_pval = rbind(pmat1.1[,13], pmat2.1[,13], pmat3.1[,13], pmat4.1[,13], pmat5.1[,13])
Graph1_Sm_PFOS_pval = rbind(pmat1.1[,14], pmat2.1[,14], pmat3.1[,14], pmat4.1[,14], pmat5.1[,14])      
Graph1_n_PFOS_pval = rbind(pmat1.1[,15], pmat2.1[,15], pmat3.1[,15], pmat4.1[,15], pmat5.1[,15])      
Graph1_n_PFOA_pval = rbind(pmat1.1[,16], pmat2.1[,16], pmat3.1[,16], pmat4.1[,16], pmat5.1[,16])



sort(pmat1.1[,10])
sort(pmat2.1[,10])
sort(pmat3.1[,10])
sort(pmat4.1[,10])
sort(pmat5.1[,10])

pmat1.2 = Con_Expo_Med_pval[[1]][1:22 ,c("APN", "GGT", "CRP", "TRIG", "CHOL")]
pmat2.2 = Con_Expo_Med_pval[[2]][1:22, c("APN", "GGT", "CRP", "TRIG", "CHOL")]
pmat3.2 = Con_Expo_Med_pval[[3]][1:22, c("APN", "GGT", "CRP", "TRIG", "CHOL")]
pmat4.2 = Con_Expo_Med_pval[[4]][1:22, c("APN", "GGT", "CRP", "TRIG", "CHOL")]
pmat5.2 = Con_Expo_Med_pval[[5]][1:22, c("APN", "GGT", "CRP", "TRIG", "CHOL")] 


pcb_pc1_link = pcb_pc2_link = pbde_pc1_link =  pbde_pc1_link = pbde_pc2_link = pfda_pc1_link = pfda_pc2_link = c()
EtFOSAA_link = FOSA_link = DDE_link = PFHpS_link = PFNA_link = PFHxS_link = MeFOSAA_link = Sm_PFOS_link = n_PFOS_link = n_PFOA_link = c()
Age_link = sex_link = BMI_link = GLFmeals_link = SaltfishMeals_link = DiabetesMed_link = c()

pcb_pc1_comb = pcb_pc2_comb = pbde_pc1_comb =  pbde_pc1_comb = pbde_pc2_comb = pfda_pc1_comb = pfda_pc2_comb = c()
EtFOSAA_comb = FOSA_comb = DDE_comb = PFHpS_comb = PFNA_comb = PFHxS_comb = MeFOSAA_comb = Sm_PFOS_comb = n_PFOS_comb = n_PFOA_comb = c()
Age_comb = sex_comb = BMI_comb = GLFmeals_comb = SaltfishMeals_comb = DiabetesMed_comb = c()


#Graph: Confounders + Exposures on Mediators
for(i in 1:5){
  pcb_pc1_comb[i] = combine_p_val(pcb_pc1[,i])
  pcb_pc2_comb[i] = combine_p_val(pcb_pc2[,i])
  pbde_pc1_comb[i] = combine_p_val(pbde_pc1[,i])
  pbde_pc2_comb[i] = combine_p_val(pbde_pc2[,i])
  pfda_pc1_comb[i] = combine_p_val(pfda_pc1[,i])
  pfda_pc2_comb[i] = combine_p_val(pfda_pc2[,i])
  EtFOSAA_comb[i] = combine_p_val(EtFOSAA[,i])
  FOSA_comb[i] = combine_p_val(FOSA[,i])
  DDE_comb[i] = combine_p_val(DDE[,i])
  PFHpS_comb[i] = combine_p_val(PFHpS[,i])
  PFNA_comb[i] = combine_p_val(PFNA[,i])
  PFHxS_comb[i] = combine_p_val(PFHxS[,i])
  MeFOSAA_comb[i] = combine_p_val(MeFOSAA[,i])
  Sm_PFOS_comb[i] = combine_p_val(Sm_PFOS[,i])
  n_PFOS_comb[i] = combine_p_val(n_PFOS[,i])
  n_PFOA_comb[i] = combine_p_val(n_PFOA[,i])
}


#Graph : Confounders on Exposures
for(i in 1:6){
  pcb_pc1_comb[i] = combine_p_val(pcb_pc1[,i])
  pcb_pc2_comb[i] = combine_p_val(pcb_pc2[,i])
  pbde_pc1_comb[i] = combine_p_val(pbde_pc1[,i])
  pbde_pc2_comb[i] = combine_p_val(pbde_pc2[,i])
  pfda_pc1_comb[i] = combine_p_val(pfda_pc1[,i])
  pfda_pc2_comb[i] = combine_p_val(pfda_pc2[,i])
  EtFOSAA_comb[i] = combine_p_val(EtFOSAA[,i])
  FOSA_comb[i] = combine_p_val(FOSA[,i])
  DDE_comb[i] = combine_p_val(DDE[,i])
  PFHpS_comb[i] = combine_p_val(PFHpS[,i])
  PFNA_comb[i] = combine_p_val(PFNA[,i])
  PFHxS_comb[i] = combine_p_val(PFHxS[,i])
  MeFOSAA_comb[i] = combine_p_val(MeFOSAA[,i])
  Sm_PFOS_comb[i] = combine_p_val(Sm_PFOS[,i])
  n_PFOS_comb[i] = combine_p_val(n_PFOS[,i])
  n_PFOA_comb[i] = combine_p_val(n_PFOA[,i])
}


for(i in 1:6){
  pcb_pc1_link[i]= ifelse(pcb_pc1_comb[i] < 0.05 * (i/(6*2)), 1, 0)
  pcb_pc2_link[i] = ifelse(pcb_pc2_comb[i] < 0.05 * (i/(6*2)), 1, 0)
  pbde_pc1_link[i] = ifelse(pbde_pc1_comb[i] < 0.05 * (i/(6*2)), 1, 0)
  pbde_pc2_link[i] = ifelse(pbde_pc2_comb[i] < 0.05 * (i/(6*2)), 1, 0)
  pfda_pc1_link[i] = ifelse(pfda_pc1_comb[i] < 0.05 * (i/(6*2)), 1, 0)
  pfda_pc2_link[i] = ifelse(pfda_pc2_comb[i] < 0.05 * (i/(6*2)), 1, 0)
  EtFOSAA_link[i] = ifelse(EtFOSAA_comb[i] < 0.05 * (i/(6*2)), 1, 0)
  FOSA_link[i] = ifelse(FOSA_comb[i] < 0.05 * (i/(6*2)), 1, 0)
  DDE_link[i] = ifelse(DDE_comb[i] < 0.05 * (i/(6*2)), 1, 0)
  PFHpS_link[i] = ifelse(PFHpS_comb[i] < 0.05 * (i/(6*2)), 1, 0)
  PFNA_link[i] = ifelse(PFNA_comb[i] < 0.05 * (i/(6*2)), 1, 0)
  PFHxS_link[i] = ifelse(PFHxS_comb[i] < 0.05 * (i/(6*2)), 1, 0)
  MeFOSAA_link[i] = ifelse(MeFOSAA_comb[i] < 0.05 * (i/(6*2)), 1, 0)
  Sm_PFOS_link[i] = ifelse(Sm_PFOS_comb[i] < 0.05 * (i/(6*2)), 1, 0)
  n_PFOS_link[i] = ifelse(n_PFOS_comb[i] < 0.05 * (i/(6*2)), 1, 0)
  n_PFOA_link[i] = ifelse(n_PFOA_comb[i] < 0.05 * (i/(6*2)), 1, 0)
}


names(pcb_pc1_link) = colnames(pcb_pc1)
names(pcb_pc2_link) = colnames(pcb_pc2)
names(pbde_pc1_link) =  colnames(pbde_pc1)
names(pbde_pc2_link) = colnames(pbde_pc2)
names(pfda_pc1_link) = colnames(pfda_pc1)
names(pfda_pc2_link) = colnames(pfda_pc2)
names(EtFOSAA_link) = colnames(EtFOSAA)
names(FOSA_link) = colnames(FOSA)
names(DDE_link) = colnames(DDE)
names(PFHpS_link) = colnames(PFHpS)
names(PFNA_link ) = colnames(PFNA)
names(PFHxS_link) = colnames(PFHxS)
names(MeFOSAA_link )= colnames(MeFOSAA)
names(Sm_PFOS_link) = colnames(Sm_PFOS)
names(n_PFOS_link) = colnames(n_PFOS)
names(n_PFOA_link) = colnames(n_PFOA)



pcb_pc1_link 
pcb_pc2_link
pbde_pc1_link 
pbde_pc2_link 
pfda_pc1_link 
pfda_pc2_link 
EtFOSAA_link 
FOSA_link 
DDE_link 
PFHpS_link
PFNA_link 
PFHxS_link
MeFOSAA_link
Sm_PFOS_link 
n_PFOS_link
n_PFOA_link 


rbind(select1.1, select2.1, select3.1, select4.1, select5.1)

pmat1 = sort(Con_Expo_Med_out2_pval[[1]][,27]) 
pmat2 = sort(Con_Expo_Med_out2_pval[[2]][,27])
pmat3 = sort(Con_Expo_Med_out2_pval[[3]][,27])
pmat4 = sort(Con_Expo_Med_out2_pval[[4]][,27])
pmat5 = sort(Con_Expo_Med_out2_pval[[5]][,27])

#FDR
select1= select2 = select3 = select4 = select5 = c()

for(i in 1:26){
  select1[i]= ifelse(Con_Expo_Med_out2_pval[[1]][,27][-1][i] < 0.05 * (i/(26*2)), 1, 0)
  select2[i]= ifelse(Con_Expo_Med_out2_pval[[2]][,27][-1][i] < 0.05 * (i/(26*2)), 1, 0)
  select3[i]= ifelse(Con_Expo_Med_out2_pval[[3]][,27][-1][i] < 0.05 * (i/(26*2)), 1, 0)
  select4[i]= ifelse(Con_Expo_Med_out2_pval[[4]][,27][-1][i] < 0.05 * (i/(26*2)), 1, 0)
  select5[i]= ifelse(Con_Expo_Med_out2_pval[[5]][,27][-1][i] < 0.05 * (i/(26*2)), 1, 0)
}

names(select2) = names(sort(Con_Expo_Med_out2_pval[[2]][,27])[-1])
names(select3) = names(sort(Con_Expo_Med_out2_pval[[3]][,27])[-1])

pvec = c(sort(Con_Expo_Med_out2_pval[[1]][,27])["DDE"], sort(Con_Expo_Med_out2_pval[[2]][,27]["DDE"]), sort(Con_Expo_Med_out2_pval[[3]][,27])["DDE"], 
         sort(Con_Expo_Med_out2_pval[[4]][,27])["DDE"], sort(Con_Expo_Med_out2_pval[[5]][,27]["DDE"]))





