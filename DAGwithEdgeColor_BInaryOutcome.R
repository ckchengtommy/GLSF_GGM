library(ggplot2)
library(ggdag)
library(dplyr)

coord_dag <- list(
  x = c(pcb_pc1 = -4, pcb_pc2 = -6, pbde_pc1 = -8, pbde_pc2 = -8, pfda_pc1 = -4, DDE = -3,  EtFOSSA = -3, 
        FOSA = -1, PFHpS= -4,  n_PFOS= 0, Sm_PFOS = -6, n_PFOA = 1, 
        Age = 20, BMI = 6+3, GLFmeals = 3+3, sex = 10+3, SaltfishMeals = 18, 
        APN = 13, TRIG = 18, CHOL = 15, CRP=17,
        Diabetes = 30),
  
  
  y = c(pcb_pc1 = 0, pcb_pc2 = -7, pbde_pc1= -16, pbde_pc2 = -3, pfda_pc1 = 6, DDE = 12,  EtFOSSA= -8, 
        FOSA = -4, PFHpS = -16,  n_PFOS = -14, Sm_PFOS=13, n_PFOA = -18,
        Age = 20, BMI = 30, GLFmeals = 28, sex = 30, SaltfishMeals = 30,
        APN = -14,  TRIG = -12, CHOL = -15, CRP = -22,
        Diabetes = 0))

my_dag <- dagify(pcb_pc1 ~ Age + GLFmeals + sex ,
                 pcb_pc2 ~ BMI, 
                 pbde_pc1 ~ Age ,
                 pfda_pc1 ~ GLFmeals + sex + SaltfishMeals, 
                 pbde_pc2 ~ pbde_pc2, 
                 EtFOSSA ~ Age, 
                 FOSA ~ Age, 
                 PFHpS ~ Age + GLFmeals + sex,
                 n_PFOS ~ Age + GLFmeals + sex, 
                 Sm_PFOS ~ GLFmeals ,
                 n_PFOA ~ GLFmeals + sex + SaltfishMeals, 
                 APN ~ sex, 
                 DDE ~ Age + sex + BMI  + GLFmeals,
                 TRIG ~ BMI,
                 CRP ~ Age + BMI,
                 GGT ~ sex, 
                 CHOL ~ sex + SaltfishMeals,
                 Diabetes ~ DDE, 
                 coords = coord_dag, 
                 labels = c(pcb_pc1 = "PCB\nPC1", Age = "AGE",
                            GLFmeals = "GFM", sex = "MALE\nSEX", 
                            pcb_pc2 = "PCB\nPC2")) %>% tidy_dagitty() %>% mutate(type = c("biomarkers", rep("confounders", 12), rep("biomarkers",2), "exposures", "outcome", rep("exposures",2), "biomarkers", 
                                                                                          rep("confounders", 7), "exposures", rep("confounders",3), "exposures", "biomarkers", rep("exposures",7), rep("confounders",9))) 
#                                                   

#                                    mutate(type = c("bio", rep("cov, 6"), "bio", "exp", "bio", rep("exp",3), "out", rep("exp", 2), 
#                                                    ))

#%>% mutate(colour = c(rep("red",4), rep("blue",4), rep("orange",5)))
#ggdag(my_dag, text = FALSE, use_labels = "labels") %>%
my_dag %>%  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
  #  geom_dag_point(size = 15) + 
  geom_dag_point(aes(colour = type), size = 20)+
  #  geom_dag_point(color = "orange", size = 18)+
  #  geom_dag_edges_link() +
  geom_dag_edges_arc(edge_color = "black", curvature = 0.1)+
  # geom_dag_edges(aes(edge_colour = colour))+
  geom_dag_text(size = 3) +
  #  geom_dag_label(aes(label = c("Age", "BMI", "SEX")))+
  theme_dag()