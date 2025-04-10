library(ggplot2)
library(ggdag)
library(dplyr)

coord_dag <- list(
  x = c(PCB_PC1 = -60, PCB_PC2 = -65, PBDE_PC1 = -75, PBDE_PC2 = -70, PFDA_PC1 = -40, DDE = -50,  EtFOSSA = -60, 
        FOSA = -55, PFHpS= -25,  n_PFOS= -10, Sm_PFOS = -55, n_PFOA = -40,
        AGE = 40, BMI = 15, GLFM = -5, MALE_SEX = 60, SFM = -20, MEDS = 80,
        APN = 13, TRIG = 45, CHOL = 30, CRP= 60, GGT = 70, 
        HA1c = 100),
  

  y = c(PCB_PC1 = -65, PCB_PC2 = 40, PBDE_PC1= -38, PBDE_PC2 = 5, PFDA_PC1 = 10, DDE = 30,  EtFOSSA= -15, 
        FOSA = -35, PFHpS = -65,  n_PFOS = -65, Sm_PFOS=13, n_PFOA = -65,
        AGE = 80, BMI = 75, GLFM = 75, MALE_SEX = 75, SFM = 75, MEDS = 55,
        APN = -50,  TRIG = -65, CHOL = -70, CRP = -50, GGT = -65,
        HA1c = 10))

merge_table = cbind(c("PCB_PC1", "PCB_PC2", "PBDE_PC1", "PBDE_PC2", "PFDA_PC1", "DDE","EtFOSSA", 
        "FOSA", "PFHpS","n_PFOS", "Sm_PFOS", "n_PFOA", 
        "AGE", "BMI", "GLFM", "MALE_SEX", "SFM", "MEDS",
        "APN",  "TRIG", "CHOL", "CRP", "GGT","HA1c"), c(rep("Exposures", 12), rep("Confounders", 6), rep("Biomarkers", 5), "Outcome"))
        

my_dag <- dagify(PCB_PC1 ~ AGE + GLFM + MALE_SEX ,
                        PCB_PC2 ~ BMI, 
                        PBDE_PC1 ~ AGE ,
                        PFDA_PC1 ~ GLFM + MALE_SEX + SFM, 
                        PBDE_PC2 ~ PBDE_PC2, 
                        EtFOSSA ~ AGE, 
                        FOSA ~ AGE, 
                        PFHpS ~ AGE + GLFM + MALE_SEX,
                        n_PFOS ~ AGE + GLFM + MALE_SEX, 
                        Sm_PFOS ~ GLFM ,
                        n_PFOA ~ GLFM + MALE_SEX + SFM, 
                        APN ~ MALE_SEX, 
                        DDE ~ AGE + MALE_SEX + BMI + GLFM + MEDS,
                        #MEDS ~ DDE, 
                        TRIG ~ BMI,
                        CRP ~ AGE + BMI,
                        GGT ~ MALE_SEX, 
                        CHOL ~ MALE_SEX + SFM,
                        HA1c ~ APN + DDE + MEDS,  
                        coords = coord_dag, 
                        labels = c(PCB_PC1 = "PCB\nPC1", AGE = "AGE",
                                   GLFM = "GFM", MALE_SEX = "MALE_SEX", 
                                   PCB_PC2 = "PCB\nPC2"))%>%tidy_dagitty()
                        # mutate(type = c("biomarkers", rep("confounders", 12), rep("biomarkers",2), rep("exposures",3), 
                        # "biomarkers", 
                        # rep("confounders", 7), "outcome", rep("confounders", 2), "exposures", rep("confounders", 3),"exposures","biomarkers",rep("exposures",7), 
                        # rep("confounders",9))) 

#                        mutate(type = my_dag$data$type)

table1 = as.data.frame(my_dag$data$name)
table2 = as.data.frame(merge_table)
colnames(table1) = c("name")
colnames(table2) = c("name", "type")
label = right_join(table1, table2, by= "name")
var_type = label$type

my_dag <- my_dag %>%mutate(type = c(var_type))

#%>% mutate(colour = c(rep("red",4), rep("blue",4), rep("orange",5)))
#ggdag(my_dag, text = FALSE, use_labels = "labels") %>%
my_dag %>%  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
#  geom_dag_point(size = 15) +
  geom_dag_point(aes(colour = type), size = 25)+
#  geom_dag_point(color = "orange", size = 18)+
#  geom_dag_edges_link() +
  geom_dag_edges_arc(edge_color = "black", curvature = 0.11)+
# geom_dag_edges(aes(edge_colour = colour))+
  geom_dag_text(size = 3.1) +
#  geom_dag_label(aes(label = c("AGE", "BMI", "MALE_SEX")))+
  theme_dag()







library(grid)  # For using the arrow function

my_dag %>%
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_point(aes(colour = type), size = 30) +
  geom_dag_edges_arc(
    aes(edge_alpha = 0.5),  # Adjust edge transparency if needed
    curvature = 0.1,  # Adjust curvature as needed
    arrow = arrow(type = "closed", length = unit(2, "mm"))
  ) +
  geom_dag_text(size = 4.0) +
  theme_dag()