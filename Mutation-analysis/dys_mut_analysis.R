opar <- par
par(opar)
mut.df <- read.csv(file.choose(), header = TRUE)
library("ggplot2")
library("dplyr")

mut.df$Reported[is.na(mut.df$Reported)] <- 0

# Defining the protein Regions as given by uniprot
mut.df$Region <- c(rep("ABD1", 249), rep("Undefined", 1177), rep("SYNM binding", 499),
          rep("Undefined", 1145), rep("SYNM binding", 249),
          rep("Zinc finger", 53), rep("SYNM binding", 55),
          rep("Undefined", 57), rep("SNTB1 binding", 53),
          rep("Undefined", 167))


# Defining the protein domains as given by uniprot. linker regions were only kept if longer than 8 amino acids 
# linkers less than 8 were split between the two sandwiching domains
mut.df$Domain <- c(rep("Nterm", 14), rep("CH1", 109), rep("CH1/CH2 linker", 14),
                 rep("CH2", 112), rep("CH2/S1 linker", 98), rep("S1", 111),
                 rep("S2", 112), rep("S3", 110), rep("S3/S4 linker", 51), rep("S4", 111), rep("S5", 105),
                 rep("S5/S6 linker", 8), rep("S6", 104), rep("S7", 110), rep("S8", 109), rep("S9", 103),
                 rep("S10", 98), rep("S11", 104), rep("S12", 108), rep("S13", 101), rep("S14", 97),
                 rep("S15", 104), rep("S15/S16 linker", 12), rep("S16", 111), rep("S17", 107), rep("S18", 109),
                 rep("S19", 107), rep("S19/S20 linker", 51), rep("S20", 104), rep("S21", 109), rep("S22", 118),
                 rep("S23", 130), rep("S24", 108), rep("S24/WW linker", 14), rep("WW", 34), rep("Cterm", 597))
mut.df$polar <- ifelse(  (mut.df$WT.AA == "S") | (mut.df$WT.AA == "T") |
                           (mut.df$WT.AA == "C") | (mut.df$WT.AA == "Y") |
                           (mut.df$WT.AA == "N") | (mut.df$WT.AA == "K") |
                           (mut.df$WT.AA == "Q") | (mut.df$WT.AA == "R") |
                           (mut.df$WT.AA == "D") | (mut.df$WT.AA == "H") |
                           (mut.df$WT.AA == "E"), 1, 0)

mut.df$charged <- ifelse(  (mut.df$WT.AA == "D") |
                             (mut.df$WT.AA == "E"), 1,
                           ifelse((mut.df$WT.AA == "K")| (mut.df$WT.AA == "R") | (mut.df$WT.AA == "H"),
                                  -1, 0))
mut.df$MW <- ifelse(mut.df$WT.AA == "G", 75.07,
ifelse(mut.df$WT.AA == "A", 89.09,
ifelse(mut.df$WT.AA == "S", 105.09,
ifelse(mut.df$WT.AA == "T", 119.1,
ifelse(mut.df$WT.AA == "C", 121.2,
ifelse(mut.df$WT.AA == "V", 117.1,
ifelse(mut.df$WT.AA == "L", 131.2,
ifelse(mut.df$WT.AA == "I", 131.2,
ifelse(mut.df$WT.AA == "M", 149.2,
ifelse(mut.df$WT.AA == "P", 115.1,
ifelse(mut.df$WT.AA == "F", 165.2,
ifelse(mut.df$WT.AA == "Y", 181.2,
ifelse(mut.df$WT.AA == "W", 204.2,
ifelse(mut.df$WT.AA == "D", 133.1,
ifelse(mut.df$WT.AA == "E", 147.1,
ifelse(mut.df$WT.AA == "N", 132.1,
ifelse(mut.df$WT.AA == "Q", 146.1,
ifelse(mut.df$WT.AA == "H", 155.2,
ifelse(mut.df$WT.AA == "K", 146.2, 174.2)))))))))))))))))))

mut.df$cp <- ifelse(mut.df$WT.AA == "V", 0.03115,
ifelse(mut.df$WT.AA == "A", 0.0196875,
ifelse(mut.df$WT.AA == "R", 0.0219875,
ifelse(mut.df$WT.AA == "N", 0.011375,
ifelse(mut.df$WT.AA == "D", 0.0126,
ifelse(mut.df$WT.AA == "C", 0.0118125,
ifelse(mut.df$WT.AA == "Q", 0.0134375,
ifelse(mut.df$WT.AA == "E", 0.015025,
ifelse(mut.df$WT.AA == "G", 0.0113125,
ifelse(mut.df$WT.AA == "H", 0.0257875,
ifelse(mut.df$WT.AA == "I", 0.037775,
ifelse(mut.df$WT.AA == "L", 0.0382875,
ifelse(mut.df$WT.AA == "K", 0.0288375,
ifelse(mut.df$WT.AA == "M", 0.030425,
ifelse(mut.df$WT.AA == "F", 0.0427125,
ifelse(mut.df$WT.AA == "P", 0.0310625,
ifelse(mut.df$WT.AA == "S", 0.014375,
ifelse(mut.df$WT.AA == "T", 0.0219875,
ifelse(mut.df$WT.AA == "W", 0.0453625,
ifelse(mut.df$WT.AA == "Y", 0.0352625,
print("There is at least one entry that isn't an amino acid")))))))))))))))))))))
mut.df$cp <- as.numeric(mut.df$cp)

mut.df$dipole <- ifelse(mut.df$WT.AA == "V", 2.692,
ifelse(mut.df$WT.AA == "A", 5.937,
ifelse(mut.df$WT.AA == "R", 37.5,
ifelse(mut.df$WT.AA == "N", 18.89,
ifelse(mut.df$WT.AA == "D", 29.49,
ifelse(mut.df$WT.AA == "C", 10.74,
ifelse(mut.df$WT.AA == "Q", 39.89,
ifelse(mut.df$WT.AA == "E", 42.52,
ifelse(mut.df$WT.AA == "G", 0,
ifelse(mut.df$WT.AA == "H", 20.44,
ifelse(mut.df$WT.AA == "I", 3.371,
ifelse(mut.df$WT.AA == "L", 3.782,
ifelse(mut.df$WT.AA == "K", 50.02,
ifelse(mut.df$WT.AA == "M", 8.589,
ifelse(mut.df$WT.AA == "F", 5.98,
ifelse(mut.df$WT.AA == "P", 7.916,
ifelse(mut.df$WT.AA == "S", 9.836,
ifelse(mut.df$WT.AA == "T", 9.304,
ifelse(mut.df$WT.AA == "W", 10.73,
ifelse(mut.df$WT.AA == "Y", 10.41,
print("There is at least one entry that isn't an amino acid")))))))))))))))))))))
mut.df$dipole <- as.numeric(mut.df$dipole)

mut.df$cp1 <- ifelse(mut.df$Mutated.AA1 == "V", 0.03115,
ifelse(mut.df$Mutated.AA1 == "A", 0.0196875,
ifelse(mut.df$Mutated.AA1 == "R", 0.0219875,
ifelse(mut.df$Mutated.AA1 == "N", 0.011375,
ifelse(mut.df$Mutated.AA1 == "D", 0.0126,
ifelse(mut.df$Mutated.AA1 == "C", 0.0118125,
ifelse(mut.df$Mutated.AA1 == "Q", 0.0134375,
ifelse(mut.df$Mutated.AA1 == "E", 0.015025,
ifelse(mut.df$Mutated.AA1 == "G", 0.0113125,
ifelse(mut.df$Mutated.AA1 == "H", 0.0257875,
ifelse(mut.df$Mutated.AA1 == "I", 0.037775,
ifelse(mut.df$Mutated.AA1 == "L", 0.0382875,
ifelse(mut.df$Mutated.AA1 == "K", 0.0288375,
ifelse(mut.df$Mutated.AA1 == "M", 0.030425,
ifelse(mut.df$Mutated.AA1 == "F", 0.0427125,
ifelse(mut.df$Mutated.AA1 == "P", 0.0310625,
ifelse(mut.df$Mutated.AA1 == "S", 0.014375,
ifelse(mut.df$Mutated.AA1 == "T", 0.0219875,
ifelse(mut.df$Mutated.AA1 == "W", 0.0453625,
ifelse(mut.df$Mutated.AA1 == "Y", 0.0352625,
print("There is at least one entry that isn't an amino acid")))))))))))))))))))))
mut.df$cp1 <- as.numeric(mut.df$cp1)
mut.df$cp1[is.na(mut.df$cp1)] <- 0
mut.df$cpabove <- ifelse(mut.df$cp > mean(mut.df$cp), 1, 0)
mut.df$size <- ifelse(mut.df$MW < 100, 0,
                      ifelse(mut.df$MW > 135, 1, 0.5))
mut.df$'Cp Change' <- ifelse(mut.df$Reported > 0, mut.df$cp1- mut.df$cp, 0)


mut.df$'Dipole Change' <- mut.df$dipole - 
                                    ifelse(mut.df$Mutated.AA1 == "V", 2.692,
                                    ifelse(mut.df$Mutated.AA1 == "A", 5.937,
                                    ifelse(mut.df$Mutated.AA1 == "R", 37.5,
                                    ifelse(mut.df$Mutated.AA1 == "N", 18.89,
                                    ifelse(mut.df$Mutated.AA1 == "D", 29.49,
                                    ifelse(mut.df$Mutated.AA1 == "C", 10.74,
                                    ifelse(mut.df$Mutated.AA1 == "Q", 39.89,
                                    ifelse(mut.df$Mutated.AA1 == "E", 42.52,
                                    ifelse(mut.df$Mutated.AA1 == "G", 0,
                                    ifelse(mut.df$Mutated.AA1 == "H", 20.44,
                                    ifelse(mut.df$Mutated.AA1 == "I", 3.371,
                                    ifelse(mut.df$Mutated.AA1 == "L", 3.782,
                                    ifelse(mut.df$Mutated.AA1 == "K", 50.02,
                                    ifelse(mut.df$Mutated.AA1 == "M", 8.589,
                                    ifelse(mut.df$Mutated.AA1 == "F", 5.98,
                                    ifelse(mut.df$Mutated.AA1 == "P", 7.916,
                                    ifelse(mut.df$Mutated.AA1 == "S", 9.836,
                                    ifelse(mut.df$Mutated.AA1 == "T", 9.304,
                                    ifelse(mut.df$Mutated.AA1 == "W", 10.73,
                                    ifelse(mut.df$Mutated.AA1 == "Y", 10.41, mut.df$dipole))))))))))))))))))))
mut.df$'MW Change' <- mut.df$MW - ifelse(mut.df$Mutated.AA1 == "G", 75.07,
ifelse(mut.df$Mutated.AA1 == "A", 89.09,
ifelse(mut.df$Mutated.AA1 == "S", 105.09,
ifelse(mut.df$Mutated.AA1 == "T", 119.1,
ifelse(mut.df$Mutated.AA1 == "C", 121.2,
ifelse(mut.df$Mutated.AA1 == "V", 117.1,
ifelse(mut.df$Mutated.AA1 == "L", 131.2,
ifelse(mut.df$Mutated.AA1 == "I", 131.2,
ifelse(mut.df$Mutated.AA1 == "M", 149.2,
ifelse(mut.df$Mutated.AA1 == "P", 115.1,
ifelse(mut.df$Mutated.AA1 == "F", 165.2,
ifelse(mut.df$Mutated.AA1 == "Y", 181.2,
ifelse(mut.df$Mutated.AA1 == "W", 204.2,
ifelse(mut.df$Mutated.AA1 == "D", 133.1,
ifelse(mut.df$Mutated.AA1 == "E", 147.1,
ifelse(mut.df$Mutated.AA1 == "N", 132.1,
ifelse(mut.df$Mutated.AA1 == "Q", 146.1,
ifelse(mut.df$Mutated.AA1 == "H", 155.2,
ifelse(mut.df$Mutated.AA1 == "K", 146.2, 174.2)))))))))))))))))))



desired.comparisons <- c(5, 21, 22, 23)
 
### Quick values
#Dystrophin
mut.df %>% 
  group_by(Effect) %>%
  summarize(sum(Reported), mean(Reported), sd(Reported), mean(cp), sd(cp), mean(MW), sd(MW))
library(ggcorrplot)
library(corrplot)
nonzeros <- which(mut.df$Reported !=0)
corrDys    <- cor(mut.df[nonzeros, desired.comparisons], use="pairwise.complete.obs")
corrDys95  <- cor_pmat(mut.df[nonzeros, desired.comparisons], conf.level=0.95)

################################### Dystrophin Correlogram #########################
ggcorrplot(corrDys, p.mat= corrDys95, type="lower", sig.level=0.05, hc.order = FALSE,
           colors=c("tomato2", "white", "springgreen3"), method="square",
           title="Correlogram for mutations in Dystrophin", lab=TRUE, 
           lab_size=4, ggtheme= theme_bw) +
  geom_rect(mapping=aes(xmin=0.5, xmax=3.5, ymin=0.5, ymax=1.5),
            color="red" , alpha=0, size=2) +
  scale_fill_gradient2("Correlation",low = "tomato2",midpoint=0 ,high="springgreen3", limits=c(-1,1))

correlation  <- cor(mut.df[, desired.comparisons], use="pairwise.complete.obs")
correlation95  <- cor_pmat(mut.df[, desired.comparisons], conf.level=0.95)
ggcorrplot(correlation, p.mat=correlation95, hc.order = FALSE, type="lower",
           lab=TRUE, lab_size=3, method="circle",
           colors=c("tomato2", "white", "springgreen3"),
           title="Correlogram for Dystrophin",
           ggtheme= theme_bw)

nonzeroplus <- which(mut.df$Reported != 0 & mut.df$Effect == "+")
corrDysplus  <- cor(mut.df[nonzeroplus, desired.comparisons], use="pairwise.complete.obs")
corrDysplus95  <- cor_pmat(mut.df[nonzeroplus, desired.comparisons], conf.level=0.95)
ggcorrplot(corrDysplus, p.mat=corrDysplus95, hc.order = FALSE, type="lower",
           lab=TRUE, lab_size=3, method="circle",
           colors=c("tomato2", "white", "springgreen3"),
           title="Correlogram for Disease causing mutations in Dystrophin",
           ggtheme= theme_bw)

nonzerominus <- which(mut.df$Reported !=0 & mut.df$Effect == "-")
corrDysminus  <- cor(mut.df[nonzerominus, desired.comparisons], use="pairwise.complete.obs")
corrDysminus95  <- cor_pmat(mut.df[nonzerominus, desired.comparisons], conf.level=0.95)
ggcorrplot(corrDysminus,p.mat=corrDysminus95, hc.order = FALSE, type="lower",
           lab=TRUE, lab_size=3, method="circle",
           colors=c("tomato2", "white", "springgreen3"),
           title="Correlogram for Asymptomatic mutations in Dystrophin",
           ggtheme= theme_bw)

diffdys <- abs(corrDysplus - corrDysminus)
ggcorrplot(diffdys, hc.order = FALSE, type="lower",
           lab=TRUE, lab_size=3, method="circle",
           colors=c("tomato2", "white", "springgreen3"),
           title="Correlogram for Difference between Asymptomatic and Symptomatic",
           ggtheme= theme_bw)




#ABD1
ABD1 <- which(mut.df$AA.pos < 247)
correlationABD1  <- cor(mut.df[ABD1, desired.comparisons], use="pairwise.complete.obs")
correlationABD195  <- cor_pmat(mut.df[ABD1, desired.comparisons], conf.level=0.95)
ggcorrplot(correlationABD1,p.mat=correlationABD195, hc.order = FALSE, type="lower",
           lab=TRUE, lab_size=3, method="circle",
           colors=c("tomato2", "white", "springgreen3"),
           title="Correlogram for ABD1",
           ggtheme= theme_bw)


nonzerosABD <- which(mut.df$Reported !=0 & mut.df$AA.pos < 247)
corrABD1  <- cor(mut.df[nonzerosABD, desired.comparisons], use="pairwise.complete.obs")
corrABD195 <- cor_pmat(mut.df[nonzerosABD, desired.comparisons], conf.level=0.95)
ggcorrplot(corrABD1,p.mat=corrABD195, hc.order = FALSE, type="lower",
           lab=TRUE, lab_size=3, method="circle",
           colors=c("tomato2", "white", "springgreen3"),
           title="Correlogram for ABD1 mutations",
           ggtheme= theme_bw) 



nonzeroplusABD <- which(mut.df$Reported != 0 & mut.df$Effect == "+" & mut.df$AA.pos < 246)
corrABDplus  <- cor(mut.df[nonzeroplusABD, desired.comparisons], use="pairwise.complete.obs")
corrABDplus95 <- cor_pmat(mut.df[nonzeroplusABD, desired.comparisons], conf.level=0.95)
ggcorrplot(corrABDplus, p.mat=corrABDplus95, hc.order = FALSE, type="lower",
           lab=TRUE, lab_size=3, method="circle",
           colors=c("tomato2", "white", "springgreen3"),
           title="Correlogram for Disease causing mutations in ABD1",
           ggtheme= theme_bw)

nonzerominusABD <- which(mut.df$Reported !=0 & mut.df$Effect == "-" & mut.df$AA.pos < 246)
corrABDminus  <- cor(mut.df[nonzerominusABD, desired.comparisons], use="pairwise.complete.obs")
corrABDminus95 <- cor_pmat(mut.df[nonzerominusABD, desired.comparisons], conf.level=0.95)
ggcorrplot(corrABDminus, p.mat=corrABDminus95, hc.order = FALSE, type="lower",
           lab=TRUE, lab_size=3, method="circle",
           colors=c("tomato2", "white", "springgreen3"),
           title="Correlogram of Asymptomatic mutations in ABD1",
           ggtheme= theme_bw)


##################################################### GOOD BELOW!!! ########
diffABD <- corrABDplus - corrABDminus
ggcorrplot(diffABD, hc.order = FALSE, type="lower",
           lab=TRUE, method="square",
           colors=c("tomato2", "white", "springgreen3"),
           title="Difference between Asymptomatic
and Symptomatic Mutations in ABD1",
           lab_size=4, ggtheme= theme_bw) +
  geom_rect(mapping=aes(xmin=0.5, xmax=1.5, ymin=0.5, ymax=1.5),
            color="red" , alpha=0, size=2) +
  scale_fill_gradient2("Correlation",low = "tomato2",midpoint=0 ,high="springgreen3", limits=c(-1,1))

#compare buried disease vs. non disease
buriedplus <- which(mut.df$Reported !=0 & mut.df$Effect == "+" & mut.df$AA.pos < 246 & mut.df$Water == "Buried")
corrburiedplus  <- cor(mut.df[buriedplus, desired.comparisons], use="pairwise.complete.obs")
corrburiedplus95 <- cor_pmat(mut.df[buriedplus, desired.comparisons], conf.level=0.95)

################################## Buried Correlogram ###########################
ggcorrplot(corrburiedplus, p.mat=corrburiedplus95, hc.order = FALSE, type="lower",
           lab=TRUE, method="square",
           title="Correlogram of Buried, Symptomatic Mutations in ABD1",
           lab_size=4, ggtheme= theme_bw) +
  geom_rect(mapping=aes(xmin=0.5, xmax=3.5, ymin=0.5, ymax=1.5),
            color="red" , alpha=0, size=2) +
  scale_fill_gradient2("Correlation",low = "tomato2",midpoint=0 ,high="springgreen3", limits=c(-1,1))


buriedminus <- which(mut.df$Reported !=0 & mut.df$Effect == "-" & mut.df$AA.pos < 246 & mut.df$Water == "Buried")
corrburiedminus  <- cor(mut.df[buriedminus, desired.comparisons], use="pairwise.complete.obs")
corrburiedminus95 <- cor_pmat(mut.df[buriedminus, desired.comparisons], conf.level=0.95)
ggcorrplot(corrburiedminus, p.mat=corrburiedminus95, hc.order = FALSE, type="lower",
           lab=TRUE, lab_size=3, method="circle",
           colors=c("tomato2", "white", "springgreen3"),
           title="Correlogram of Asymptomatic buried mutations in ABD1",
           ggtheme= theme_bw)
burieddif <- corrburiedplus - corrburiedminus
#################### Buried Difference Correlogram #########################
ggcorrplot(burieddif, hc.order = FALSE, type="lower",
           lab=TRUE, method="square",
           colors=c("tomato2", "white", "springgreen3"),
           title="Difference Between Buried Mutations in ABD1",
           lab_size=4, ggtheme= theme_bw) +
  geom_rect(mapping=aes(xmin=0.5, xmax=3.5, ymin=0.5, ymax=1.5),
            color="red" , alpha=0, size=2) +
  scale_fill_gradient2("Correlation",low = "tomato2",midpoint=0 ,high="springgreen3", limits=c(-1,1))



#compare exposed disease vs. non disease
expoplus <- which(mut.df$Reported !=0 & mut.df$Effect == "+" & mut.df$AA.pos < 246 & mut.df$Water == "Exposed")
correxpoplus  <- cor(mut.df[expoplus, desired.comparisons], use="pairwise.complete.obs")
correxpoplus95 <- cor_pmat(mut.df[expoplus, desired.comparisons], conf.level=0.95)
ggcorrplot(correxpoplus, p.mat=correxpoplus95, hc.order = FALSE, type="lower",
           lab=TRUE, lab_size=3, method="circle",
           colors=c("tomato2", "white", "springgreen3"),
           title="Correlogram of Exposed diseased mutations in ABD1",
           ggtheme= theme_bw)

expominus <- which(mut.df$Reported !=0 & mut.df$Effect == "-" & mut.df$AA.pos < 246 & mut.df$Water == "Exposed")
correxpominus  <- cor(mut.df[expominus, desired.comparisons], use="pairwise.complete.obs")
ggcorrplot(correxpominus, hc.order = FALSE, type="lower",
           lab=TRUE, lab_size=3, method="circle",
           colors=c("tomato2", "white", "springgreen3"),
           title="Correlogram of Asymptomatic exposed mutations in ABD1",
           ggtheme= theme_bw)
expodif <- correxpoplus - correxpominus
ggcorrplot(expodif, hc.order = FALSE, type="lower",
           lab=TRUE, lab_size=3, method="circle",
           colors=c("tomato2", "white", "springgreen3"),
           title="Difference between exposed +/- mutations in ABD1",
           ggtheme= theme_bw)


waterdif <- correxpoplus - corrburiedplus
ggcorrplot(waterdif, hc.order = FALSE, type="lower",
           lab=TRUE, lab_size=3, method="circle",
           colors=c("tomato2", "white", "springgreen3"),
           title="Difference between diseased exp/buried in ABD1",
           ggtheme= theme_bw)

# are ? same as + ???
var.test(mut.df$`change in cp`[(mut.df$Effect == "?" & mut.df$Water == "Exposed")], mut.df$`change in cp`[expoplus], alternative = "two.sided")
t.test(mut.df$`change in cp`[(mut.df$Effect == "?" & mut.df$Water == "Exposed")], mut.df$`change in cp`[expoplus], var.equal = TRUE) #different

var.test(mut.df$`change in cp`[(mut.df$Effect == "?" & mut.df$Water == "Buried")], mut.df$`change in cp`[expoplus], alternative = "two.sided")
t.test(mut.df$`change in cp`[(mut.df$Effect == "?" & mut.df$Water == "Buried")], mut.df$`change in cp`[expoplus], var.equal = TRUE) # different

# are ? same as - ?? 
var.test(mut.df$`change in cp`[(mut.df$Effect == "?" & mut.df$Water == "Exposed")], mut.df$`change in cp`[expominus], alternative = "two.sided") # different
t.test(mut.df$`change in cp`[(mut.df$Effect == "?" & mut.df$Water == "Exposed")], mut.df$`change in cp`[expominus], var.equal = FALSE) # different

var.test(mut.df$`change in cp`[(mut.df$Effect == "?" & mut.df$Water == "Buried")], mut.df$`change in cp`[buriedminus], alternative = "two.sided")
t.test(mut.df$`change in cp`[(mut.df$Effect == "?" & mut.df$Water == "Buried")], mut.df$`change in cp`[buriedminus], var.equal = TRUE) 

var.test(mut.df$`change in cp`[(mut.df$Effect == "?" & mut.df$Water == "Buried")], mut.df$`change in cp`[buriedplus], alternative = "two.sided")
t.test(mut.df$`change in cp`[(mut.df$Effect == "?" & mut.df$Water == "Buried")], mut.df$`change in cp`[buriedplus], var.equal = TRUE) 

var.test(mut.df$`change in cp`[(mut.df$Effect == "?" & mut.df$Water == "Buried")], mut.df$`change in cp`[expominus], alternative = "two.sided")
t.test(mut.df$`change in cp`[(mut.df$Effect == "?" & mut.df$Water == "Buried")], mut.df$`change in cp`[expominus], var.equal = FALSE) 

var.test(mut.df$`change in cp`[(mut.df$Effect == "?" & mut.df$Water == "Exposed")], mut.df$`change in cp`[buriedminus], alternative = "two.sided")
t.test(mut.df$`change in cp`[(mut.df$Effect == "?" & mut.df$Water == "Exposed")], mut.df$`change in cp`[buriedminus], var.equal = TRUE) 

var.test(mut.df$`change in cp`[(mut.df$Effect == "?" & mut.df$Water == "Buried")], mut.df$`change in cp`[expoplus], alternative = "two.sided")
t.test(mut.df$`change in cp`[(mut.df$Effect == "?" & mut.df$Water == "Buried")], mut.df$`change in cp`[expoplus], var.equal = TRUE) 

var.test(mut.df$`change in cp`[(mut.df$Effect == "?" & mut.df$Water == "Buried")], mut.df$`change in cp`[(mut.df$Effect == "?" & mut.df$Water == "Exposed")], alternative = "two.sided")
t.test(mut.df$`change in cp`[(mut.df$Effect == "?" & mut.df$Water == "Buried")], mut.df$`change in cp`[buriedminus], var.equal = TRUE) 

# change in cp
var.test(mut.df$`change in cp`[buriedminus], mut.df$`change in cp`[buriedplus], alternative = "two.sided")$p.value
t.test(mut.df$`change in cp`[buriedminus], mut.df$`change in cp`[buriedplus], var.equal = TRUE)$p.value
var.test(mut.df$`change in cp`[expoplus], mut.df$`change in cp`[buriedplus], alternative = "two.sided")$p.value
t.test(mut.df$`change in cp`[expoplus], mut.df$`change in cp`[buriedplus], var.equal = TRUE)$p.value
var.test(mut.df$`change in cp`[expominus], mut.df$`change in cp`[expoplus], alternative = "two.sided")$p.value
t.test(mut.df$`change in cp`[expominus], mut.df$`change in cp`[expoplus], var.equal = FALSE)$p.value
var.test(mut.df$`change in cp`[expominus], mut.df$`change in cp`[buriedminus], alternative = "two.sided")$p.value
t.test(mut.df$`change in cp`[expominus], mut.df$`change in cp`[buriedminus], var.equal = FALSE)$p.value
var.test(mut.df$`change in cp`[expominus], mut.df$`change in cp`[buriedplus], alternative = "two.sided")$p.value
t.test(mut.df$`change in cp`[expominus], mut.df$`change in cp`[buriedplus], var.equal = FALSE)$p.value
var.test(mut.df$`change in cp`[buriedminus], mut.df$`change in cp`[expoplus], alternative = "two.sided")$p.value
t.test(mut.df$`change in cp`[buriedminus], mut.df$`change in cp`[expoplus], var.equal = TRUE)$p.value



# WT cp
var.test(mut.df$cp[buriedminus], mut.df$cp[buriedplus], alternative = "two.sided")$p.value
var.test(mut.df$cp[expoplus], mut.df$cp[buriedplus], alternative = "two.sided")$p.value
var.test(mut.df$cp[expominus], mut.df$cp[expoplus], alternative = "two.sided")$p.value
var.test(mut.df$cp[expominus], mut.df$cp[buriedminus], alternative = "two.sided")$p.value
var.test(mut.df$cp[expominus], mut.df$cp[buriedplus], alternative = "two.sided")$p.value
var.test(mut.df$cp[buriedminus], mut.df$cp[expoplus], alternative = "two.sided")$p.value

t.test(mut.df$cp[buriedminus], mut.df$cp[buriedplus], var.equal = TRUE)$p.value
t.test(mut.df$cp[expoplus], mut.df$cp[buriedplus], var.equal = TRUE)$p.value
t.test(mut.df$cp[expominus], mut.df$cp[expoplus], var.equal = TRUE)$p.value
t.test(mut.df$cp[expominus], mut.df$cp[buriedminus], var.equal = TRUE)$p.value
t.test(mut.df$cp[expominus], mut.df$cp[buriedplus], var.equal = TRUE)$p.value
t.test(mut.df$cp[buriedminus], mut.df$cp[expoplus], var.equal = TRUE)$p.value


# WT dipole
var.test(mut.df$dipole[buriedminus], mut.df$dipole[buriedplus], alternative = "two.sided")$p.value
t.test(mut.df$dipole[buriedminus], mut.df$dipole[buriedplus], var.equal = TRUE)$p.value
var.test(mut.df$dipole[expoplus], mut.df$dipole[buriedplus], alternative = "two.sided")$p.value
t.test(mut.df$dipole[expoplus], mut.df$dipole[buriedplus], var.equal = FALSE)$p.value
var.test(mut.df$dipole[expominus], mut.df$dipole[expoplus], alternative = "two.sided")$p.value
t.test(mut.df$dipole[expominus], mut.df$dipole[expoplus], var.equal = TRUE)$p.value
var.test(mut.df$dipole[expominus], mut.df$dipole[buriedminus], alternative = "two.sided")$p.value
t.test(mut.df$dipole[expominus], mut.df$dipole[buriedminus], var.equal = TRUE)$p.value
var.test(mut.df$dipole[expominus], mut.df$dipole[buriedplus], alternative = "two.sided")$p.value
t.test(mut.df$dipole[expominus], mut.df$dipole[buriedplus], var.equal = TRUE)$p.value
var.test(mut.df$dipole[buriedminus], mut.df$dipole[expoplus], alternative = "two.sided")$p.value
t.test(mut.df$dipole[buriedminus], mut.df$dipole[expoplus], var.equal = TRUE)$p.value





# MW
var.test(mut.df$MW[buriedminus], mut.df$MW[buriedplus], alternative = "two.sided")$p.value
t.test(mut.df$MW[buriedminus], mut.df$MW[buriedplus], var.equal = TRUE)$p.value
var.test(mut.df$MW[expoplus], mut.df$MW[buriedplus], alternative = "two.sided")$p.value
t.test(mut.df$MW[expoplus], mut.df$MW[buriedplus], var.equal = FALSE)$p.value
var.test(mut.df$MW[expominus], mut.df$MW[expoplus], alternative = "two.sided")$p.value
t.test(mut.df$MW[expominus], mut.df$MW[expoplus], var.equal = TRUE)$p.value
var.test(mut.df$MW[expominus], mut.df$MW[buriedminus], alternative = "two.sided")$p.value
t.test(mut.df$MW[expominus], mut.df$MW[buriedminus], var.equal = TRUE)$p.value
var.test(mut.df$MW[expominus], mut.df$MW[buriedplus], alternative = "two.sided")$p.value
t.test(mut.df$MW[expominus], mut.df$MW[buriedplus], var.equal = TRUE)$p.value
var.test(mut.df$MW[buriedminus], mut.df$MW[expoplus], alternative = "two.sided")$p.value
t.test(mut.df$MW[buriedminus], mut.df$MW[expoplus], var.equal = TRUE)$p.value





### Compare ABD to Dys
ABD.DYS.muts <- abs(corrDys- corrABD1)
ggcorrplot(ABD.DYS.muts, hc.order = FALSE, type="lower",
           lab=TRUE, lab_size=3, method="circle",
           colors=c("tomato2", "white", "springgreen3"),
           title="Correlogram of Difference between ABD1 and Dys mutations",
           ggtheme= theme_bw)
# compare diseased mutations
ABD.DYS.dis.muts <- abs(corrDysplus - corrABDplus)
ggcorrplot(ABD.DYS.dis.muts, hc.order = FALSE, type="lower",
           lab=TRUE, lab_size=3, method="circle",
           colors=c("tomato2", "white", "springgreen3"),
           title="Correlogram of Difference disease causing mutations in ABD and Dys",
           ggtheme= theme_bw)

mut.df %>% 
  filter(AA.pos < 246) %>%
  group_by(Effect) %>%
  summarize(sum(Reported), mean(Reported), sd(Reported), mean(cp), sd(cp), mean(MW), sd(MW))


# Graph: number of reported mutations vs position colored by region
mut.df %>%
  filter(Reported >0) %>%
  ggplot(aes(x=AA.pos, y=Reported, col= Region, width=Reported)) +
  geom_line(alpha=0.5)
# Graph: number of reported mutations vs position colored by Domain
mut.df %>%
  filter(Reported >0) %>%
  ggplot(aes(x=AA.pos, y=Reported, fill=Domain)) +
  geom_area(alpha=0.5)

mut.df %>%
  filter(Water != "") %>%
  ggplot(aes(x=AA.pos, y=Water)) +
  geom_point()

mut.df %>%
  filter(Water != "") %>%
  ggplot(aes(x=Water, y=cp, fill=Water)) +
  geom_violin()

# Graph: The density of different amino acids that have been reportely mutated in a given region. So How many different amino acid
# Substitutions have been documented.
mut.df %>% 
  ggplot(aes(x=AA.pos, fill=Region)) +
  geom_density() +
  labs(x = "Position", y = "Mutations per Length")
################# original graph of mutations vs position#######
mut.df %>% 
  filter(Reported > 0 & Effect == "+") %>%
  ggplot(aes(x=AA.pos, y=Reported, col=Region)) +
  geom_point() +
  labs(x = "Position", y = "Symptomatic Mutations") +
  theme_bw() +
  theme(legend.position = "none") +
  geom_rect(xmin=-50, xmax=350, ymin=0.5, ymax=4.6, col='red', alpha=0, size=1) +
  geom_rect(xmin=3100, xmax=3600, ymin=0.5, ymax=3.6, col='red', alpha=0, size=1)
  

mut.df %>% 
  filter(Reported > 0 & Effect == "+") %>%
  ggplot(aes(x=AA.pos, fill=Region)) +
  geom_density(bw=25) +
  labs(x = "Position", y = "Mutations per Length")

mut.df %>% 
  filter(Reported > 0 & Effect == "+") %>%
  ggplot(aes(x=AA.pos, fill=Region, alpha = 0.2)) +
  geom_density(bw=25) +
  labs(x = "Position", y = "Mutations per Length") +
  scale_alpha_continuous(guide = FALSE)

library(raster)
mut.df %>% 
  filter(Reported > 1.5, Effect == "+") %>%
  ggplot(aes(x=AA.pos, y=Reported, col= Region, size=Reported)) +
  geom_point() +
  scale_size_continuous(guide=FALSE) +
  geom_rect(xmin=-50, xmax=320, ymin=1.5, ymax=5.1, col='red', alpha=0, size=1) +
  geom_rect(xmin=-50, xmax=320, ymin=1.5, ymax=5.1, col='red', alpha=0, size=1)
  
library(zoo)
diseases <- mut.df %>%
  filter(Effect == "+")
smooth <- rollapply(diseases$Reported, width = 4, mean , by = 4)
smooth.pos <- diseases$AA.pos[seq(1,length(diseases$AA.pos) - 3,
                                  by=4)]
smooth.reg <- diseases$Region[seq(1,length(diseases$AA.pos) - 3,
                                  by=4)]
df1 <- data.frame(smooth.pos, smooth, smooth.reg)


############################ This will be the original plot #################
ggplot(df1, aes(x=smooth.pos, y=smooth, col=smooth.reg)) +
  geom_point(alpha=0.7, size=3) +
  theme_bw() +
  scale_color_discrete(name = "Region") +
  labs(x= "Position", y="Number of Mutations") +
  geom_rect(xmin=-50,xmax=500, ymin=1.1, ymax=2.4,col='red',size=1, alpha=0) +
  geom_rect(xmin=3150, xmax=3500, ymin=0.9, ymax=1.6, col='red', size=1, alpha=0)
  

###### Average mutation density plot
Regions <- c("ABD1", "Undefined", "SYNM", "SNTB1",  "Zinc Finger")



zinlen <- length(unique(filter(mut.df, Region == "Zinc finger")$AA.pos))
ABDlen <- length(unique(filter(mut.df, Region == "ABD1")$AA.pos))
SNTB1len <- length(unique(filter(mut.df, Region == "SNTB1 binding")$AA.pos))
SYNMlen <- length(unique(filter(mut.df, Region == "SYNM binding")$AA.pos))
Undefinedlen <- length(unique(filter(mut.df, Region == "Undefined")$AA.pos))

zintot <-sum(filter(mut.df, Region == "Zinc finger")$Reported)
ABDtot <-sum(filter(mut.df, Region == "ABD1")$Reported)
SNTB1tot <-sum(filter(mut.df, Region == "SNTB1 binding")$Reported)
SYNMtot <-sum(filter(mut.df, Region == "SYNM binding")$Reported)
Undefinedtot <-sum(filter(mut.df, Region == "Undefined")$Reported)

zin <-zintot / zinlen
ABD <- ABDtot / ABDlen
SNTB1 <- SNTB1tot / SNTB1len
SYNM <- SYNMtot / SYNMlen
Undefined <- Undefinedtot / Undefinedlen


Density <- c(ABD, Undefined,  SYNM, SNTB1 ,   zin)
Cases    <- c(ABDtot, Undefinedtot, SYNMtot,SNTB1tot ,   zintot)
blah <- data.frame(Regions, Density, Cases)
ggplot(blah, aes(x=Regions, y= Density, col = Regions, size= Cases)) +
  geom_point() +
  labs(x= "Region",y = "Mutations per Length") +
  scale_x_discrete(limits=Regions) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.55)) + 
  geom_text(aes(label=Cases, size = 75),hjust=0.5, vjust=-1.5) +
  scale_y_continuous(limits = c(0, 0.33)) +
  scale_size_continuous(guide = FALSE)

#############################Disease causing ############
zintotdis <-sum(filter(mut.df, Region == "Zinc finger" & Effect == "+")$Reported)
ABDtotdis <-sum(filter(mut.df, Region == "ABD1" & Effect == "+")$Reported)
SNTB1totdis <-sum(filter(mut.df, Region == "SNTB1 binding" & Effect == "+")$Reported)
SYNMtotdis <-sum(filter(mut.df, Region == "SYNM binding" & Effect == "+")$Reported)
Undefinedtotdis <-sum(filter(mut.df, Region == "Undefined" & Effect == "+")$Reported)

zindis <-zintotdis / zinlen
ABDdis <- ABDtotdis / ABDlen
SNTB1dis <- SNTB1totdis / SNTB1len
SYNMdis <- SYNMtotdis / SYNMlen
Undefineddis <- Undefinedtotdis / Undefinedlen


Density <- c(ABDdis, Undefineddis, SYNMdis,  SNTB1dis, zindis)
Cases    <- c(ABDtotdis, Undefinedtotdis,SYNMtotdis,SNTB1totdis, zintotdis)
blah <- data.frame(Regions, Density, Cases)

###### Density plot comparison ########

ggplot(blah, aes(x=Regions, y= Density, col = c("#00B0F6","#E76BF3","#A3A500", "#00BF7D" ,"#F8766D"), size= Cases)) +
  geom_point() +
  labs(x= "Region",y = "Mutations per Amino Acid") +
  scale_x_discrete(limits=Regions) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.55)) + 
  geom_text(aes(label=Cases, size = 50),hjust=0.5, vjust=-1) +
  scale_y_continuous(limits = c(0, 0.25)) +
  scale_size_continuous(guide = FALSE) +
  theme_bw() + 
  theme(legend.position = "none", axis.text.x = element_text(size=10))+
  scale_x_discrete(limits=c("ABD1","Undefined","SYNM","Zinc Finger", "SNTB1")) +
  theme(axis.text.x = element_text(colour = c("#F8766D","#00B0F6","#00BF7D","#E76BF3","#A3A500")))
  geom_rect(mapping=aes(xmin=0.75, xmax=1.25, ymin=0.175, ymax=0.24),
            color="red" , alpha=0, size=2) +
  geom_rect(mapping=aes(xmin=3.75, xmax=4.25, ymin=0.2, ymax=0.24),
            color="red" , alpha=0, size=2)

Density <- c(0.19583333, 0.02871755, 0.01496259, 0.20833333, 0)
Cases    <- c(47, 73,12,10, 0)
blah <- data.frame(Regions, Density, Cases)

p.zin.dis <- zintotdis/ zintot
p.ABD.dis <- ABDtotdis/ ABDtot
p.SNTB1.dis <- SNTB1totdis/ SNTB1tot
p.SYNM.dis <- SYNMtotdis/ SYNMtot
p.Undefined.dis <- Undefinedtotdis/ Undefinedtot
percent.bad <- c(p.zin.dis, p.ABD.dis, p.SNTB1.dis, p.SYNM.dis, p.Undefined.dis)


mut.df %>% 
  ggplot(aes(x=AA.pos, fill=Domain)) +
  geom_density() +
  labs(x = "Position", y = "Mutations per Length")


mut.df %>%
  filter(Reported >0) %>%
  ggplot(aes(x=AA.pos, fill=Region)) +
  geom_density() +
  facet_wrap(~ Effect)

mut.df %>%
  filter(Reported >0) %>%
  ggplot(aes(x=AA.pos, y=Reported, col=Region)) +
  geom_point() +
  facet_wrap(~ Effect)


mut.df %>%
  filter(Reported > 0) %>%
  ggplot(aes(x=AA.pos, fill=Region)) +
  geom_histogram() +
  facet_wrap(~ Effect)

mut.df %>%
  filter(Reported > 0) %>%
  ggplot(aes(x=AA.pos)) +
  geom_histogram(style="fill") +
  facet_wrap(~ Effect)

mut.df %>%
  filter(Reported >0) %>%
  ggplot(aes(x=cp, y=Reported, col=Region )) +
  geom_line() +
  facet_wrap(~ Effect)

library("plotly")
#All dystrophin
mut.df%>%
  filter(Reported > 0) %>%
  plot_ly(x= ~cp, y=~Domain, z=~Reported, color=~Effect) 

# just ABD1
mut.df%>%
  filter(Reported > 0) %>%
  filter(Water != "") %>%
  plot_ly(x= ~cp, y=~Water, z=~Reported, color=~Effect)
#########################################################################SUPER IMPORTANT GRAPH BELOW ####
mut.df%>%
  filter(Reported > 0) %>%
  filter(Water != "") %>%
  plot_ly(x= ~`change in cp`, y=~Water, z=~Reported, color=~Effect, size=~Water, sizes=c(100, 500),
          colors=c("tomato2", "orange", "springgreen3"))
#########################################################################################################
mut.df %>%
  filter(Reported > 0) %>% 
  filter(AA.pos < 241) %>%
  plot_ly(x= ~cp, y=~Domain, z=~Reported, color=~Effect)

# Mutations vs MW
mut.df %>%
  filter(Reported >0) %>%
  ggplot(aes(x=MW)) +
  geom_area(aes(y=Reported, fill=Region)) +
  facet_wrap(~ Effect)
mut.df %>%
  filter(Reported >0) %>%
  ggplot(aes(x=Region, y=Reported, col=MW, size=3)) +
  geom_line() +
  facet_wrap(~ Effect)

# Mutations vs Cp
mut.df %>%
  filter(Reported >0) %>%
  ggplot(aes(x=cp)) +
  geom_area(aes(y=Reported, fill=Region)) +
  facet_wrap(~ Effect)
mut.df %>%
  filter(Reported >0) %>%
  ggplot(aes(x=Region, y=Reported, col=cp)) +
  geom_bar(stat="identity", width=.5)+
  facet_wrap(~ Effect)

  
  




ggplot(data= df.unique, aes(x=Position)) +
  geom_histogram(bins=10, fill=1, alpha=0.3)
#Broad strokes
ggplot(data= df.unique, aes(x=Position)) +
  geom_density(fill=1, alpha=0.3)
# finer separation
ggplot(data= df.unique, aes(x=Position)) +
  geom_density(fill=1, alpha=0.3, bw=70)
# position and effect
mut.df %>%
  filter(Reported > 0) %>%
  ggplot(aes(x=AA.pos, fill=Effect)) +
  geom_density(alpha=0.3, bw=250)
mut.df %>%
  filter(Reported > 0) %>%
  ggplot(aes(x=AA.pos, fill=as.factor(cpabove))) +
  geom_density(alpha=0.3, bw=100) +
  facet_wrap(~Effect)

mut.df %>%
  filter(Reported > 0) %>%
  ggplot(aes(x=AA.pos, fill=as.factor(polar))) +
  geom_density(alpha=0.3, bw=100) +
  facet_wrap(~Effect)


mut.df %>%
  filter(Reported > 0) %>%
  group_by(Effect) %>%
  summarize(mean(Reported),
            sum(Reported), length(Reported),
            sd(Reported))

ggplot(data=mut.df, aes(x=Reported, fill=Effect)) +
  geom_density(alpha=0.3)


ggplot(mut.df, aes(x= Index, y= act, col=Domain, size=act, cex.axis=2)) +
  geom_line() +
  xlab("Position") +
  ylab("number of mutations") +
  theme(legend.position="none", axis.text=element_text(size=20), axis.title=element_text(size=24,face="bold"))

ggplot(mut.df, aes(x= Index, y= allmuts, col=Domain, size=allmuts, cex.axis=2)) +
  geom_line() +
  xlab("Position") +
  ylab("number of mutations") +
  theme(legend.position="none", axis.text=element_text(size=20), axis.title=element_text(size=24,face="bold"))

ggplot(mut.df, aes(x= Index, y= allmuts, col=Region, size=allmuts)) +
  geom_line() +
  guides(size=FALSE) +
  xlab("Position") +
  ylab("Number of mutations") +
  theme(legend.position=c(0.1,0.8), legend.text=element_text(size=15),
        axis.text=element_text(size=20), axis.title=element_text(size=24,face="bold"))

ggplot(mut.df, aes(x= Region, y= allmuts, col=Domain, size=allmuts)) +
  geom_point() +
  ggtitle("Mutation positions") +
  xlab("Region") +
  ylab("number of mutations") +
  theme(legend.position="none")

ggplot(mut.df, aes(x= WT.AA, y= allmuts, col=Region, size=allmuts)) +
  geom_point() +
  theme(legend.position="none")

A <- sum(mut.df$Reported.1.[mut.df$WT.AA=="A"],mut.df$Reported.2.[mut.df$WT.AA=="A"],
         mut.df$Reported.3.[mut.df$WT.AA=="A"])
C <- sum(mut.df$Reported.1.[mut.df$WT.AA=="C"],mut.df$Reported.2.[mut.df$WT.AA=="C"],
         mut.df$Reported.3.[mut.df$WT.AA=="C"])
D <- sum(mut.df$Reported.1.[mut.df$WT.AA=="D"],mut.df$Reported.2.[mut.df$WT.AA=="D"],
         mut.df$Reported.3.[mut.df$WT.AA=="D"])
E <- sum(mut.df$Reported.1.[mut.df$WT.AA=="E"],mut.df$Reported.2.[mut.df$WT.AA=="E"],
         mut.df$Reported.3.[mut.df$WT.AA=="E"])
f <- sum(mut.df$Reported.1.[mut.df$WT.AA=="F"],mut.df$Reported.2.[mut.df$WT.AA=="F"],
         mut.df$Reported.3.[mut.df$WT.AA=="F"])
G <- sum(mut.df$Reported.1.[mut.df$WT.AA=="G"],mut.df$Reported.2.[mut.df$WT.AA=="G"],
         mut.df$Reported.3.[mut.df$WT.AA=="G"])
H <- sum(mut.df$Reported.1.[mut.df$WT.AA=="H"],mut.df$Reported.2.[mut.df$WT.AA=="H"],
         mut.df$Reported.3.[mut.df$WT.AA=="H"])
I <- sum(mut.df$Reported.1.[mut.df$WT.AA=="I"],mut.df$Reported.2.[mut.df$WT.AA=="I"],
         mut.df$Reported.3.[mut.df$WT.AA=="I"])
K <- sum(mut.df$Reported.1.[mut.df$WT.AA=="K"],mut.df$Reported.2.[mut.df$WT.AA=="K"],
         mut.df$Reported.3.[mut.df$WT.AA=="K"])
L <- sum(mut.df$Reported.1.[mut.df$WT.AA=="L"],mut.df$Reported.2.[mut.df$WT.AA=="L"],
         mut.df$Reported.3.[mut.df$WT.AA=="L"])
M <- sum(mut.df$Reported.1.[mut.df$WT.AA=="M"],mut.df$Reported.2.[mut.df$WT.AA=="M"],
         mut.df$Reported.3.[mut.df$WT.AA=="M"])
N <- sum(mut.df$Reported.1.[mut.df$WT.AA=="N"],mut.df$Reported.2.[mut.df$WT.AA=="N"],
         mut.df$Reported.3.[mut.df$WT.AA=="N"])
P <- sum(mut.df$Reported.1.[mut.df$WT.AA=="P"],mut.df$Reported.2.[mut.df$WT.AA=="P"],
         mut.df$Reported.3.[mut.df$WT.AA=="P"])
Q <- sum(mut.df$Reported.1.[mut.df$WT.AA=="Q"],mut.df$Reported.2.[mut.df$WT.AA=="Q"],
         mut.df$Reported.3.[mut.df$WT.AA=="Q"])
R <- sum(mut.df$Reported.1.[mut.df$WT.AA=="R"],mut.df$Reported.2.[mut.df$WT.AA=="R"],
         mut.df$Reported.3.[mut.df$WT.AA=="R"])
S <- sum(mut.df$Reported.1.[mut.df$WT.AA=="S"],mut.df$Reported.2.[mut.df$WT.AA=="S"],
         mut.df$Reported.3.[mut.df$WT.AA=="S"])
t <- sum(mut.df$Reported.1.[mut.df$WT.AA=="t"],mut.df$Reported.2.[mut.df$WT.AA=="t"],
         mut.df$Reported.3.[mut.df$WT.AA=="t"])
V <- sum(mut.df$Reported.1.[mut.df$WT.AA=="V"],mut.df$Reported.2.[mut.df$WT.AA=="V"],
         mut.df$Reported.3.[mut.df$WT.AA=="V"])
W <- sum(mut.df$Reported.1.[mut.df$WT.AA=="W"],mut.df$Reported.2.[mut.df$WT.AA=="W"],
         mut.df$Reported.3.[mut.df$WT.AA=="W"])
Y <- sum(mut.df$Reported.1.[mut.df$WT.AA=="Y"],mut.df$Reported.2.[mut.df$WT.AA=="Y"],
         mut.df$Reported.3.[mut.df$WT.AA=="Y"])

Adens <- A/sum(mut.df$WT.AA=="A")
Cdens <- C/sum(mut.df$WT.AA=="C")
Ddens <- D/sum(mut.df$WT.AA=="D")
Edens <- E/sum(mut.df$WT.AA=="E")
Fdens <- f/sum(mut.df$WT.AA=="F")
Gdens <- G/sum(mut.df$WT.AA=="G")
Hdens <- H/sum(mut.df$WT.AA=="H")
Idens <- I/sum(mut.df$WT.AA=="I")
Kdens <- K/sum(mut.df$WT.AA=="K")
Ldens <- L/sum(mut.df$WT.AA=="L")
Mdens <- M/sum(mut.df$WT.AA=="M")
Ndens <- N/sum(mut.df$WT.AA=="N")
Pdens <- P/sum(mut.df$WT.AA=="P")
Qdens <- Q/sum(mut.df$WT.AA=="Q")
Rdens <- R/sum(mut.df$WT.AA=="R")
Sdens <- S/sum(mut.df$WT.AA=="S")
Tdens <- t/sum(mut.df$WT.AA=="T")
Vdens <- V/sum(mut.df$WT.AA=="V")
Wdens <- W/sum(mut.df$WT.AA=="W")
Ydens <- Y/sum(mut.df$WT.AA=="Y")
den <- c(Adens, Cdens, Ddens, Edens, Fdens, Gdens, Hdens, Idens, Kdens, Ldens, Mdens, Ndens, Pdens, Qdens, Rdens,
         Sdens, Tdens, Vdens, Wdens, Ydens)
name <- c("A","C", "D", "E", "F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")


peram2 <- data.frame(den,name)
peram2 <- peram2[order(peram2$den),]
peram2$polar <- ifelse(  (peram2$name == "S") | (peram2$name == "T") |
                           (peram2$name == "C") | (peram2$name == "Y") |
                           (peram2$name == "N") | (peram2$name == "K") |
                           (peram2$name == "Q") | (peram2$name == "R") |
                           (peram2$name == "D") | (peram2$name == "H") |
                           (peram2$name == "E"), 1, 0)
peram$polar <- peram$polar+1
ggplot(peram, aes(x=name, y=den)) +
  geom_bar(stat="identity", fill=peram$polar) +
  theme(legend.position="none") +
  xlab("Amino Acid") +
  ylab("Mutations per residue") +
  scale_x_discrete(limits=peram$name)

ggplot(peram2, aes(x=name, y=den)) +
  geom_bar(stat="identity", fill=peram2$polar) +
  theme(legend.position="none") +
  xlab("Amino Acid") +
  ylab("Mutations per residue") +
  scale_x_discrete(limits=peram2$name)


peram <- peram[order(peram$name),]
peram2 <- peram2[order(peram2$name),]
peram$dif <- peram$den - peram2$den
peram$dif <- peram$dif[order(peram$dif)]
ggplot(peram, aes(x=name, y=dif)) +
  geom_bar(stat="identity", fill=peram2$polar) +
  theme(legend.position="none") +
  xlab("Amino Acid") +
  ylab("Mutations per residue") +
  scale_x_discrete(limits=peram$name)


unchart.regs         <- sum(mut.df$Region == "Undefined")
ABD.length           <- sum(mut.df$Region == "ABD")
SNTB1.binding.length <- sum(mut.df$Region == "SNTB1 binding")
SYNM.binding.length  <- sum(mut.df$Region == "SYNM binding")
zinc.finger.length   <- sum(mut.df$Region == "Zinc finger")

mutvsregion   <- table(mut.df$Region, mut.df$allmuts)
ABDmuts   <- sum(mutvsregion[1,-1])
SNTB1muts       <- sum(mutvsregion[2,-1])
SYNMmuts     <- sum(mutvsregion[3,-1])
unchartmuts      <- sum(mutvsregion[4,-1])
zincmuts      <- sum(mutvsregion[5,-1])
totalmuts     <- sum(unchartmuts, ABDmuts, SNTB1muts, SYNMmuts, zincmuts)


unchart.mut.p.length <- unchartmuts/unchart.regs
ABD.mut.p.length     <- ABDmuts/ABD.length
SNTB1.mut.p.length   <- SNTB1muts/SNTB1.binding.length
SYNM.mut.p.length    <- SYNMmuts/SYNM.binding.length
zinch.mut.p.length   <- zincmuts/zinc.finger.length
mut.p.length <- c(zinch.mut.p.length, ABD.mut.p.length, unchart.mut.p.length,  SYNM.mut.p.length, SNTB1.mut.p.length)

par(mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
plot(c(1,2,3,4,5),mut.p.length, type="n", axes=FALSE, ylim=c(0, 0.25))
points(c(1,2,3,4,5), mut.p.length, pch=c(16,16,15,16,1), 
       cex=c(zincmuts, ABDmuts, unchartmuts, SYNMmuts, 10)/10,
       col=c("black", "blue", "red", "black", "black"))
axis(2, col="grey", col.axis="black", at=seq(0,0.25, 0.05))
box(col="black")
mtext("Zinc finger", side = 1, cex = 1.5, line = 1, col = "black", adj=0)
mtext("ABD1",        side = 1, cex = 1.5, line = 1, col = "black", adj=0.25)
mtext("Undefined
Regions",   side = 1, cex = 1.5, line = 2, col = "black", adj=0.48)
mtext("SYNM",        side = 1, cex = 1.5, line = 1, col = "black", adj=0.76)
mtext("SNTB1",       side = 1, cex = 1.5, line = 1, col = "black", adj=1)
text(1.1, zinch.mut.p.length+0.005, labels="n=10", cex=2)
text(2, ABD.mut.p.length+0.007, labels="n=40", cex=2)
text(3, unchart.mut.p.length+0.015, labels="n=139", cex=2)
text(4, SYNM.mut.p.length+0.005, labels="n=36", cex=2)
text(4.9, SNTB1.mut.p.length+0.005, labels="n=0", cex=2)
mtext("Mutation density", side = 2, cex = 2, line = 2, col = "black")

par(opar)
colnames(mut.df) <- c("Region", "Domain", "Position", "WT AA", "AA1","AA2", "AA3",
                      "Cases 1", "Cases 2", "Cases 3","Surface or Buried","Polar", "charged", "MW", "Cp","Cp1", "Cp2", "Cp3", "Above average Cp?",
                      "total #", "AA size")
corrDys  <- cor(mut.df[nonzeros, c(12, 14, 15, 19, 20,22)])
corrplot.mixed(corr=corrDys, lower="number",
               upper="ellipse", upper.col=terrain.colors(10), lower.col="black", tl.col="blue", tl.cex=1.5)
mtext("Dystrophin correlations", side=3, line = 0, cex=2, adj=0.6)
ABD.df <- head(mut.df, n=241)
nonzerosABD <- which(ABD.df$`total #` !=0)
corrABD <- cor(ABD.df[nonzerosABD, c(12, 14, 15, 19, 20,22)])
corrplot.mixed(corr=corrABD, lower="number", 
               upper="ellipse", upper.col=terrain.colors(10), lower.col="black", tl.col="blue", tl.cex=1.5)
mtext("ABD1 correlations", side=3, line = 0, cex=2, adj=0.6)

corrdif <- corrDys - corrABD
corrplot.mixed(corr=corrdif, lower="number", 
               upper="ellipse", upper.col=terrain.colors(10), lower.col="black", tl.col="blue")
mtext("Difference", side=3, line = -8, cex=2, adj=0.6)

ABD.mut <- ABD.df[nonzerosABD,]
ABD.buried <- ABD.mut[ABD.mut$`Surface or Buried`=="Buried",]
ABD.surf   <- ABD.mut[ABD.mut$`Surface or Buried`== "Surface",]
t.test(ABD.surf$Cp, ABD.buried$Cp)
t.test(ABD.surf$`AA size`, ABD.buried$`AA size`)
t.test(ABD.surf$MW, ABD.buried$MW)
t.test(ABD.surf$Polar, ABD.buried$Polar)

mut.df$avedcp <- ((mut.df$Cp1*mut.df$`Cases 1` + mut.df$Cp2*mut.df$`Cases 2` +
                    mut.df$Cp3*mut.df$`Cases 3`)/(mut.df$`total #`))- mut.df$Cp
mut.df$avedcp[is.na(mut.df$avedcp)] <- 0 
corrSurf <- cor(ABD.surf[, c(12, 14, 15, 20,22)])
corrplot.mixed(corr=corrSurf, lower="number", 
               upper="ellipse", upper.col=terrain.colors(10), lower.col="black", tl.col="blue", tl.cex=1.5)
mtext("Surface Mutations", side=3, line = -4, cex=2, adj=0.6)
corrBury<- cor(ABD.buried[, c(12, 14, 15, 20,22)])
corrplot.mixed(corr=corrBury, lower="number", 
               upper="ellipse", upper.col=terrain.colors(10), lower.col="black", tl.col="blue", tl.cex=1.5)
mtext("Buried Mutations", side=3, line = -4, cex=2, adj=0.6)

bury.surf <- corrBury-corrSurf
corrplot.mixed(corr=bury.surf, lower="number", 
               upper="ellipse", upper.col=terrain.colors(10), lower.col="black", tl.col="blue", tl.cex=1.5)
mean(mut.df$avedcp)
sd(mut.df$avedcp)

mean(ABD.buried$avedcp)
sd(ABD.buried$avedcp)

mean(ABD.surf$avedcp)
sd(ABD.surf$avedcp)

t.test(ABD.surf$Cp, ABD.buried$Cp)
t.test(ABD.surf$avedcp, ABD.buried$avedcp)

dis.df <- read.csv(file.choose(), header = TRUE)

ggplot(dis.df, aes(x= N.aligned)) + 
  geom_point(aes(y = Mean.with.GS, color = "red")) + 
  geom_point(aes(y = Mean.without.GS, color = "blue")) +
  ggtitle("Disorder predictions") +
  xlab("Amino acid") +
  ylab("Disorder Index") +
  xlim(0, 50) +
  theme(legend.position=c(0.5,0.9)) + 
  labs(color = "Legend") +
  scale_color_manual(labels = c("Disorder without GS", "Disorder with GS"), values = c("blue", "red"))
