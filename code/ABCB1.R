library(tidyverse)
library(haplo.stats)

phenodata <- read_csv("abcb1_haplotyoes.csv")
phenodata <- phenodata %>% select(ID, label)

data_reshape <- read_csv("abcb1_haplotyoes.csv")

reshape <- data_reshape %>% 
  separate(rs1128503, c("rs1128503_1"  , "rs1128503_2")) %>% 
  separate(rs2032582, c("rs2032582_1"   , "rs2032582_2")) %>% 
  separate(rs1045642, c("rs1045642_1"  , "rs1045642_2")) 

reshape <- reshape[,3:8]

reshape <- as.data.frame(lapply(reshape, unlist))

label <-c("rs1128503","rs2032582","rs1045642")

set.seed(123)

save.em <- haplo.em(geno = reshape, locus.label = label)

save.em$haplotype

summary <- summary(save.em)

summary <- summary %>% filter(posterior > 0.8)

summary_final <- summary %>% 
  mutate(hap1 = recode(hap1,
                       `1`  = "C/G/C",
                       `2`  = "C/G/T",
                       `5`  = "T/G/C",
                       `8`  = "T/nonG/T"),
         hap2 = recode(hap2,
                       `1`  = "C/G/C",
                       `2`  = "C/G/T",
                       `3`  = "C/nonG/C",
                       `7`  = "T/nonG/C",
                       `8` =  "T/nonG/T"))

final <- merge(summary_final, phenodata, all.y = TRUE, by.x = "subj.id", by.y = "label")

final <- final %>% select(ID, hap1, hap2, posterior)

write_csv(final, "output_ABCB1_v1.csv")

final_diplotype <- final %>% unite("diplotype", hap1:hap2, sep = "_", remove = FALSE)
