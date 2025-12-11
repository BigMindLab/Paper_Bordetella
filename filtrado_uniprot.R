
setwd("C:/Users/jalme/OneDrive/Escritorio/Immunoinformatics")

install.packages("writexl")
library(writexl)
library(openxlsx)
library(stringr)
library(dplyr)
table <- read.xlsx("uniprotkb_bordetella_AND_reviewed_true_2025_04_17 (2).xlsx")


## SIGNAL PEPTIDE: TRUE NEGATIVE
true_negative <- table[grepl("^SUBCELLULAR LOCATION: Cytoplasm",
                             table$`Subcellular.location.[CC]`), ]
signal_peptide_TN <- table[grepl("^SUBCELLULAR LOCATION: Cytoplasm",
                             table$`Subcellular.location.[CC]`) 
                       & grepl("^Bordetella", table$Organism), ]

write_xlsx(true_negative, "signal_peptide_TN.xlsx") #865 
write_xlsx(signal_peptide_TN, "SP_bortedella_TN.xlsx") #863

filtrado <- true_negative[!grepl("Secreted", true_negative$`Subcellular.location.[CC]`
                                 , ignore.case = TRUE), ]

## SIGNAL PEPTIDE: TRUE POSITIVE
signal_peptide <- table[grepl("^SIGNAL", table$Signal.peptide) & 
                          !is.na(table$`Subcellular.location.[CC]`), ]

signal_peptide_1 <- table[grepl("^SIGNAL", table$Signal.peptide) & 
                          !is.na(table$`Subcellular.location.[CC]`) &
                          grepl("^Bordetella", table$Organism), ] #83



write_xlsx(signal_peptide, "signal_peptide_TP.xlsx") #91 
write_xlsx(signal_peptide_1, "SP_bortedella_TP.xlsx") #83


## TRANSMEMBRANE : TRUE POSITIVE
transmembrane_TP <- table[grepl("^TRANSMEM",
                             table$Transmembrane), ]
transmembrane_TP_1 <- table[grepl("^TRANSMEM",
                                table$Transmembrane) &
                              grepl("^Bordetella", table$Organism), ] #149

transmembrane_TP_1$UniProt <- str_count(transmembrane_TP_1$Transmembrane, "TRANSMEM")
write_xlsx(transmembrane_TP_1, "transmembrane_UniProt_TP.xlsx") #149

write_xlsx(transmembrane_TP, "transmembrane_TP.xlsx") #155 
write_xlsx(transmembrane_TP_1, "transmembrane_bortedella_TP.xlsx")

## TRANSMEMBRANE : TRUE NEGATIVE
transmembrane_TN <- table[!grepl("^TRANSMEM", 
                                 table$Transmembrane, ignore.case = TRUE) 
                          & !is.na(table$`Subcellular.location.[CC]`), ]

transmembrane_TN_1 <- table[!grepl("^TRANSMEM", 
                                 table$Transmembrane, ignore.case = TRUE) 
                          & !is.na(table$`Subcellular.location.[CC]`) 
                          & !grepl("membrane", table$`Subcellular.location.[CC]`, 
                                   ignore.case = TRUE), ]
transmembrane_TN_2 <- table[!grepl("^TRANSMEM", 
                                   table$Transmembrane, ignore.case = TRUE) 
                            & !is.na(table$`Subcellular.location.[CC]`) 
                            & !grepl("membrane", table$`Subcellular.location.[CC]`, 
                                     ignore.case = TRUE)& grepl("^Bordetella", table$Organism), ]

write_xlsx(transmembrane_TN, "transmembrane_TN.xlsx") #901
write_xlsx(transmembrane_TN_2, "transmembrane_bortedella_TN.xlsx")#894

nrow(transmembrane_TP)



filtro <- true_negative[true_negative$Length >= 200, ]

filtrado <- true_negative[!grepl("Secreted", true_negative$`Subcellular.location.[CC]`
                                 , ignore.case = TRUE), ]

write_xlsx(true_negative, "tabla.xlsx")


library(biomaRt)
ensembl <- useMart("ensembl")

listDatasets(ensembl)
