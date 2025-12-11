

setwd("C:/Users/jalme/OneDrive/Escritorio/Immunoinformatics")

##BIOCONDUCTOR
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

##Biostrings
library(Biostrings)

##verificar la longitud de las secuencias
cds <- readDNAStringSet("selected_genes.faa")
table(width(cds) %% 3)

##verificar que al traducir aparezca solo un codón stop al final
aa <- translate(cds, genetic.code = getGeneticCode("11"))
table(grepl("\\*", as.character(aa)))


# Buscar letras que no sean ACGTN
bad_idx <- vcountPattern("[^ACGTN]", cds, fixed=FALSE) > 0
table(bad_idx)

#código genético de bacteria
gc11 <- getGeneticCode("11")


##ejemplo 
translate(DNAString("acgtncagc"), genetic.code=gc11,
          if.fuzzy.codon="X", no.init.codon=TRUE)

# Marca codones ambiguos como X (seguro)
aa <- translate(cds, genetic.code = gc11,
                if.fuzzy.codon = "X", # cualquier codón con N -> 'X'
                no.init.codon = TRUE)

# Excluir el último aa (donde debe ir el stop final)
has_stop11 <- vcountPattern("*", subseq(aa, 1, width(aa)-1), fixed=TRUE) > 0  
sum(has_stop11)


length(aa)                # cuántas proteínas
names(aa)[1:10]           # primeros IDs
aa[1:3]                   # primer/ as 3 como AAStringSet (muestra resumen)
as.character(aa[1:3])     # las secuencias como texto

##buscar y ver solo group:
aa_group <- aa[ grepl("^group_", names(aa)) ]
length(aa_group)
as.character(aa_group[1:3])


##al menos un X
# aa es tu AAStringSet
x_per_seq <- vcountPattern("X", aa)   # nº de X por secuencia
sum(x_per_seq > 0)                    # cuántas tienen algún X
mean(x_per_seq > 0) * 100             # % con X
summary(x_per_seq)                    # distribución de X por secuencia

#tipo de objeto
class(aa)

#nombre de las secuencias
names(aa)



##Descargar archivo fasta
writeXStringSet(aa, filepath = "proteinas.faa", format = "fasta", width = 60)




