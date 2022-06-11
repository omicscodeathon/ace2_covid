install.packages("enrichR")

# By default human genes are selected otherwise select your organism of choice.
library(enrichR)
listEnrichrSites()
setEnrichrSite("Enrichr") # Human genes
websiteLive <- TRUE

#  find the list of all available databases from Enrichr.
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)

library(knitr)
if (websiteLive) kable(head(dbs[c(1:6),-4]))

# query enrichr for downregulated genes
dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021")
if (websiteLive) {
  enriched <- enrichr(c("ZNF460", "RNR2", "ZBED6", 
                        "RNR1", "SNORA70G", "RPL19P12", "NACA4P", "RPL23P8", "FOSB", "NACA2", "PRNCR1", 
                        "FAAH2", "SNORD141A", "LOC105374010", "PCNP", "LOC102723752", "LOC105374216", 
                        "LOC105369180", "LOC107986254", "TSPYL4", "RPS18P9", "CAMK4", "TULP4", "FSIP2", 
                        "MIR4426", "LOC105371085", "SHPRH", "INPP4B", "LOC101927840", "SNORA3A", "POMK",
                        "HMCN1", "SNORA3B", "SPATA5", "PARD6B", "RFX3", "ELMO1-AS1", "SNORD12C", "UBR1", 
                        "ARID5B", "CCDC18-AS1", "MED10", "LOC105377148", "LOC105379752", "H3C1", 
                        "IMMP2L", "LOC107987043", "MRS2P2", "NCBP2AS2", "LOC105374724", "NLGN1", "KLHL28", 
                        "MAPKAPK5-AS1", "PHC3", "PHACTR2-AS1", "LRP1B", "SNHG1", "ETNK1", "MAP1B", 
                        "PKHD1", "RAD54L2", "TRNT", "LOC374443", "RPS13", "RO60", "LOC107984315", "MUC19", 
                        "FAM172A", "LOC114224", "SLC7A6", "LATS1", "METTL17", "PHF11", "CD69", "FOS", 
                        "LOC105373269", "FOXP1", "NBEAL1", "ARHGAP15-AS1", "LOC105376364", "LOC102724146", "CDH8", 
                        "MGAT4C","UXT", "TTC19", "MGC4859", "GAS5", "LOC105372118", "LOC105379524", "ULK4", 
                        "HELZ", "DARS1", "EIF4A2", "SNHG15", "RIMOC1", "NKAIN3", "PIK3C2G", "SBDS", "LOC101927817"), dbs)
}

# view the results table (Mol. Funct.)
if (websiteLive) enriched[["GO_Molecular_Function_2021"]]

if (websiteLive) {
  x <- head(enriched[["GO_Molecular_Function_2021"]])
  x[,1] <- gsub("GO:", "GO_", x[,1])
  kable(x)
}

# Plot Enrichr GO-MF output
if (websiteLive) plotEnrich(enriched[[1]], showTerms =15, numChar = 30, 
                            y = "Count", orderBy = "P.value", 
                            title = "GO Molecular Function"
)

# view the results table (Cell comp..)
if (websiteLive) enriched[["GO_Cellular_Component_2021"]]

if (websiteLive) {
  x <- head(enriched[["GO_Cellular_Component_2021"]])
  x[,1] <- gsub("GO:", "GO_", x[,1])
  kable(x)
}

# Plot Enrichr GO-CC output
if (websiteLive) plotEnrich(enriched[[2]], showTerms = 15, numChar = 30, 
                            y = "Count", orderBy = "P.value", 
                            title = "GO Cellular Component"
)

# view the results table (Bio process)
if (websiteLive) enriched[["GO_Biological_Process_2021"]]

if (websiteLive) {
  x <- head(enriched[["GO_Biological_Process_2021"]])
  x[,1] <- gsub("GO:", "GO_", x[,1])
  kable(x)
}

# Plot Enrichr GO-BP output
if (websiteLive) plotEnrich(enriched[[3]], showTerms = 15, numChar = 30, 
                            y = "Count", orderBy = "P.value", 
                            title = "GO Biological Process"
)

