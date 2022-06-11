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
  enriched <- enrichr(c("RNR2", "TRNT", "LOC105376364", "CDH8", "SNORD27", "MDGA2", "LOC105370777", "ANO3",
                        "SNX9-AS1", "TEX15", "OTOGL", "GABRB2", "LOC105379080", "LOC105375153", "FSIP2", "LOC101929278",
                        "PTPRQ", "PKHD1", "ADGRL3", "OR4K2", "KCNH5", "CSRNP3", "LOC105378798", "LOC107987107", "ASB15",
                        "LOC105370236", "FGF12", "LINC02600", "LOC105377567", "DLG2-AS2",  "RNR1", "ST7-OT4", "SI", "LOC112268150",
                        "OR56A1", "LOC105375149", "LOC105374235", "HMCN2", "LOC105373703", "LOC105375097", "SYT14", "CNTN6",
                        "LOC112268125", "LOC105373648", "GABRA2", "PPP1R2C", "LOC105378747", "LINC02822", "LOC105376960", "LOC105373776",
                        "GLI2", "LOC105377445", "LOC105378340", "LOC105370781", "LOC105374473", "LOC105375997", "LOC105377992",
                        "MAPK4", "MROH9", "CNTN5", "TSG1", "MIR4485", "KCNMB2-AS1", "LOC105372031", "CNTN3", "PKIA-AS1", "GABRG1",
                        "CLVS2", "LRRIQ1", "LOC107986020", "LOC101928849", "STXBP5L", "MPPED2-AS1", "LRRC3B-AS1", "MIR181A1HG", "LOC102724146",
                        "MIR623", "TENT5C-DT", "LINC02265", "CACNA2D1-AS1", "RGS1", "LOC105377280", "OR2M3", "LOC105376228", "OR5AS1",
                        "C8orf34", "MICOS10-DT", "IGFBP-AS1", "LOC105377017", "TRDN", "LOC112267875", "LOC107984894", "SYNJ2-IT1",
                        "LINC01794", "LOC107984142", "SH3GL2", "LOC105369456", "CFAP69", "LOC105377977", "LOC112268015"), dbs)
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

