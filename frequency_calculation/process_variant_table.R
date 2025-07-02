### Add id names based on position information
# We received a CSV table from the server with the header: c("CHROM", "Position", "DP", "DP4", "REF", "ALT", "BamID")
modern = comperative_alkan

## Combine chromosome and position information
modern$POS = paste(modern$CHROM, modern$POS, sep = "_")

## Extract the id information from the Risk_allele table into modern
colnames(Risk_allele)[9] = "POS"
# Risk_allele$POS = paste("chr", Risk_allele$Chromosome, "_" ,Risk_allele$POS, sep = "" )
# snip_id = Risk_allele[,c(9,1)]

modern$SNP_ID = NA
for (i in 1:nrow(modern)) {
  modern[i,8] = snip_id[which(modern[i,2] == snip_id[,1]),2]    # Compares the position column of modern with the position column of snip_id and, if they match, retrieves the id from the second column of snip_id
}

## Add risk allele and minor allele information
modern$Risk_Allele = NA
modern$Minor_Allele = NA

for (i in 1:nrow(modern)) {
  modern[i,9] = Risk_allele[which(modern[i,2] == Risk_allele[,9]),3] 
  modern[i,10] = Risk_allele[which(modern[i,2] == Risk_allele[,9]),4] 
}

## Split the DP4 column into 4 separate columns
library(stringr)
modern_DP4 = as.data.frame(str_split_fixed(modern$DP4, ",", 4)); head(modern_DP4)

modern = as.data.frame(cbind(modern$POS, modern$DP, modern_DP4, 
                             modern$REF, modern$ALT, modern$name, modern$SNP_ID,
                             modern$Risk_Allele, modern$Minor_Allele))
colnames(modern) = c("POS", "DP", "Ref1","Ref2","Alt1","Alt2","REF","ALT","BAMID", "SNP_ID", "Risk_Allele", "Minor_Allele")

## Add phenotype information
modern$Phenotype = NA

for (i in 1:nrow(modern)) {
  modern[i,13] = Risk_allele[which(modern[i,1] == Risk_allele[,9]), 7]
}

### Risk allele, Total_RAC column
modern$Ref1 = as.numeric(modern$Ref1)
modern$Ref2 = as.numeric(modern$Ref2)
modern$Alt1 = as.numeric(modern$Alt1)
modern$Alt2 = as.numeric(modern$Alt2)

modern$Total_RAC = ifelse(modern$Risk_Allele == modern$REF,
                          modern$Ref1 + modern$Ref2,
                          modern$Alt1 + modern$Alt2)

## For neutral alleles where Risk_Allele is NA, perform the comparison using Minor_Allele
for (i in 1:nrow(modern)) {
  if(modern$Risk_Allele[i] == "NA"){
    if(modern$Minor_Allele[i] == modern$REF[i]){
      modern$Total_RAC[i] = modern$Ref1[i] + modern$Ref2[i]
    } else {
      modern$Total_RAC[i] = modern$Alt1[i] + modern$Alt2[i]
    }
  }
}

### Add group information
modern$Group = "Modern"

modern_df = modern

### Final form
modern_df = as.data.frame(cbind(modern_df$BAMID, modern_df$Group, modern_df$POS, modern_df$SNP_ID, modern_df$Phenotype, modern_df$Total_RAC, modern_df$DP))
colnames(modern_df) = c("SampleName", "Group", "POS", "SNP_ID", "Phenotype", "R", "T")
