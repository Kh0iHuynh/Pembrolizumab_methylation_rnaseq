library(tidyverse)
library(MethylResolver)
#####
# turn of scientific notation
#####
options(scipen = 999)

####
# read in normalized count
# convert SRR col name to row values
####
metadata <- read.table("deseq_metadata_withouthyphen.csv",header = TRUE)

samplecount <- as.data.frame(read.csv('../methylresolver_data.txt', header = TRUE, sep = '\t',row.names =1))

# Subset df1 to keep only columns matching df2$Sample_ID
samplecount_filtered <- samplecount[, colnames(samplecount) %in% metadata$Sample_ID]

# Deconvolution with default signature and calculating absolute fractions with default RF model:
#MethylResolver(methylMix = samplecount)

# Deconvolution with default signature and calculating absolute fractions with default RF model while 
# specifying a particular alpha value:
#MethylResolver(methylMix = samplecount, alpha = 0.5)
MethylResolver(methylMix = samplecount_filtered, alpha = 0.5, absolute = FALSE)
# Deconvolution with default signature and only calculating relative fractions:
#MethylResolver(methylMix = samplecount, absolute = FALSE)

# Specify your own signature matrix and RF model for calculating absolute fractions:
#MethylResolver(methylMix = samplecount, methylSig = MethylSig, purityModel = RFmodel)
