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
samplecount <- as.data.frame(read.csv('methylresolver_data.txt', header = TRUE, sep = '\t',row.names =1))


# Deconvolution with default signature and calculating absolute fractions with default RF model:
#MethylResolver(methylMix = samplecount)

# Deconvolution with default signature and calculating absolute fractions with default RF model while 
# specifying a particular alpha value:
MethylResolver(methylMix = samplecount, alpha = 0.5)

# Deconvolution with default signature and only calculating relative fractions:
#MethylResolver(methylMix = samplecount, absolute = FALSE)

# Specify your own signature matrix and RF model for calculating absolute fractions:
#MethylResolver(methylMix = samplecount, methylSig = MethylSig, purityModel = RFmodel)
