## set up analysis environment
library(zip)
library(RnBeads.hg19)
library(RnBeads.hg38)
library(RnBeads)
# Directory where your data is located
data.dir <- "/project/modrek_1129/KhoiHuynh/methylation/Methylation_Array_Data/Organized_IDAT_Files/chow"
idat.dir <- file.path(data.dir, "idat")
print(idat.dir)
# Directory where the output should be written to
analysis.dir <- "./"
# Directory where the report files should be written to
report.dir <- file.path(analysis.dir, "reports")
sample.annotation <- file.path(analysis.dir, "sample_annotation_test.csv")
#get hg38 annotation


#set rnbead option

rnb.options(
	identifiers.column                = "Sample_ID",
	assembly                          = "hg19",
	export.to.bed = TRUE,
	export.to.trackhub = "bigWig",
	export.to.csv = TRUE,
	export.to.ewasher = TRUE,
	filtering.low.coverage.masking    = TRUE,
	filtering.greedycut               = FALSE,
	filtering.missing.value.quantile  = 0.5,
	filtering.high.coverage.outliers  = TRUE,
	inference=TRUE,
	differential.comparison.columns   = c("Sex","Age","Overall_Survival","RANO","Relapses_Count","OS_Status","IFN-G_response","KPS")
	#differential.site.test.method="refFreeEWAS"
)
# optionally disable some parts of the analysis to reduce runtime
# optionally disable some parts of the analysis to reduce runtime
rnb.options(
	exploratory.intersample           = FALSE,
	exploratory.region.profiles       = character(0),
	differential.report.sites         = FALSE
)
#rnb.run.analysis(dir.reports=report.dir, sample.sheet=sample.annotation, data.dir=idat.dir,data.type="infinium.idat.dir")
rnb.run.analysis(
	dir.reports=report.dir,
	sample.sheet=sample.annotation,
	data.dir=idat.dir,
	data.type="idat.dir"
)
