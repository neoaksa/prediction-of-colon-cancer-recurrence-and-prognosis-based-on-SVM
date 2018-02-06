#2.1 MicroArrayData
source("https://bioconductor.org/biocLite.R")
biocLite("hgu133plus2.db")
library(hgu133plus2.db) 
biocLite("hgu133plus2cdf")
library(hgu133plus2cdf) 

library(affy)
library(limma)


library("GSEABase")
library(GEOquery)

# store all CEL file names in the variable fn
setwd("//fileserv/Profiles$/tweedya/Desktop/project 1/data/GSE17537_RAW")
fn<-dir()  
str(fn) 

rawAffyData <- ReadAffy(filenames=fn, cdfname="hgu133plus2cdf")  

#2.2 Data Normalization
# to use the justRMA function we need to make a character vector of CEL files 
eset <- rma(rawAffyData)

# extract the data matrix from the expression set
exprs <- exprs(eset) 

# get the gene symbols   
# The platform is hgu133a. This links the probe ids with the gene names.
# List all available options in the package
ls("package:hgu133plus2.db")

# Replace the probe ids with gene symbols (3 steps)
# 1. Assign the row names (probe ids) of expression matrix to variable sid 
sid<-rownames(exprs)

# 2. extract the gene names
sym <- unlist(mget(sid,hgu133plus2SYMBOL,ifnotfound=NA))  

# 3. replace with gene symbols
rownames(exprs)<-sym
colnames(exprs)<-  c("NoRecur","NoRecur","NoRecur","NoRecur","Recur","NoRecur","Recur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","Recur","NoRecur","NoRecur","NoRecur","Recur","NoRecur","Recur","Recur","NoRecur","Recur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","Recur","Recur","NoRecur","NoRecur","Recur","NoRecur","NoRecur","Recur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","Recur","NoRecur","Recur","Recur","Recur","Recur","Recur","Recur","Recur","NoRecur","NoRecur")

# Assign column names into the variable levels 
levels<-colnames(exprs)  

# Create design matrix for model 
design<-model.matrix(~0+factor(levels))  

# Changing the column names of design 
colnames(design)<-c("NoRecur","Recur") 

# Fitting a linear model to the expression matrix as per the design matrix 
fit <- lmFit(exprs, design) 

# Creating a contrast matrix as contrast between the two classes NoRecur and Recur
contrast.matrix <-makeContrasts(NoRecur-Recur ,levels=design)

# Applying the contrast.fit 
fit <- contrasts.fit(fit, contrast.matrix) 

# eBayes statistic 
fit <- eBayes(fit) 

# find the top genes 
siggenes <- topTable(fit, coef=1, adjust= "none", sort.by= "logFC",number=60000)
#head(siggenes)
#not sure this is the right filter criteria. 
logsel = (abs(siggenes$logFC) < .5 | abs(siggenes$logFC) > 2)

# find top genes using adjusted p-value .05
selected  <- p.adjust(fit$p.value,method="none") <0.05
esetSel <- eset [(logsel & selected), ]
dim(esetSel)

#2.4. Protein-protein interaction (PPI) network construction

#2.5. Network-based neighborhood scoring analysis
