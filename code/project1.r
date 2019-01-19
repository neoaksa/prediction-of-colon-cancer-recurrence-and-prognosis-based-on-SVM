#2.1 MicroArrayData
source("https://bioconductor.org/biocLite.R")
# biocLite("hgu133plus2.db")
library(hgu133plus2.db) 
# biocLite("hgu133plus2cdf")
library(hgu133plus2cdf) 

library(affy)
library(limma)

library("GSEABase")
library(GEOquery)

# store all CEL file names in the variable fn
setwd("//fileserv/Profiles$/tweedya/Desktop/project 1/data/GSE17537_RAW")
#setwd("//Users/jietao/Google Drive/gvsu/course/CIS635/projects1/data/GSE17537_RAW")
fn<-dir()  
str(fn) 

dim(rawAffyData)
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

setwd("//fileserv/Profiles$/tweedya/Desktop/project 1/data/")
gse <- getGEO("GSE17537", destdir = getwd())
recurrance_data1 <- pData(gse[[1]])$characteristics_ch1.5
recurrance_data2 <- pData(gse[[1]])$characteristics_ch1.6

no_recur1 = grep("(: no recurrence)",recurrance_data1)
recur1 = grep("(: recurrence)",recurrance_data1)

no_recur2 = grep("(: no recurrence)",recurrance_data2)
recur2 = grep("(: recurrence)",recurrance_data2)

colnames(exprs)[no_recur1] <- "NoRecur"
colnames(exprs)[no_recur2] <- "NoRecur"
colnames(exprs)[recur1] <- "Recur"
colnames(exprs)[recur2] <- "Recur"
summary(as.factor(colnames(exprs)))

#colnames(exprs)<-  c("NoRecur","NoRecur","NoRecur","NoRecur","Recur","NoRecur","Recur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","Recur","NoRecur","NoRecur","NoRecur","Recur","NoRecur","Recur","Recur","NoRecur","Recur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","Recur","Recur","NoRecur","NoRecur","Recur","NoRecur","NoRecur","Recur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","Recur","NoRecur","Recur","Recur","Recur","Recur","Recur","Recur","Recur","NoRecur","NoRecur")

# Assign column names into the variable levels 
levels<-colnames(exprs)  

# Create design matrix for model 
design<-model.matrix(~0+factor(levels))  
?model.matrix
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
siggenes <- topTable(fit, coef=1, adjust= "none", sort.by= "none",number=60000)
head(siggenes)
logsel = (abs(siggenes$logFC) > 0.7 & siggenes$adj.P.Val<0.05)

esetSel <- eset [(logsel ), ]
dim(esetSel)
head(exprs(esetSel))
# heatmap 
heatmap(exprs(esetSel))
# volcanopot
volcanoplot(fit, coef=1, highlight=10, names=fit$genes$ID,
            xlab="Log Fold Change", ylab="Log Odds")

#2.4. Protein-protein interaction (PPI) network construction
biocLite("STRINGdb")
library("STRINGdb")
# 9606 for human
string_db=STRINGdb$new(version="10",species=9606,score_threshold=0,input_directory="")
# mapping to the DEG
myDF=string_db$map(siggenes[(logsel),],"ID",removeUnmappedRows=TRUE)
# set color for up and down regulated
myDF_reg <- string_db$add_diff_exp_color(subset(myDF, logFcColStr="logFC"))
payload_id <- string_db$post_payload(myDF_reg$STRING_id, colors = myDF_reg$color)
# plot
string_db$plot_network(myDF$STRING_id,payload_id=payload_id)

#2.5. Network-based neighborhood scoring analysis
dim(myDF)

#2.7 Support vector machine classification
#extract the expression set
colnames(exprs(esetSel))
x.train  <- t(exprs(esetSel))
#extract the pheno data 
y.train <- as.factor(colnames(exprs))

# Create learning and test set
x.learn    <- x.train
y.learn    <- y.train

#get these from GSE17538 (Just GSE17536 because the other one ues a different chip. We can add it in seperatly later if we need to)
setwd("//fileserv/Profiles$/tweedya/Desktop/project 1/data/GSE17536_RAW")
fn<-dir()  
str(fn) 

valid_rawAffyData <- ReadAffy(filenames=fn, cdfname="hgu133plus2cdf")  

valid_eset <- rma(valid_rawAffyData)


# extract the data matrix from the expression set
#need to select only genes that are in the other set
#replace with gene symbols
#only grab the genes we're interested in 
valid_exprs <- exprs(valid_eset)
rownames(valid_exprs) <-unlist(mget(sid,hgu133plus2SYMBOL,ifnotfound=NA))  
valid_exprs_subset <- valid_exprs[(logsel ), ]
#check the dimensions
dim(valid_exprs_subset)

#need to label which samples are recur vs not recur 
gse_17536 <- getGEO(filename = "//fileserv/Profiles$/tweedya/Desktop/project 1/data/GSE17536_series_matrix.txt.gz")
valid_recur_data <- pData(gse_17536[1])$characteristics_ch1.7

no_recur = grep("(: no recurrence)",valid_recur_data)
recur = grep("(: recurrence)",valid_recur_data)
na = grep("(: NA)",valid_recur_data)

colnames(valid_exprs_subset)[no_recur] <- "NoRecur"
colnames(valid_exprs_subset)[recur] <- "Recur"
colnames(valid_exprs_subset)[na] <- "NA"
summary(as.factor(colnames(valid_exprs_subset)))
#col_vals <- c("NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","Recur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","Recur","NoRecur","Recur","NoRecur","NA","Recur","NoRecur","NoRecur","Recur","NoRecur","NoRecur","Recur","NoRecur","NoRecur","Recur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","Recur","NoRecur","NoRecur","Recur","Recur","NoRecur","Recur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NA","NoRecur","NoRecur","NoRecur","Recur","NoRecur","NoRecur","NoRecur","NoRecur","Recur","NoRecur","Recur","Recur","NoRecur","NoRecur","Recur","Recur","Recur","Recur","Recur","NoRecur","NoRecur","Recur","Recur","Recur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","Recur","NoRecur","Recur","Recur","Recur","NoRecur","NoRecur","NA","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","Recur","NoRecur","NoRecur","Recur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","NoRecur","Recur","Recur","NoRecur","NoRecur","NoRecur","Recur","NoRecur","NoRecur","NoRecur","NA","NoRecur","NA","NA","NA","NA","NA","NoRecur","NA","NA","NA","NA","Recur","NA","NA","NA","NA","Recur","NoRecur","NA","NA","NA","NA","NA","NoRecur","NA","NA","NA","NA","NA","Recur","NA","NA","NoRecur","NoRecur","Recur","NA","NA","NA")
#colnames(valid_exprs_subset) <-  col_vals
valid_exprs_subset

#Filter out the NA Values 
sel_cols <- colnames(valid_exprs_subset) != "NA"
x.valid    <- t(valid_exprs_subset[,sel_cols])
y.valid    <- as.factor(colnames(valid_exprs_subset[,sel_cols]))

dim(x.learn)
dim(x.valid)

library("e1071")
#help(svm)

svm.error <- matrix(0, nrow = 3, ncol = 3)  ## error will be saved
for(cost in 0:2)     ## cost = 1, 2, 4
{
  for(gamma in (-1):1)   ## gamma = 2^gamma/ncol(x.learn)
  {
    i              <- cost+1   ## index for error matrix
    j              <- gamma+2
    svm.fit        <- svm(x.learn, y.learn, cost = 2^cost,    
                          gamma = 2^gamma/ncol(x.learn))
    svm.predic     <- predict(svm.fit, newdata = x.valid)
    svm.error[i,j] <- mean(svm.predic != y.valid)
  }
}

#2.8 Validation 
## Done with svm - you may print the error matrix
svm.error

