#suppose you have a already normalized values from the expression datasets, how to obtain the model fitting
#and applying differential expression. I was asked to do this before an interview.
#R based analysis
library(Biobase)
library(affy)
library(limma)
library(DataEditR) # package for datascience
#Visualizing the data and cutting the tables in dataeditR and reading the matrix into the frame and 
#then converting into an expression set
read <- read.table("normalised_data_curated.csv", row.names = 1, header = T)
data_edit(read)
View(read)
dim(read)
data <- as.matrix(read)
eset <- new("ExpressionSet", exprs=data)
#######################
#Making a design matrix of the expression matrix and defining the model and the contrasts.
Estimate the fold changes and standard errors by fitting a linear model for each gene
design <- model.matrix(~ 0+factor(c(0,0,0,1,1,1,2,2,2)))
colnames(design) <- c("group1", "group2", "group3")
fit <- lmFit(eset, design)
contrast.matrix <- makeContrasts(group2-group1, group3-group2, group3-group1, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
# Apply empirical Bayes smoothing to the standard errors.
fit2B <- eBayes(fit2) 
fitdataframe <- as.data.frame(fit2B)
head(fitdataframe)
# to see the top 10 differentially expressed genes for coef 1
topTable(fit2B, coef=1, adjust="BH") for the group2-group1 
# to see the top 10 differentially expressed genes for coef 2
topTable(fit2B, coef=2, adjust="BH") for the group3-group2 
# showing the top 100 differentially expressed genes for the coef 1
topTable(fit2B, coef=1, number = 100, adjust="BH") 
# exporting all the differentially expressed from the coef 1
topall1 <- topTable(fit2B, coef=1, number=Inf) 
# exporting all the differentially expressed from the coef 2
topall2 <- topTable(fit2B, coef=2, number=Inf) 
# exporting all the differentially expressed from the coef 2
topall3 <- topTable(fit2B, coef=3, number=Inf) 
results <- decideTests(fit2B)
summary(results)
vennDiagram(results)
volcanoplot(fit2)
#writing the csv file for the topall1, coef 1
write.csv(topall1,file="~/Desktop/group2-group1_diff_exp_all1.csv") 
#writing the csv file for the topall1, coef 2
write.csv(topall2,file="~/Desktop/group2-group1_diff_exp_all2.csv") 
#writing the csv file for the topall1, coef 3
write.csv(topall3,file="~/Desktop/group2-group1_diff_exp_all3.csv") 
# filtering the 0.05 p-value genes for the coef 1
all_1_0.5 <- topall1[topall1$P.Value < 0.05,] 
# writing the 0.05 p-value genes for the coef 1
write.csv(topall_1_0.5,file="~/Desktop/group2_group1_exp_p_0.05.csv") 
# filtering the 0.01 p-value genes for the coef 1
all_1_0.1 <- topall1[topall1$P.Value < 0.01,] 
# writing the 0.05 p-value genes for the coef 1
write.csv(all_1_0.1,file="~/Desktop/group2_group1_exp_p_0.01.csv") 
# filtering the 0.05 p-value genes for the coef 2
all_2_0.5 <- topall2[topall2$P.Value < 0.05,] 
# writing the 0.05 p-value genes for the coef 2
write.csv(all_2_0.5,file="~/Desktop/group3_group2_exp_p_0.05.csv") 
# filtering the 0.01 p-value genes for the coef 2
all_2_0.1 <- topall2[topall2$P.Value < 0.01,] 
# writing the 0.05 p-value genes for the coef 2
write.csv(all_2_0.1,file="~/Desktop/group3_group2_exp_p_0.01.csv") 
# filtering the 0.05 p-value genes for the coef 3
all_3_0.5 <- topall3[topall3$P.Value < 0.05,] 
# writing the 0.05 p-value genes for the coef 3
write.csv(all_3_0.5,file="~/Desktop/group3_group1_exp_p_0.05.csv") 
# filtering the 0.01 p-value genes for the coef 3
all_3_0.1 <- topall3[topall3$P.Value < 0.01,] 
# writing the 0.05 p-value genes for the coef 3
write.csv(all_1_0.1,file="~/Desktop/group3_group1_exp_p_0.01.csv") 
### for enrichment, i took the 0.1 from all the coef #####
Yekutieli (FDR under dependency)
