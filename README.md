### Data Inspection of `fang_et_al_genotypes.txt`

```sh
genotype <- as.data.frame(read.table("fang_et_al_genotypes.txt", sep="\t",header=TRUE))
```
1. View the no of rows and cloumn
```sh
dim(genotype)
```
2. Shows first 10 elements of `fang_et_al_genotypes.txt`
```sh
sapply(genotype, class)[1:10]
```
3. Shows first 10 column names of `fang_et_al_genotypes.txt`
```sh
colnames(genotype)[1:10]
```
4. Shows structure of `fang_et_al_genotypes.txt`
```sh
str(genotype)
```
5. Shows the no of groups in `fang_et_al_genotypes.txt`
```sh
genotype %>% group_by(Group) %>% count()
```
From inspecting `fang_et_al_genotypes.txt` file, I found
1. There are 2782 rows and 986 columns
2. Observed first 10 elements of the file
3. Observed first 10 columns of the file
4. Observed the structure of the file
5. There are 16 groups 
### Data Inspection of `snp_position.txt`
```sh 
snp <- as.data.frame(read.table("snp_position.txt", sep="\t",header=TRUE))
```
1. View the no of rows and cloumn
```sh
dim(snp)
```
2. Shows first 10 elements of `snp_position.txt`
```sh
sapply(snp, class)[1:10]
```
3. Shows first 10 column names of `snp_position.txt`
```sh
colnames(snp)[1:10]
```
4. Shows structure of `snp_position.txt`
```sh
str(snp)
```
5. overview of the data's structure in `snp_position.txt`
```sh
glimpse(snp)
```
From inspecting `snp_position.txt` file, I found
1. There are 983 rows and 15 columns
2. Observed first 10 elements of the file
3. Observed first 10 columns of the file
4. Observed the structure of the file
5. Overviewed the data structure of the file

### Data Processing
1. Loading the following packages and files to start working (if required)
```sh
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)
genotype <- as.data.frame(read.table("fang_et_al_genotypes.txt", sep="\t",header=TRUE))
snp <- as.data.frame(read.table("snp_position.txt", sep="\t",header=TRUE))
```
2. Converting all 'unknown' and 'multiple' to NA in ```snp_position.txt```
```sh
snp[snp == "unknown"] <- NA
snp[snp == "multiple"] <- NA
```
3. Selecting only SNP_ID, Chromosome, Position columns from dataframe and converting Chromosome and Position cloumn to numeric. Filters out any rows with 'NA'.
```sh
snpnew <- snp %>% dplyr::select(SNP_ID, Chromosome, Position) %>% mutate(Chromosome=as.numeric(Chromosome), Position=as.numeric(Position))%>% filter_all(all_vars(. != "NA"))
```
4. Subsetting maize and teosinte genotypes by the Group column
```sh
maize <- genotype[which(genotype$Group=="ZMMIL" | genotype$Group =="ZMMLR" | genotype$Group == "ZMMMR"),] 
```
```sh
teosinte <-genotype[which(genotype$Group=="ZMPBA" | genotype$Group =="ZMPIL" | genotype$Group == "ZMPJA"),]
```
5. selecting all columns of the maize data frame except the second and third columns.
```sh
maize <- maize[,c(-2,-3)]
```
6. Transposing maize data frame  , , 
```sh
maize <- t(maize)
```
7. Adding new cloumn with row names
```sh
maize <- cbind(rownames(maize),maize) 
```
8. Removing the row names from maize data frame
```sh
rownames(maize) <- NULL 
```
9. Setting the column names of the maize data frame to the first row values
```sh
colnames(maize) <- maize[1,] 
```
10. Removing the first row of the maize data frame
```sh
maize <- maize[-1,] 
```
11. Converting to data frame
```sh
maize <- as.data.frame(maize) 
```
12. Renaming first coulmn to 'SNP_ID'
```sh
colnames(maize)[1] <- "SNP_ID"
```
13. Merging the two data frames (snpnew and maize) based on SNP_ID and arranging them by Chromosome first then by Position.
```sh
maizemerge <- merge(snpnew, maize, by = "SNP_ID") 
maizemerge <- maizemerge %>% arrange(Chromosome,Position)
```
14. Follwing the same instructions for teosinte
```sh
teosinte <- teosinte[,c(-2,-3)] 
teosinte <- t(teosinte) 
teosinte <- cbind(rownames(teosinte),teosinte) 
rownames(teosinte) <- NULL 
colnames(teosinte) <- teosinte[1,] 
teosinte <- teosinte[-1,] 
teosinte <- as.data.frame(teosinte) 
colnames(teosinte)[1] <- "SNP_ID" 
teosintemerge <- merge(snpnew, teosinte, by = "SNP_ID") 
teosintemerge <- teosintemerge %>% arrange(Chromosome,Position)
```
15. Creating files for all 10 Chromosomes with increasing and decreasing chromosome position , replacing missing values with ?/- for maize
```sh
chr <- 1:10 
for (i in chr) { 
files_inc <- maizemerge[maizemerge$Chromosome == i,] 
files_inc[files_inc == "?/?"] <- "?"
if (i < 10) { write.table(files_inc, file = paste("Maize_Chr0",i,"_increase.txt",sep=""),row.names = FALSE,sep = 
"\t",quote = FALSE) } 
else {write.table(files_inc, file = 
paste("Maize_Chr",i,"_increase.txt",sep=""),row.names = FALSE, sep = 
"\t",quote = FALSE)} 
files_dec <- maizemerge[maizemerge$Chromosome == i,] 
files_dec[files_dec == "?/?"] <- "-" 
files_dec <- files_dec %>% arrange(desc(Chromosome),desc(Position)) 
if (i < 10) { write.table(files_dec, file = paste("Maize_Chr0",i,"_decrease.txt",sep=""),row.names = FALSE,sep = "\t",quote = FALSE) } 
else {write.table(files_dec, file = paste("Maize_Chr",i,"_decrease.txt",sep=""),row.names = FALSE, sep = 
"\t",quote = FALSE)} }
```
16. Creating files for all 10 Chromosomes with increasing and decreasing chromosome position , replacing missing values with ?/- for teosinte
```sh
chr <- 1:10 
for (i in chr) { 
files_inc <- teosintemerge[teosintemerge$Chromosome == i,] 
files_inc[files_inc == "?/?"] <- "?" 
if (i < 10) { write.table(files_inc, file= paste("Teosinte_Chr0",i,"_increase.txt",sep=""),row.names = FALSE,sep = 
"\t",quote = FALSE) } 
else {write.table(files_inc, file = paste("Teosinte_Chr",i,"_increase.txt",sep=""),row.names = FALSE, sep = "\t",quote = FALSE)}
files_dec <- teosintemerge[teosintemerge$Chromosome == i,] 
files_dec[files_dec == "?/?"] <- "-" 
files_dec <- files_dec %>% arrange(desc(Chromosome),desc(Position)) 
if (i < 10) { write.table(files_dec, file = paste("Teosinte_Chr0",i,"_decrease.txt",sep=""),row.names = FALSE,sep = "\t",quote = FALSE) } 
else {write.table(files_dec, file = paste("Teosinte_Chr",i,"_decrease.txt",sep=""),row.names = FALSE, sep = 
"\t",quote = FALSE)} 
}
```
### Plotting
```sh
snp %>% 
dplyr::select(SNP_ID, Chromosome, Position) %>%
drop_na() %>%
mutate(Chromosome = as.double(Chromosome)) %>%
ggplot() +
geom_bar(aes(x = Chromosome, fill = factor(Chromosome)), alpha = 0.7) + scale_x_continuous(breaks = 1:10) +
scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#C14B00", "#00BFC4"))
```
#### Distribution of SNP markers

Installing Packages
```sh
install.packages("viridis")
```
```sh
library(viridis)
```

```sh
sample_size = snpnew %>% group_by(Chromosome) %>% summarize(num=n())
snpnew %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(Chromosome, "\n", "n=", num)) %>%
  ggplot(aes(x = myaxis, y = Position, fill = as.character(Chromosome))) +
  geom_violin(width = 1.4) +
  geom_boxplot(width = 0.1, color = "grey", alpha = 0.2) +
  scale_fill_viridis(discrete = TRUE) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 11)
  ) +
  ggtitle("The distribution of SNP position on each chromosome") +
  xlab("Chromosome")
```
