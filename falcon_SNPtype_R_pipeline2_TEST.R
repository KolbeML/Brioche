library(reshape)
library(adegenet)

#This function edits data from a given input csv file, and saves it to a new csv file.
edit_data <- function(inputFile, outputFile) {
  data <- read.csv(inputFile, header =T, skip=15) #import and get rid of lines at top
  keeps <- c("Assay","Name","Converted") #keeps columns with locus name, ID of individual (woodrat or falcon, etc.), snp genotype
  data2 <- data[keeps]
  
  data2 <- data2[- grep("NTC",data2$Converted),] #gets rid of the negative control
  data2$Converted <- gsub(':','_',data2$Converted) 
  data2$Converted <- gsub('A','1',data2$Converted) #switches out nucleotides for different numbers
  data2$Converted <- gsub('C','2',data2$Converted)
  data2$Converted <- gsub('G','3',data2$Converted)
  data2$Converted <- gsub('T','4',data2$Converted)
  data2$Converted <- gsub('No 2all','NA_NA',data2$Converted) #downstream programs like NA_NA better
  write.csv(data3,outputFile,row.names=FALSE) #exports data
}

# This function merges two given CSV files by a given column 
# header - input the name of your first column (value in cell A1?)
# - and returns the result.
merge_csvs <- function(file1, file2, sortBy){
  set1 <- read.csv(file1,header=T) 
  set2 <- read.csv(file2,header=T)
  return(merge(set1,set2,by=sortBy))
}

edit_data("run1_set1.csv","plate1_set1_transposed.csv") #insert file names as appropriate
edit_data("run1_set2.csv","plate1_set2_transposed.csv")
edit_data("run2_set1.csv","plate2_set1_transposed.csv")
edit_data("run2_set2.csv","plate2_set2_transposed.csv")

# Get your "file" (dataframe) load down to 2 before proceeding to the fullData step
combined1 <- merge_csvs("plate1_set1_transposed.csv", "plate1_set2_transposed.csv", "Name");
combined2 <- merge_csvs("plate2_set1_transposed.csv", "plate2_set2_transposed.csv", "Name")

fullData <- rbind(combined2, combined1) # merges the two remaining dataframes - must be merged in this order
write.csv(fullData,"merged_plates.csv",row.names=FALSE)


#----Kolbe's Merge Script----#

mydf <- read.csv("merged_plates.csv", header=F) # replace w file name
badrows <- c("placeholder-unique") #I initialized this with a value that will hopefully never actually be in a file
badcols <- readLines("TESTbadcols.txt") # replace w appropriate file name. should contain column headers on individual lines

print(badcols)

cols <- ncol(mydf) #number columns in the file
rows <- nrow(mydf) #number rows in the file
acceptable_NANA <- 10 # replace w your upper bound for acceptable number of NA_NA's

numBadCols = 0
for(i in 1:cols - 1)
{
  if (mydf[1, cols - i] %in% badcols)
  {
    cat("column ", cols - i + 1, " was in the bad columns list. its header was ", as.character(mydf[1, cols - i + 1]), "\n")
    mydf[cols - i + 1] <- NULL
    numBadCols <- numBadCols + 1
  }
  else
    cat("column ", cols - i + 1, " is a good column, its header is ", as.character(mydf[1, cols - i + 1]),"\n")
}
cat("Number of columns removed: ", numBadCols)

cols <- ncol(mydf)

numBadRows <- 0
for (i in 1:rows) #go through entire file, counting # NA_NAs per row. 
{
  count <- 0
  for (j in 1:cols)
  {
    if (mydf[i,j] == ("NA_NA"))
      count <- count + 1
  }
  cat("number of NANA in row ", i , " is ", count, "\n")
  
  if (count > acceptable_NANA) #if too many NA_NA...
  {
    badrows <- append(badrows, as.character(mydf[i, 1])) #...add first cell in row to badrows
    cat("Row ", i , " had too many NA_NAs. It's header was ", as.character(mydf[i,1]), "\n")
    numBadRows <- numBadRows + 1
  }
}
cat("Total number of bad rows: ", numBadRows)

mydf <- mydf[!(as.character(mydf$V1) %in% badrows), ] #remove all rows starting with an item in badrows from the dataframe.
write.table(mydf, file="edited_merged_plates.csv", sep=",", col.names = F, row.names = F)

#----End of Kolbe's Merge Script----#


combined <- read.csv("combined_greaterthan10.csv",header=TRUE)

library(splitstackshape)
library(allelematch)

cSplit(combined,2:164, sep="_") -> output

falcon_allelematch <- amDataset(output,indexColumn ="Name",ignoreColumn= 1:1, missingCode="NA")

unique_falcon_allelematch <- amUnique(falcon_allelematch, alleleMismatch=7)

summary(unique_falcon_allelematch, html="unique_falcon_allelematch.html")
summary(unique_falcon_allelematch, csv="falcon_allelematch_allresults.csv")
summary(unique_falcon_allelematch, csv="falcon_allelematch_uniqueindiv.csv",uniqueOnly=TRUE)

uniqueindiv = read.csv("falcon_allelematch_uniqueindiv_edited.csv")
siteinfo = read.csv("falcon_field_data_for_R.csv")
total <- merge(siteinfo,uniqueindiv,by="DNA_ID")
write.csv(total,"submission3_allelematch_uniqueindiv_siteinfo.csv")

library(related)
input <- readgenotypedata("submission3_forRelated.txt")
outfile <- coancestry(input$gdata,lynchli=1,quellergt=1,ritland=1, wang=1,lynchrd=1)
write.csv(outfile$relatedness,"outfile.csv")

#To calculate LD for loci:
#Use PGDSpider to generate a ped file from a GDA file.
library(snpStats)
sample <- read.pedfile("sub4_forSNPstats.ped", snps="pedsnps.txt")
ldDprime_r2 = ld(x=sample$genotypes, depth=(ncol(sample$genotypes)-1), stats="R.squared", symmetric=TRUE)
d = as.matrix(ldDprime_r2)
write.csv(d,"d.csv")

Linkage (d > 0.20)

1612964_1419399_ASPN_selection X 1612964_1456785_ECM2_nonsyn -> d = 0.27 (remove 1612964_1419399_ASPN_selection)
1615281_1321851_AGA_selection X 1612891_582304_AIM1_selection -> d = 0.23 (remove 1612891_582304_AIM1_selection)

#Pairwise locus Fst values from diveRsity in R:
#Export genepop file from GenAlEx, population names can not have spaces.

diff_stats <- diffCalc(infile="genepop_for_diveRsity.txt",pairwise=TRUE, fst=TRUE)
diff_stats
write.csv(diff_stats$pw_locus,"diff_stats.csv")
write.csv(diff_stats$pairwise,"diff_stats_pairwise.csv")

#Alternatively for 95% CIs:
diff_stats <- diffCalc(infile="genepop_for_diveRsity.txt",pairwise=TRUE, fst=TRUE, bs_pairwise=TRUE, boots =999)


#Identifying outlier loci with BAYESCAN:

#Used PGDSpider to convert genepop file to Bayescan file to test for selection based on FST values.
#Used Bayescan to test for selection based on F statistics.

#Parameters:


#To look at output from Bayescan:
#Copy and paste plot_R.r script into R.
mydata2=read.table(file.choose(),colClasses="numeric")
plot_bayescan(mydata2,FDR=0.05)

Outlier locus = 102 = X1613580_160905_A2ML1_nonsyn_1 = SNP42 in Set2 for Fluidigm = GenAlEx also found significant departures from HW for this marker

#Reviewed results in Fluidigm SNP genotyping analysis software and they look pretty solid!  Clustering looks good.  Seems safe to move forward.
