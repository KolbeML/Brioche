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

cols <- ncol(mydf)

cSplit(combined, 2:cols, sep="_") -> output

