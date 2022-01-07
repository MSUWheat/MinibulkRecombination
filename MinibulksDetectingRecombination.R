########################################################################################################
#Determining Number of Recombination Events Per Chromosome Arm in Minibulk Populations
########################################################################################################
library(ABHgenotypeR)
library(dplyr)
########################################################################################################
#MSU18000092 Example
########################################################################################################
rm(list = ls(all = TRUE))
setwd('/Users/Olson/Documents/Research/Minibulks/GenomeAnalysis/MSU18000092')
########################################################################################################
#Adjust wheat chromsome names to numeric
########################################################################################################
x<-read.table('MSU18000092.csv', sep=",", header=FALSE) #No header name so that marker names can be changed along with chromosome names
x<-as.data.frame(sapply(x,gsub,pattern="1A",replacement="11"))
x<-as.data.frame(sapply(x,gsub,pattern="2A",replacement="21"))
x<-as.data.frame(sapply(x,gsub,pattern="3A",replacement="31"))
x<-as.data.frame(sapply(x,gsub,pattern="4A",replacement="41"))
x<-as.data.frame(sapply(x,gsub,pattern="5A",replacement="51"))
x<-as.data.frame(sapply(x,gsub,pattern="6A",replacement="61"))
x<-as.data.frame(sapply(x,gsub,pattern="7A",replacement="71"))
x<-as.data.frame(sapply(x,gsub,pattern="1B",replacement="12"))
x<-as.data.frame(sapply(x,gsub,pattern="2B",replacement="22"))
x<-as.data.frame(sapply(x,gsub,pattern="3B",replacement="32"))
x<-as.data.frame(sapply(x,gsub,pattern="4B",replacement="42"))
x<-as.data.frame(sapply(x,gsub,pattern="5B",replacement="52"))
x<-as.data.frame(sapply(x,gsub,pattern="6B",replacement="62"))
x<-as.data.frame(sapply(x,gsub,pattern="7B",replacement="72"))
x<-as.data.frame(sapply(x,gsub,pattern="1D",replacement="13"))
x<-as.data.frame(sapply(x,gsub,pattern="2D",replacement="23"))
x<-as.data.frame(sapply(x,gsub,pattern="3D",replacement="33"))
x<-as.data.frame(sapply(x,gsub,pattern="4D",replacement="43"))
x<-as.data.frame(sapply(x,gsub,pattern="5D",replacement="53"))
x<-as.data.frame(sapply(x,gsub,pattern="6D",replacement="63"))
x<-as.data.frame(sapply(x,gsub,pattern="7D",replacement="73"))
x<-as.matrix(x); colnames(x)<-x[1,];x<-x[-1,];x<-as.data.frame(x) #Add back the headers and make it a dataframe

########################################################################################################
#Count the number of markers per chromosome and drop chromosomes with 2 or less markers
########################################################################################################
cm<-as.data.frame(t(x))
ct<-as.data.frame(table(cm$V1))
ct<-tabyl(cm$V1)

#Specify the S chromosome number to drop, if needed
#x <- x %>% select(-contains("S43"))

#write.table(x,"MSU18000092.1.csv",sep=",", row.names=FALSE)
#write.table(ct,"markersPerArm.csv",sep=",", row.names=FALSE)

########################################################################################################
#correctStretches with ABHgenotypeR
########################################################################################################
genotypes <- readABHgenotypes("MSU18000092.1.csv",  nameA = "Female", nameB = "Male", readPos = TRUE)
ErrCorr2Genotypes <- correctStretches(genotypes, maxHapLength = 2)
plotGenos(ErrCorr2Genotypes)
writeABHgenotypes(ErrCorr2Genotypes, outfile = "MSU18000092EXPORT.csv")
chr<-as.matrix(ErrCorr2Genotypes$chrom)
colnames(chr)<-"Chr"

########################################################################################################
#Create a ABH genotype dataframe from the genotypes list >>> geno
########################################################################################################
geno<-cbind(chr, t(ErrCorr2Genotypes$ABHmatrix))
geno[geno=="A"]<-1
geno[geno=="B"]<-0
geno<-as.data.frame(geno) 

########################################################################################################
#Create a lag row for each chromosome in geno >>> out
########################################################################################################
out<-geno %>% group_by(Chr) %>% mutate_all(lag, n = 1) #Adds a blank row before each chromosome

########################################################################################################
#Drop the Chr column from geno and out
########################################################################################################
geno<-geno[,-1]
out<-out[-1]

########################################################################################################
#Convert geno and out from factor to numeric
########################################################################################################
geno<- geno %>%
  mutate_all(funs(as.numeric(as.character(.)))) 
out<- out %>%
  mutate_all(funs(as.numeric(as.character(.))))

########################################################################################################
#Subtract 'geno' (original genotype dataframe) from 'out' (lag dataframe) to identify recombination 
########################################################################################################
dif<-geno-out #subracts the dataframes elementwise
dif2<-dif^2#square everything to make values positive, Non-zero values are recombination events
dif2<-cbind(ErrCorr2Genotypes$marker_names, chr, dif2) #add back the chromosome column
dif3<-dif2[complete.cases(dif2), ]#drop rows with NAs

########################################################################################################
#Return mean recombination per chromosome
########################################################################################################
chr<-dif2[,2]
chr<-as.data.frame(chr)
cc<-chr %>% count(chr)
chrs<-length(cc$chr)#get the number of chromosomes

#divide column sums by the number of chromosomes to get the mean number of recombination events per line
difSums<-as.data.frame(colSums(dif3[,-1:-2])/chrs)
overall<-colMeans(difSums)#Average recombination events per chromosome for population MSU18000092
roverall<-as.data.frame(overall)

#write.table(roverall, "MSU18000092.overall.recombination.csv", sep=",", row.names=FALSE)

########################################################################################################
#Return recombinanation events by chromosome
########################################################################################################
cdifs<-aggregate(. ~ Chr, dif3[,-1], sum)#calculate per chromosome recombination events in each line
cdifs<-cdifs[,-1]#remove the chromosome colum for the row means calculation
cmeans<-as.data.frame(rowMeans(cdifs)) #average number of recombination events per chromosome across all lines
cmeans<-cbind(cc$chr, cmeans)
colnames(cmeans)<-c("Chr", "Rec")
#write.table(cmeans, "MSU18000092.chromosome.recombination.csv", sep=",", row.names=FALSE)

save.image("MinibulksDetectingRecombination.RData")
