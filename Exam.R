#R script from EXAM

# Load packages
library(cliqueMS)
library(xcms)
library(enviPat)
library(RColorBrewer)
library(ggplot2)
library(magrittr)
library(plotly)
library(DT)


# XCMS Analysis 
## Load data
### List .mzXML files contained in the working directory

path.mzXML <- "C:/Users/Bernat/Desktop/UNI/3r 2010-2021/Omicas/Practicas/PractFINAL" 

mzXML.files <- list.files(path.mzXML, pattern= ".mzXML", 
                         recursive=T, full.names = T)



#2.1 How many replicates were measured for each cell line? ¿Cuántas réplicas se midieron para cada línea celular?

### Phenodata dataframe from the working directory structure

pheno <- phenoDataFromPaths(mzXML.files)
pheno$sample_name <- rownames(pheno) 
pheno$sample_group <- pheno$class


#The phenoDataFromPaths function builds a data.frame representing the experimental design from the folder structure in which the files of the experiment are located.
#Each row is a serum sample.

### Reading data .mzxML files with readMSData method from the MSnbase package.
raw_data <- readMSData(files = mzXML.files,
                       pdata = new("NAnnotatedDataFrame", pheno),
                       mode = "onDisk")

#mode = "onDisk": reads only spectrum header from files, but no data which is retrieved on demand allowing to handle very large experiments with high memory demand.
#It does not upload it to RAM, only when I tell it to (on demand). It reads only the metadata, not the mass specs, which is what it weighs. 

##2.1
#### **Number of experimental groups and number of samples per group: ** 
table(pheno$sample_group)
#As it is shown, there are 6 experimental groups: CTR, DIA_0, DIA_18, PIO_0, PIO_18 and QC, with 14, 5, 6, 6, 6, and 8 samples respectively. 

##2.2 & 2.3
#### **How much time did the entire chromatographic run span?**
#In order to see the entire chromatographic run time we plot TIC: 

## We define colors according to experimental groups
group_colors <- paste0(brewer.pal(length(levels(pData(raw_data)$sample_group)),
                                  "Dark2"), "60")
names(group_colors) <- levels(pData(raw_data)$sample_group)

## We plot TIC using chromatogram() function 
## raw_data is the s4 object
TICs <- chromatogram(raw_data, aggregationFun = "sum") 

# Plot of TIC according to experimental groups
plot(TICs, col = group_colors[raw_data$sample_group], main = "TIC")
#We note the entire chromatographic run takes 600 seconds, that is to say, 10 minutes. 

#.	2.4 Extract Adenine ion chromatogram ([M+H]+) from raw data. 
  #Which sample does apparently hold higher levels of Adenine?. 
  #Plot the Adenine extracted ion chromatogram for this sample (0.25) 

cw_onlinexcms_default <- CentWaveParam(ppm = 15, 
                                       peakwidth = c(5, 10),snthresh = 6, 
                                       prefilter = c(3, 100),mzCenterFun = "wMean", 
                                       integrate = 1L,mzdiff = -0.001, fitgauss = FALSE, 
                                       noise = 0, verboseColumns = FALSE, 
                                       roiList = list(), firstBaselineCheck = TRUE, 
                                       roiScales = numeric())

## rt and m/z range of the Adenine peak area 
M <- 135.054495185 #monoisotopic weight of Adenine 
H <- 1.007276 #proton weight
MH <- M + H #theoretical weight = 136.0618

## Adenine peak mz range width according to XCMS parameters
error <- cw_onlinexcms_default@ppm #errors in ppm 
mmu_min <- MH-(error*MH/1e6) #proton weight of Adenine +- a ppm error: 15 ppm 
mmu_max <- MH+(error*MH/1e6)
mzr <- c(mmu_min, mmu_max) #m/z range

## Adenine peak rt range width according to XCMS parameters
apex.Adenine <- 263 #retention time ~ 263 s
rt.window <- cw_onlinexcms_default@peakwidth[2] #max rt width in seconds
rtr <- c(apex.Adenine -rt.window, apex.Adenine +rt.window)

#### **Extracted Ion Chromatograph (EIC) of Adenine from raw data** 
chr_Adenine  <- chromatogram(raw_data, mz = mzr, rt = rtr)
plot(chr_Adenine, col = group_colors[chr_Adenine$sample_group], lwd = 2)

#This is the Extracted Ion Chromatogram of Adenine.
#It plots the chromatogram in a mass range and in a retention time range.
#We note that we have 45 peaks, one for each sample, and colored according to their experimental group. 
#We are in a mass range that goes from a minimum mass range (mmu_min) of 150,0561 to a maximum mass range (mmu_max) of 150,0606. 
#The retention time is 293 s, because when we know that a metabolite is present in the samples, we know which is its the retention time.
#We observe that the maximum intensity at the apex (maxo) is 80000 approximately. 


#2.4 Which sample does apparently hold higher levels of Adenine?. 
#HCC70  green
#MCF7   orange --> EL PIC MÉS ALT DEL EIC 
#MDA231  purple
#MDA436 pink 
#MDA468 green fluix 
#QC yellow

#Mirem a quina mostra fa referencia a cada color.
group_colors
#veiem en hexadecimal que el taronja és la més alta, que correspon a la mosytra MCF7
#per tant, MCF7 parently hold higher levels of Adenine

#2.4 Plot the Adenine extracted ion chromatogram for this sample (0.25)
rtr2 <- filterRt(rt=c(apex.adenine, apex.adenine+1))
chr_most_adenine_<-chromatogram(raw_data,mz = mzr,rt = rtr2)
plot(chr_most_adenine)

spMS1.Lmethionine<- raw_data %>%
  filterRt(rt = c(apex.Lmethionine, apex.Lmethionine+1)) %>%
  filterFile(1) %>%
  spectra
plot(spMS1.Lmethionine[[1]])

#Extract and plot the XIC for Serine
data %>%
  filterRt(rt = c(175, 189)) %>%
  filterMz(mz = c(106.02, 106.07)) %>%
  chromatogram(aggregationFun = "max") %>%
  plot() 
## Plot raw data
raw_data %>%
  filterRt(rt = rtr) %>%
  filterMz(mz = mzr) %>%
  plot(type = "XIC")



#Exercise 3: 
#Run a full XCMS workflow. Explain sequentially each one of the steps in your own words. 
#Use centwave peaks detection, obiwarp retention time alignment and peak density for correspondence.

## Step1: Chromatographic peak detection through centwave
CentWaveParam()
#it tell us the default parameters (already defined)


#We must modify these parameters according what it is described in the on-line xcms parameters
#We define xcms online parameters 

cw_onlinexcms_default <- CentWaveParam(ppm = 15, 
                                       peakwidth = c(5, 20),snthresh = 6, 
                                       prefilter = c(3, 100),mzCenterFun = "wMean", 
                                       integrate = 1L,mzdiff = -0.001, fitgauss = FALSE, 
                                       noise = 0, verboseColumns = FALSE, 
                                       roiList = list(), firstBaselineCheck = TRUE, 
                                       roiScales = numeric())

xdata <- findChromPeaks(raw_data, param = cw_onlinexcms_default) 

#* findChromPeaks(): to detect chromatographic peaks. We will have an s4 object that stores the peaks it has found for each sample.
#* The findChromPeaks() function will give us a result that will be a list of all the peaks it finds with the centWave function in the 45 samples we have in raw_data. 
#* xData is the s4 object that contains all the xcms processing results. 
#* When the findChromPeaks() finishes running, the results will be stored in the variable xData

#3.1 How many peaks were detected on average per sample? 

#### ** Number of peaks detected on average per sample:**   
xdata  
#* Chromatographic peak detection:
#  + Method we have used: centWave
#+ On average 20678 chromatographic peaks per sample.




## Step 2: Alignment 
#There can be a little bit of movement of the RT from sample to sample.
#These minimums alignments can be corrected and we do use adjustRtime() wrapper function. 
#We have different algorithms to adjust RT. We are using obiwarp. 
#ObiwarpParam() warps the (full) data to a reference sample.

#Settings for the alignment:
  #We are going to use the online xcms parameters: profStep= 1 --> binSize = 1 
xdata <- adjustRtime(xdata, param = ObiwarpParam(binSize = 1))
#Sample number 16 used as center sample
#It uses sample number 16 as a sample against which you align the remaining ones. 
#It is done for all the samples. 

#3.2 Was retention time deemed necessary?. Plot difference between raw and adjusted retention times (1).
#We inspect difference between raw and adjusted retention times.
plotAdjustedRtime(xdata, col = group_colors[xdata$sample_group]) 
#EXPLICACIÓ 
#We  represent which are the RT applied for each sample.
  #Each line represents a sample. 
  #Between 100 and 300 s we have adjusted more or less like 7s above and 2 s below. 
#cada color és la diferència entre la mostra amb el rawdata i la mateixa mostra amb el rtime adjusted 
#Alignment results
#Shown is the difference between raw and adjusted retention times and the hook peaks that were used for the alignment (shown as points).
#The difference between raw and adjusted retention time should be reasonable. 
#In our example it is mostly below one second, which is OK since the samples were measured within a short time period 
#and differences are thus expected to be small.
#-------------------

## We use adjustedRtime parameter to access raw/adjusted retention times
par(mfrow = c(1, 2), mar = c(4, 4.5, 0.9, 0.5))
plot(chromatogram(xdata, mz = mzr,rt = rtr, adjustedRtime = FALSE), 
     col = group_colors[xdata$sample_group], peakType = "none", lwd = 2,
     main = "Before alignment")
plot(chromatogram(xdata, mz = mzr, rt = rtr),
     col = group_colors[xdata$sample_group], peakType = "none", lwd = 2,
     main = "After alignment")

#Abans del alignment ja estava bastant ben alineat Després del alignment, veiem que està millor alineat. Per tant, podem concloure observant el plot
#que no és estrictament necessari alinearho pero que si ho alinem queda millor i que al alinearho it does not performs worse after the alignment. 
#so qe can work OK with the alignment. 
#The base peak chromatograms are nicely aligned after retention time adjustment.



## Step 3: Correspondence 
#Correspondence group peaks or signals from the same ion across samples. 

### Defining peak density parameters

#Setting on-line xcms default parameters for peak density using PeakDensityParam().

#We define the online xcms peak density parameters: onlinexcms_pdp (go to online xcms and we use the parameters we see in the alignment section).

onlinexcms_pdp <- PeakDensityParam(sampleGroups = pheno$class,
                                   minFraction = 0.8, bw = 5, 
                                   binSize = 0.015)


#Performance of the correspondence analysis using optimized peak density settings on on-line xcms.
#In order to see how correspondence works we use groupChromPeaks().
xdata <- groupChromPeaks(xdata, param = onlinexcms_pdp) 
## xdata contains my features. 

### Correspondence Results 
#Extracting the feature definitions as data frame.
#We can see the correspondence results with featureDefinitions().
feature.info <- as.data.frame(featureDefinitions(xdata))

#3.3 How many features were finally detected? 
#### **Number of detected features:**
dim(feature.info)
#We have found 27398  features.
xdata  
#* Correspondence:
#  + 27398  features have been identified 




## Step 4: Missing Values 

#The feature matrix we will create later (fmat) contains NA values for samples in which no chromatographic peak was detected in the feature's m/z-rt region. We must fill in missing peaks in order to fill in the empty areas, and reduce the number of NA values in our matrix. 
## We get missing values before filling in peaks
apply(featureValues(xdata, filled = FALSE), MARGIN = 2,
      FUN = function(z) sum(is.na(z)))


#We use the fillChromPeaks method to fill in intensity data for such missing values from the original files. 
xdata <- fillChromPeaks(xdata)


#EXERCICSE 4 
#4.1 Filter detected features using an intensity minimum threshold = 3000. How many features out of the initially detected are retained? 

### Get intensity values for each feature--------------------------------------
D <- as.data.frame(featureValues(xdata, value = "maxo", method = "maxint"))
colnames(D) <- pData(xdata)$sample_name

## Compute the mean intensities for each group-----------------------
class <- pData(xdata)$class
median.intensities <- as.data.frame(t(apply(D,
                                            1, function(x) tapply(x, class, median, na.rm = TRUE))))
median.intensities$QC <- NULL #remove QC group
## Stablish intensity threshold value--------------------------------
thresholdvalue <- 3000
#Getting number of features with mean intensity above certain 
#threshold counts in at least one of the groups except QC group------
idx_i <- names(which(apply(median.intensities,
                           1, function(x) any(x>thresholdvalue)==TRUE)==TRUE))
cat(paste(length(idx_i), "out of", dim(median.intensities)[1]), 
    "above the threshold intensity value")
#Sabemos que en nuestro caso el threshold tiene que ser de 10.000 
#Miramos que features están en una intensidad media por encima de 10.000, porque si no ya no me las voy a quedar. 
#Nos dice que sólo 5679 de las 26591 están en esta intensidad media por encima. 
## 5679 out of 26591 above the threshold intensity value
#11883 out of 27398 above the threshold intensity value

#4.2 How many features out of the initially detected meet the 80% rule? Explain this rule in your own words. What is the rationale behind it?
### Compute number of samples per experimental group and minimum percent ---
n.samplespergroup <- table(pData(xdata)$class) 
percent <- 0.8 #apply the desired percentage
samples.min <- round(n.samplespergroup*percent,0)

samples_per_group <- feature.info[,c(grep("npeaks", 
                                          colnames(feature.info))+1:length(samples.min))]
samples_per_group_QC <- samples_per_group[colnames(samples_per_group) != "QC"]
## remove QCs
samples.min_QC <- samples.min[names(samples.min)!= "QC"]
resta <- apply(samples_per_group_QC, 1,  function(x) 
{any(x-as.vector(samples.min_QC)>0)}
)
idx_s <- rownames(D)[which(resta==TRUE)]

cat(length(idx_s), "out of", dim(D)[1], 
    "found in 80% of the samples of at least one of the experimental groups")
#20027 out of 27398 found in 80% of the samples of at least one of the experimental groups

#MARTA: 80% rule ens diu que ens quedarem amb aquelles features que estiguin almenys en un 80% de les mostres d'algun dels 5 grups experimentals
samples.min

#Em quedaré amb la feature si està representada en almenys 4 de HCC70, o en 5 dels MCF7, o en 5 dels MDA231, o en 4 de MDA436, o en 4 de MDA468. 
#Per exemple, si una feature l'he trobat en 4 de HCC70 i en cap altre lloc, me la quedaré perquè ja compliex que al menys està en un 80% de les mostres d'un dels grups experimentals.
#Em desfaig d'asquelles features que no compleixin aquesta regla. 

#Hi ha alguma raó per triar aquesta regla del 80%

#4.3 How many features out of the initially detected hold a peak width below 40 seconds? 
#que features tienen una anchura del pico menor que 40 segundos. 
#se puede sacar de la matriz de features 
#si me refiero a features estoy refiriendo a l amatriz de featurs 

#calculem rtmax - rtmin de totes les files i afegim la resta al final del dataframe 
feature.info$peakwidth40 = feature.info[,6] - feature.info[,5] 
#calculem el number of features out of the initially detected that hold a peak width below 40 seconds
num_features_peakwidth_below40s <- length(feature.info$peakwidth40[(feature.info$peakwidth40)<40])
#There are 27376 features out of the initially detected that hold a peak width below 40 seconds

#4.4 How many features out of the initially detected enclose higher biological than analytical variation?
## small function to calculate RSD%---------------------------------
RSD <- function(x){100*sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE)}

## Define  QC idx--------------------------------------------------
QChits <- which(class == "QC")

## Compute RSDs across samples and QCs----------------------------
RSD_samples <- apply(D[,-QChits],1, RSD)
RSD_QC <- apply(D[, QChits],1, RSD)

## Filter those features with RSD(QCs) > RSD (Samples)------------
idx_qc <- names(RSD_samples)[which(RSD_QC< RSD_samples)]

cat(paste(length(idx_qc), "out of", dim(median.intensities)[1]), 
    "hold higher biological than analytical variation")
#24121 features out of 27398 hold higher biological than analytical variation

#4.5 Select features that meet all filtering criteria. How many features resulted from the intersection of all filters?
# Merge intensity  sample representation and QCs criteria---------
data2stats <- D[intersect(intersect(idx_i, idx_s),idx_qc),-c(grep("QC", class))]
colnames(data2stats) <- pData(xdata)$sample_name[which(pData(xdata)$class!="QC")]
dim(data2stats)
  
#De las 27398 originales me he quedado con 9867. El producto de aplicar todo el proceso de el filtering, me ha disminuido mucho 
#el número de features que yo llevo para mi estadística. Ha pasado de 26591 a 3953, es decir, solo 3953 features nos interesan 
#para llevarlas a hacer experimentos MS/MS, o cumplen unos criterios estadísticos que hacen que podamos obtener MS/MS con sentido de nuestros datos. 

#4.5 Perform a PCA exploratory analysis considering them. Explain main trends in the data in terms of PC1 vs PC2 scores plot. 
##  Row-wise normalzation------------------------------------------------
D.norm <- as.data.frame(apply(data2stats, 1, function(x) (x/max(x, na.rm = TRUE))))
D.norm[is.na(D.norm)] <- 0
##  Compute PCA----------------------------------------------------------
pca_data <- prcomp(D.norm,retx=TRUE, center=TRUE, scale=FALSE)
summary(pca_data)$importance[,c(1:3)]


## We declare new class BRCA1 and non-BRCA1------------------------
cl <- as.character(pData(xdata)$class)

#BRCA1: HCC70, MDA-MB-436 and MDA-MB-468
cl[grep("4|C7", cl)] <- "BRCA1"
#nonBRCA1: MCF-7 and MDA-MB-231
cl[grep("MC|231", cl)] <- "nonBRCA1"





## Scores data
scores <- data.frame(pca_data$x[, c("PC1", "PC2")])
scores$class <- cl[-c(grep("QC", class))]
scores$lab <- rownames(D.norm)

scores.plot <- ggplot(data = scores, aes(x = PC1, y = PC2, colour = class, label=lab)) + 
  geom_point(alpha = I(0.7), size = 4) + 
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlab(paste("PC1 (", round(summary(pca_data)$importance[2,1], 2) * 100, "%)"))+
  ylab(paste("PC2 (", round(summary(pca_data)$importance[2,2], 2) * 100, "%)"))+
  stat_ellipse() + 
  theme_bw()
ggplotly(scores.plot)
#EXPLICAR PCA 

















#4.6 Impute NA values in the matrix resulting from filtered features using the minimum value in this matrix. 
#Compare between non-BRCA1-like (MCF-7 and MDA-MB-231) and BRAC-1 like (HCC70, MDA-MB-436 and MDA-MB-468) cell lines 
#computing Fold Changes and the Student's t-test. Correct for multiple testing using FDR (false discovery rate). 
#Consider FC>5 and p-corrected values<0.001 for statistical significance. How many features wese considered to be significantly 
#varying between the two phenotypes?. Plot the corresponding volcano plot and comment on it.

#We get the minimmum value in data2stast matrtix (We see the minimum value in the matrix is 200)
#We replace NA values in data2stats with the minimum value in the matrix.
data2stats[is.na(data2stats)]<- min(data2stats, na.rm=T)



#Compare between non-BRCA1-like (MCF-7 and MDA-MB-231) and BRAC-1 like (HCC70, MDA-MB-436 and MDA-MB-468) cell lines 
#computing Fold Changes and the Student's t-test. Correct for multiple testing using FDR (false discovery rate). 
#Consider FC>5 and p-corrected values<0.001 for statistical significance. How many features wese considered to be significantly 
#varying between the two phenotypes?. Plot the corresponding volcano plot and comment on it.

#remove the QC of the cl vector
cl_no <- cl[cl!="QC"] 
## Perform an ANOVA comparison for the non-BRCA1-like (MCF-7 and MDA-MB-231) and BRAC-1 like (HCC70, MDA-MB-436 and MDA-MB-468) cell lines ---------------------------------------------
gr <- as.factor(cl_no)
pm <- matrix(ncol=1,nrow=dim(data2stats)[1])
#convert pm to a matrix 
#as.matrix(pm)

for(i in 1:nrow(data2stats)){ 
  aov.out <- aov(as.numeric(data2stats[i,]) ~ gr)
  multcomp <- TukeyHSD(aov.out)
  pm[i,]  <- as.matrix(multcomp$gr[,"p adj"])
}
rownames(pm) <- rownames(data2stats)
colnames(pm) <- rownames(multcomp$gr)

#Adjust for multiple testing using false discovery rate
p.val.adj.nonBRCA1.BRCA1 <- p.adjust(pm[,"nonBRCA1-BRCA1"],"fdr")

## Create a function to compute FC-----------------------------------------------
fc.test <- function(df, classvec, class.case, class.control) {
  coln <- names(df)
  case <- df[which(coln == class.case)]
  control <- df[which(coln == class.control)]
  logFC <- log2(case/control)
  FC <- case/control
  FC2 <- -control/case
  FC[FC < 1] <- FC2[FC < 1]
  fc.res <- c(FC, logFC)
  names(fc.res) <- c("FC", "logFC")
  return(fc.res)
}

# Calculate median intensities
median.intensities <- t(apply(data2stats, 1, tapply, cl_no, median))

# Calculate FC for 'nonBRCA1-BRCA1' groups
fc.nonBRCA1_BRCA1 <- as.data.frame(t(apply(median.intensities, 1, function(x)
  fc.test(x, classvec = cl_no, class.case = "nonBRCA1", class.control = "BRCA1"))))
colnames(fc.nonBRCA1_BRCA1) <- c("FC_nonBRCA1_BRCA1", "logFC_nonBRCA1_BRCA1")


idx_nonBRCA1_BRCA1 <- which(abs(fc.nonBRCA1_BRCA1$FC)> 5 & 
                              p.val.adj.nonBRCA1.BRCA1 < 0.001)

cat(paste(length(idx_nonBRCA1_BRCA1), "out of", dim(data2stats)[1]), 
    " features significantly varying in nonBRCA1 vs BRCA1")
#2029 out of 9867  features significantly varying in nonBRCA1 vs BRCA1



# Draw Volcano plot for either comparisons
RES <- cbind.data.frame(median.intensities, 
                        fc.nonBRCA1_BRCA1, 
                        p.val.adj.nonBRCA1.BRCA1)
RES$labels <- rownames(RES)

p1 <- ggplot(data = RES, 
             aes(x = RES$logFC_nonBRCA1_BRCA1, 
                 y = -log10(RES$p.val.adj.nonBRCA1.BRCA1),
                 label = labels)) +
  geom_point(alpha = 0.4, size = 1.75) + theme(legend.position = "none") +
  xlim(-10, 10) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = log2(2), linetype = "dashed") + 
  geom_vline(xintercept = -log2(2), linetype = "dashed") + 
  labs(x = "log2(FC)", y = "-log10(p.adj)", 
       title = "nonBRCA1 vs BRCA1") +
  theme_bw()

ggplotly(p1)






#4.7 Download the HMDB and adducts tables and search those significant features using exact mass. 
    #Summarise results to MS/MS. Save the final results summary to a ".csv" file


## Read adduct information
positive <- read.table("C:/Users/Bernat/Desktop/UNI/3r 2010-2021/Omicas/Practicas/PractFINAL/ID/Positive_ESI.txt", header = T, stringsAsFactors = F)
# Read HMDB database
hmdb <- read.table("C:/Users/Bernat/Desktop/UNI/3r 2010-2021/Omicas/Practicas/PractFINAL/ID/HMDB_17_09.txt", header = T, sep = "\t", stringsAsFactors = F, 
                   fill = T, quote = "")





# list of significant features' m/z
input.mz <- feature.info[rownames(RES), "mzmed"]
input.rt <- feature.info[rownames(RES), "rtmed"]
# ppm error to perform the search. Instrument dependant (e.g. orbitrap 2-6;
# qTOF 8-12)
ppm <- 10

hmdbhits <- lapply(input.mz, function(x) {
  search_vect <- (x + positive$AdductMass)  # here we are using positive ionization
  h <- lapply(search_vect, function(x) {
    mass_range <- c(x - ppm * (x/1e+06), x + ppm * (x/1e+06))
    a <- which(hmdb$MonoisotopicMass > mass_range[1] & hmdb$MonoisotopicMass < 
                 mass_range[2])
  })
  
  h2 <- sapply(1:length(h), function(x) {
    y <- h[[x]]
    if (length(y) > 0) {
      r <- paste(hmdb$Name[y], sep = "", collapse = ";")
    } else {
      r <- "No hit"
    }
    return(r)
  })
  
  return(h2)
})
hmdbhits <- do.call("rbind", hmdbhits)
colnames(hmdbhits) <- positive$Adduct

# Summarize results  to MSMS
RESULTS2MSMS <- cbind.data.frame(RES,input.mz, input.rt, hmdbhits)

#We save the final results summary to a ".csv" file
write.csv(RESULTS2MSMS, file="results2MSMS.csv")

#4.8 Plot the levels of the metabolic feature with a rtmed~221 seconds and mzmed = 166.0713. 
#Which are its associated putative identifications considering the [M+H]+ adduct? 
#Is there more than a plausible identification? If yes, which ones?

#Subset the dataframe RESULTS2MSMS to get the feature with rtmed~221 seconds and mzmed = 166.0713
#-	Mi input.mz será las mzmed de todas las features que me han salido significantes
#-	Input.rt es el retention time (rtmed)

My_feature_subset <-RESULTS2MSMS[which.min(abs(221-RESULTS2MSMS$input.rt)),]
My_feature_subset
View(My_feature_subset)
#hem obtingut la feature FT01844

#Which are its associated putative identifications considering the [M+H]+ adduct? 
#Is there more than a plausible identification? If yes, which ones?

  #per veure que hi ha al N+H
    My_feature_subset[,"M+H", drop=FALSE]
  #FT01844 7-Methylguanine;3-Methylguanine;1-Methylguanine;N2-Methylguanine
  
  #per saber quin dels anteriors és mirem el boxplot de la figura 3 de l'articiule i veiem que
  #en 1-Methylguanine els nonbRCA1 tenen la 1-Methylguanine més elevada que els que tenen mutacions. 
  
  
  
  # Plot boxplot for the feature putatively identified as Methionine Sulfoxide
  #Boxplot amb la feature que hem trobat 
  boxplot(as.numeric(D["FT01844",])~pData(xdata)$class,
          xlab="class",
          ylab="FT01844")
  #per saber quin dels anteriors és mirem el boxplot de la figura 3 de l'articiule i veiem que
  #en 1-Methylguanine els nonbRCA1 tenen la 1-Methylguanine més elevada que els que tenen mutacions.
  #Podem observar que en el nostre boxplot passa el mateix, per tant serà la 1-Methylguanine
  #### In order to make a boxplot with ggplot2, we must convert D (matrix) to data frame
  
  
  
  
  
  
  
  
  
  
  
  
  
  
#---------SHORTCUTS-------#
save.image("Exam.Rdata")
setwd("C:/Users/Bernat/Desktop/UNI/3r 2010-2021/Omicas/Practicas/PractFINAL")
load("Exam.Rdata") 

