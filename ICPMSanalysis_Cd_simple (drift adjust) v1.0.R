#Simple script for ICP-MS data analysis of cadmium (single metal, with drift adjust) 
#Latest version of RStudio recommended


#R SET UP

#install then load tidyverse (contains all packages needed such as dplyer and ggplot2) and basicTrendline
#Make sure these are loaded each time this is used
if(!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)
if(!require("basicTrendline")) install.packages("basicTrendline")
library(basicTrendline)

#Set working directory - change user name and file path as required
setwd("C:/Users/tomas/Desktop/R working directory ICPMS")


#PREPARING EXCEL DOCUMENTS

#Prepare Sample_run_intensities csv files by adding:
#...Sequence_position and Sequence_name column titles, underscores to metals (Cd 111 becomes Cd_111),
#...and removing columns for acquisition time, dataset file and method file
#...and deleting all worksheets except the one with both Intensities and RSDs
#Save this worksheet as a .csv file (Sample_run_XYZ) in same folder as working directory/project/script


#STEP ONE: Load relevant .csv file into RStudio

#Change file names and sample names (standards e.g. "1 ppb standard" or "10 ppb metal standard") in script below as required
#Create data frames for .csv file (Sample_run_XYZ)
#Read as csv (read.csv) or excel file (read_excel from readxl package) 
Intensities <- read.csv("Sample_run_HeavyMetalstest_modified.csv")


#STEP TWO: Rearrange 'Intensities' data frame to move RSD values from separate rows to new column
#Can always skip STEP TWO if this has already been done manually in Excel

#Retrieve sequences positions without 'NA's in new data frame (seq_pos)
seq_pos <- Intensities$Sequence_position
seq_pos <- seq_pos[!is.na(seq_pos)]

#Create empty dataframe (Intensities_Cd) and start loop command
Cadmium_analysis <- data.frame()
for (n in seq_pos) {
  #Get row number(idx) corresponding to sequence position
  idx <- which(Intensities$Sequence_position == n)
  
  #Retrieve RSD values from Cd_111 column in Intensities (one row below 'idx')
  df_rsd <- Intensities[idx+1, "Cd_111"]
  #Use function to remove % sign
  df_rsd <- lapply(df_rsd, function(x) { 
    gsub("%", "", x)
  })
  #Convert to data frame
  df_rsd <- as.data.frame(df_rsd, check.names = FALSE)
  colnames(df_rsd) <- "RSD_percent"
  
  #Use column bind (cbind) to assemble newly named columns into a new data frame (df)
  #combine Seq_pos, Sample name, Cd_111 intensity, RSD values
  df <- cbind(data.frame(Sequence_position = n, 
                         Sample_name = Intensities[idx,"Sample_name"]), 
              "Cd_111 (cps)" = Intensities$Cd_111[idx], 
              "RSD_percent" = df_rsd,
              "In_115 (cps)" = Intensities$In_115[idx])
  
  #At each iteration, stack df to new data frame (Intensities_Cd)
  Cadmium_analysis <- rbind(Cadmium_analysis, df)
}
#Look at resulting Cd dataframe with RSD as new column
View(Cadmium_analysis)


#STEP THREE: Account for machine sampling variation using In_115 internal standard

#Convert Cd and In values from characters to numeric using as.numeric function.
#Add new columns to Cadmium_analysis for: 
#Infactor (all In values arbitrarily divided by 200000) 
#Inadjust (Cd values divided by Infactor to adjust for amount of sample analyzed)
Cadmium_analysis$Z01_Infactor <- as.numeric(Cadmium_analysis$In_115) / 200000
Cadmium_analysis$Z02_Inadjust <- as.numeric(Cadmium_analysis$Cd_111) / Cadmium_analysis$Z01_Infactor


#STEP FOUR: Continue Cd analysis by adjusting samples for baseline drift using 1 ppb standards

#Create data frame with only 1 ppb standard samples (Cadmium_standards)
Cadmium_standards <- filter(Cadmium_analysis, Sample_name == "1 ppb standard")
#Use trendline function to create scatter plot for 1 ppb standards with trendline and equation
trendline(Cadmium_standards$Sequence_position, Cadmium_standards$Z02_Inadjust,model = "line2P", 
          Pvalue.corrected = TRUE, linecolor = "blue", lty = 1, lwd = 3, 
          show.equation = TRUE, show.Rpvalue = TRUE, Rname = 1, Pname = 0, xname = "x", yname = "y",
          yhat = FALSE, summary = TRUE, text.col = "black", main="Baseline Drift (1 ppb standards)", 
          xlab = "Sequence position", ylab = "Signal (Cd)", CI.fill = FALSE, CI.lty = "blank")
#Delete from here to remove drift adjust. Remove Z03 drift adjust appearances later and replace with Z02
#lm function allows use of m (gradient, coef[2]) in equation to apply drift adjustment
#This gives new column (Z03_Cd_driftadjust) in Cadmium_analysis
Cadmium_driftplot <- lm(Cadmium_standards$Z02_Inadjust ~ Cadmium_standards$Sequence_position)
Cadmium_analysis$Z03_Cd_driftadjust <- Cadmium_analysis$Z02_Inadjust - coef(Cadmium_driftplot)[2] * 
  (Cadmium_analysis$Sequence_position - Cadmium_standards[1,1])


#STEP FIVE: Continue Cd analysis by plotting blanks and deducting blank average from all values

#Make data frame with only Blank samples (Cadmium_blanks) and plot to check visually for outliers
Cadmium_blanks <- filter(Cadmium_analysis, Sample_name == "Blank")
plot(Cadmium_blanks$Sequence_position, Cadmium_blanks$Z03_Cd_driftadjust, 
     main="Trend of blank values", xlab = "Sequence position", ylab = "Signal (Cd)")
#Take average of all Blank values and use this to blank correct all other values
#This gives new column (Z04_Cd_blankadjust) in Cadmium_analysis
Cadmium_analysis$Z04_Cd_blankadjust <- Cadmium_analysis$Z03_Cd_driftadjust - mean(Cadmium_blanks[["Z03_Cd_driftadjust"]])


#STEP SIX: Plot Cd standard curve and use gradient to convert all sample Cd signals from cps to ppb

#Group (under Cadmium_standard curve) then rename all standard curve samples and convert to numeric
Cadmium_standardcurve <- filter(Cadmium_analysis, Sample_name == "0.0004 ppb metal standard" 
                             | Sample_name == "0.002 ppb metal standard" | Sample_name == "0.02 ppb metal standard"
                             | Sample_name == "0.1 ppb metal standard" | Sample_name == "0.5 ppb metal standard"
                             | Sample_name == "1 ppb metal standard" | Sample_name == "5 ppb metal standard" 
                             | Sample_name == "10 ppb metal standard")
Cadmium_stdconccharacters <- gsub(" ppb metal standard", "", Cadmium_standardcurve$Sample_name) 
Cadmium_standardcurve$Concentration <- as.numeric(Cadmium_stdconccharacters)
#Plot trendline and use gradient to convert signal to ppb for all samples
trendline(Cadmium_standardcurve$Concentration, Cadmium_standardcurve$Z04_Cd_blankadjust,model = "line2P", 
          Pvalue.corrected = TRUE, linecolor = "blue", lty = 1, lwd = 3, 
          show.equation = TRUE, show.Rpvalue = TRUE, Rname = 1, Pname = 0, xname = "x", yname = "y",
          yhat = FALSE, summary = TRUE, text.col = "black", main="Cd standard curve", 
          xlab = "Cd concentration (ppb)", ylab = "Signal", CI.fill = FALSE, CI.lty = "blank")
Cadmium_stdplot <- lm(Cadmium_standardcurve$Z04_Cd_blankadjust ~ Cadmium_standardcurve$Concentration)
#This gives new column (Z05_Cd_ppb) in Cadmium_analysis
Cadmium_analysis$Z05_Cd_ppb <- Cadmium_analysis$Z04_Cd_blankadjust / coef(Cadmium_stdplot)[2]


#STEP SEVEN: Calculate Cd concentration in original samples by accounting for 1/50 dilution during sample prep/digestion

#Sample_Cd_ppb (column Z07) - Adjust sample Cd values by dilution factor from digestion and sample prep
Cadmium_analysis$Z07_Sample_Cd_ppb <- Cadmium_analysis$Z05_Cd_ppb * 50


#STEP EIGHT: Apply RSDs and summarise results for Cd (original sample and SD in ppb, detection limit)

#Multiply all Sample_Cd_ppb values by numeric equivalent of (RSD/100) to get error
Cadmium_analysis$Z08Cd_RSD_ppb <- Cadmium_analysis$Z07_Sample_Cd_ppb * (as.numeric(Cadmium_analysis$RSD_percent) / 100)

#Use 3 * STD of Blanks to get detection limit (in ppb) to appear below results
Cadmium_detectionlimit_ppb <- (3 * sd(Cadmium_blanks[["Z03_Cd_driftadjust"]])) / coef(Cadmium_stdplot)[2]

#Arrange results summary in new data frame (Cadmium_results)
#Include Sample_name, Cd_ppb detected in measured samples, Cd_ppb concentration in original sample, RSD_percent, RSD_ppb
Cadmium_results <- select(Cadmium_analysis, Sample_name, Z07_Sample_Cd_ppb, RSD_percent, Z08Cd_RSD_ppb, Z05_Cd_ppb)

#Rename column titles to give output of results in console
Cadmium_results <- data.frame(Sample_name = Cadmium_results$Sample_name, 
                              Cd_ppb_original_sample = Cadmium_results$Z07_Sample_Cd_ppb,
                              RSD_percent = Cadmium_results$RSD_percent,
                              RSD_ppb = Cadmium_results$Z08Cd_RSD_ppb, 
                              Cd_ppb_detected_in_measured_sample = Cadmium_results$Z05_Cd_ppb, 
                              Detection_limit_ppb = Cadmium_detectionlimit_ppb)
View(Cadmium_results)


#STEP NINE: Save data to file

#Write results to file, first creating Results folder in working directory
if(!dir.exists("Results")) dir.create("Results")
#Make csv file of Cadmium_results dataframe
write.csv(Cadmium_results, file = file.path("Results", "Cadmium_results.csv"), row.names = FALSE)


#######################################
#END OF DATA ANALYSIS
#######################################


#View step by step data analysis
View(Cadmium_analysis)
#View raw data (csv files) to bring dataframes up in tabs and check them
View(Intensities)




#Clear plots
dev.off()
#Clear environment
rm(list = ls())
#Clear console
cat("\f")