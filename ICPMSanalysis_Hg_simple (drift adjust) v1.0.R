#Simple script for ICP-MS data analysis of mercury (single metal, with drift adjust) 
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
#...Sequence_position and Sequence_name column titles, underscores to metals (Hg 202 becomes Hg_202),
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

#Create empty dataframe (Intensities_Hg) and start loop command
Mercury_analysis <- data.frame()
for (n in seq_pos) {
  #Get row number(idx) corresponding to sequence position
  idx <- which(Intensities$Sequence_position == n)
  
  #Retrieve RSD values from Hg_202 column in Intensities (one row below 'idx')
  df_rsd <- Intensities[idx+1, "Hg_202"]
  #Use function to remove % sign
  df_rsd <- lapply(df_rsd, function(x) { 
    gsub("%", "", x)
  })
  #Convert to data frame
  df_rsd <- as.data.frame(df_rsd, check.names = FALSE)
  colnames(df_rsd) <- "RSD_percent"
  
  #Use column bind (cbind) to assemble newly named columns into a new data frame (df)
  #combine Seq_pos, Sample name, Hg_202 intensity, RSD values
  df <- cbind(data.frame(Sequence_position = n, 
                         Sample_name = Intensities[idx,"Sample_name"]), 
              "Hg_202 (cps)" = Intensities$Hg_202[idx], 
              "RSD_percent" = df_rsd,
              "Ir_193 (cps)" = Intensities$Ir_193[idx])
  
  #At each iteration, stack df to new data frame (Intensities_Hg)
  Mercury_analysis <- rbind(Mercury_analysis, df)
}
#Look at resulting Hg dataframe with RSD as new column
View(Mercury_analysis)


#STEP THREE: Account for machine sampling variation using Ir_193 internal standard

#Convert Hg and Ir values from characters to numeric using as.numeric function.
#Add new columns to Mercury_analysis for: 
#Irfactor (all Ir values arbitrarily divided by 200000) 
#Iradjust (Hg values divided by Irfactor to adjust for amount of sample analyzed)
Mercury_analysis$Z01_Irfactor <- as.numeric(Mercury_analysis$Ir_193) / 200000
Mercury_analysis$Z02_Iradjust <- as.numeric(Mercury_analysis$Hg_202) / Mercury_analysis$Z01_Irfactor


#STEP FOUR: Continue Hg analysis by adjusting samples for baseline drift using 1 ppb standards

#Create data frame with only 1 ppb standard samples (Mercury_standards)
Mercury_standards <- filter(Mercury_analysis, Sample_name == "1 ppb standard")
#Use trendline function to create scatter plot for 1 ppb standards with trendline and equation
trendline(Mercury_standards$Sequence_position, Mercury_standards$Z02_Iradjust,model = "line2P", 
          Pvalue.corrected = TRUE, linecolor = "blue", lty = 1, lwd = 3, 
          show.equation = TRUE, show.Rpvalue = TRUE, Rname = 1, Pname = 0, xname = "x", yname = "y",
          yhat = FALSE, summary = TRUE, text.col = "black", main="Baseline Drift (1 ppb standards)", 
          xlab = "Sequence position", ylab = "Signal (Hg)", CI.fill = FALSE, CI.lty = "blank")
#Delete from here to remove drift adjust. Remove Z03 drift adjust appearances later and replace with Z02
#lm function allows use of m (gradient, coef[2]) in equation to apply drift adjustment
#This gives new column (Z03_Hg_driftadjust) in Mercury_analysis
Mercury_driftplot <- lm(Mercury_standards$Z02_Iradjust ~ Mercury_standards$Sequence_position)
Mercury_analysis$Z03_Hg_driftadjust <- Mercury_analysis$Z02_Iradjust - coef(Mercury_driftplot)[2] * 
  (Mercury_analysis$Sequence_position - Mercury_standards[1,1])


#STEP FIVE: Continue Hg analysis by plotting blanks and deducting blank average from all values

#Make data frame with only Blank samples (Mercury_blanks) and plot to check visually for outliers
Mercury_blanks <- filter(Mercury_analysis, Sample_name == "Blank")
plot(Mercury_blanks$Sequence_position, Mercury_blanks$Z03_Hg_driftadjust, 
     main="Trend of blank values", xlab = "Sequence position", ylab = "Signal (Hg)")
#Take average of all Blank values and use this to blank correct all other values
#This gives new column (Z04_Hg_blankadjust) in Mercury_analysis
Mercury_analysis$Z04_Hg_blankadjust <- Mercury_analysis$Z03_Hg_driftadjust - mean(Mercury_blanks[["Z03_Hg_driftadjust"]])


#STEP SIX: Plot Hg standard curve and use gradient to convert all sample Hg signals from cps to ppb

#Group (under Mercury_standard curve) then rename all standard curve samples and convert to numeric
Mercury_standardcurve <- filter(Mercury_analysis, Sample_name == "0.0004 ppb metal standard" 
                             | Sample_name == "0.002 ppb metal standard" | Sample_name == "0.02 ppb metal standard"
                             | Sample_name == "0.1 ppb metal standard" | Sample_name == "0.5 ppb metal standard"
                             | Sample_name == "1 ppb metal standard" | Sample_name == "5 ppb metal standard" 
                             | Sample_name == "10 ppb metal standard")
Mercury_stdconccharacters <- gsub(" ppb metal standard", "", Mercury_standardcurve$Sample_name) 
Mercury_standardcurve$Concentration <- as.numeric(Mercury_stdconccharacters)
#Plot trendline and use gradient to convert signal to ppb for all samples
trendline(Mercury_standardcurve$Concentration, Mercury_standardcurve$Z04_Hg_blankadjust,model = "line2P", 
          Pvalue.corrected = TRUE, linecolor = "blue", lty = 1, lwd = 3, 
          show.equation = TRUE, show.Rpvalue = TRUE, Rname = 1, Pname = 0, xname = "x", yname = "y",
          yhat = FALSE, summary = TRUE, text.col = "black", main="Hg standard curve", 
          xlab = "Hg concentration (ppb)", ylab = "Signal", CI.fill = FALSE, CI.lty = "blank")
Mercury_stdplot <- lm(Mercury_standardcurve$Z04_Hg_blankadjust ~ Mercury_standardcurve$Concentration)
#This gives new column (Z05_Hg_ppb) in Mercury_analysis
Mercury_analysis$Z05_Hg_ppb <- Mercury_analysis$Z04_Hg_blankadjust / coef(Mercury_stdplot)[2]


#STEP SEVEN: Calculate Hg concentration in original samples by accounting for 1/50 dilution during sample prep/digestion

#Sample_Hg_ppb (column Z07) - Adjust sample Hg values by dilution factor from digestion and sample prep
Mercury_analysis$Z07_Sample_Hg_ppb <- Mercury_analysis$Z05_Hg_ppb * 50


#STEP EIGHT: Apply RSDs and summarise results for Hg (original sample and SD in ppb, detection limit)

#Multiply all Sample_Hg_ppb values by numeric equivalent of (RSD/100) to get error
Mercury_analysis$Z08Hg_RSD_ppb <- Mercury_analysis$Z07_Sample_Hg_ppb * (as.numeric(Mercury_analysis$RSD_percent) / 100)

#Use 3 * STD of Blanks to get detection limit (in ppb) to appear below results
Mercury_detectionlimit_ppb <- (3 * sd(Mercury_blanks[["Z03_Hg_driftadjust"]])) / coef(Mercury_stdplot)[2]

#Arrange results summary in new data frame (Mercury_results)
#Include Sample_name, Hg_ppb detected in measured samples, Hg_ppb concentration in original sample, RSD_percent, RSD_ppb
Mercury_results <- select(Mercury_analysis, Sample_name, Z07_Sample_Hg_ppb, RSD_percent, Z08Hg_RSD_ppb, Z05_Hg_ppb)

#Rename column titles to give output of results in console
Mercury_results <- data.frame(Sample_name = Mercury_results$Sample_name, 
                              Hg_ppb_original_sample = Mercury_results$Z07_Sample_Hg_ppb,
                              RSD_percent = Mercury_results$RSD_percent,
                              RSD_ppb = Mercury_results$Z08Hg_RSD_ppb, 
                              Hg_ppb_detected_in_measured_sample = Mercury_results$Z05_Hg_ppb, 
                              Detection_limit_ppb = Mercury_detectionlimit_ppb)
View(Mercury_results)


#STEP NINE: Save data to file

#Write results to file, first creating Results folder in working directory
if(!dir.exists("Results")) dir.create("Results")
#Make csv file of Mercury_results dataframe
write.csv(Mercury_results, file = file.path("Results", "Mercury_results.csv"), row.names = FALSE)


#######################################
#END OF DATA ANALYSIS
#######################################


#View step by step data analysis
View(Mercury_analysis)
#View raw data (csv files) to bring dataframes up in tabs and check them
View(Intensities)




#Clear plots
dev.off()
#Clear environment
rm(list = ls())
#Clear console
cat("\f")