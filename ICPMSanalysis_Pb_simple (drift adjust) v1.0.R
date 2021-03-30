#Simple script for ICP-MS data analysis of lead (single metal, with drift adjust) 
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
#...Sequence_position and Sequence_name column titles, underscores to metals (Pb 208 becomes Pb_208),
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

#Create empty dataframe (Intensities_Pb) and start loop command
Lead_analysis <- data.frame()
for (n in seq_pos) {
  #Get row number(idx) corresponding to sequence position
  idx <- which(Intensities$Sequence_position == n)
  
  #Retrieve RSD values from Pb_208 column in Intensities (one row below 'idx')
  df_rsd <- Intensities[idx+1, "Pb_208"]
  #Use function to remove % sign
  df_rsd <- lapply(df_rsd, function(x) { 
    gsub("%", "", x)
  })
  #Convert to data frame
  df_rsd <- as.data.frame(df_rsd, check.names = FALSE)
  colnames(df_rsd) <- "RSD_percent"
  
  #Use column bind (cbind) to assemble newly named columns into a new data frame (df)
  #combine Seq_pos, Sample name, Pb_208 intensity, RSD values
  df <- cbind(data.frame(Sequence_position = n, 
                         Sample_name = Intensities[idx,"Sample_name"]), 
              "Pb_208 (cps)" = Intensities$Pb_208[idx], 
              "RSD_percent" = df_rsd,
              "Ir_193 (cps)" = Intensities$Ir_193[idx])
  
  #At each iteration, stack df to new data frame (Intensities_Pb)
  Lead_analysis <- rbind(Lead_analysis, df)
}
#Look at resulting Pb dataframe with RSD as new column
View(Lead_analysis)


#STEP THREE: Account for machine sampling variation using Ir_193 internal standard

#Convert Pb and Ir values from characters to numeric using as.numeric function.
#Add new columns to Lead_analysis for: 
#Irfactor (all Ir values arbitrarily divided by 200000) 
#Iradjust (Pb values divided by Irfactor to adjust for amount of sample analyzed)
Lead_analysis$Z01_Irfactor <- as.numeric(Lead_analysis$Ir_193) / 200000
Lead_analysis$Z02_Iradjust <- as.numeric(Lead_analysis$Pb_208) / Lead_analysis$Z01_Irfactor


#STEP FOUR: Continue Pb analysis by adjusting samples for baseline drift using 1 ppb standards

#Create data frame with only 1 ppb standard samples (Lead_standards)
Lead_standards <- filter(Lead_analysis, Sample_name == "1 ppb standard")
#Use trendline function to create scatter plot for 1 ppb standards with trendline and equation
trendline(Lead_standards$Sequence_position, Lead_standards$Z02_Iradjust,model = "line2P", 
          Pvalue.corrected = TRUE, linecolor = "blue", lty = 1, lwd = 3, 
          show.equation = TRUE, show.Rpvalue = TRUE, Rname = 1, Pname = 0, xname = "x", yname = "y",
          yhat = FALSE, summary = TRUE, text.col = "black", main="Baseline Drift (1 ppb standards)", 
          xlab = "Sequence position", ylab = "Signal (Pb)", CI.fill = FALSE, CI.lty = "blank")
#Delete from here to remove drift adjust. Remove Z03 drift adjust appearances later and replace with Z02
#lm function allows use of m (gradient, coef[2]) in equation to apply drift adjustment
#This gives new column (Z03_Pb_driftadjust) in Lead_analysis
Lead_driftplot <- lm(Lead_standards$Z02_Iradjust ~ Lead_standards$Sequence_position)
Lead_analysis$Z03_Pb_driftadjust <- Lead_analysis$Z02_Iradjust - coef(Lead_driftplot)[2] * 
  (Lead_analysis$Sequence_position - Lead_standards[1,1])


#STEP FIVE: Continue Pb analysis by plotting blanks and deducting blank average from all values

#Make data frame with only Blank samples (Lead_blanks) and plot to check visually for outliers
Lead_blanks <- filter(Lead_analysis, Sample_name == "Blank")
plot(Lead_blanks$Sequence_position, Lead_blanks$Z03_Pb_driftadjust, 
     main="Trend of blank values", xlab = "Sequence position", ylab = "Signal (Pb)")
#Take average of all Blank values and use this to blank correct all other values
#This gives new column (Z04_Pb_blankadjust) in Lead_analysis
Lead_analysis$Z04_Pb_blankadjust <- Lead_analysis$Z03_Pb_driftadjust - mean(Lead_blanks[["Z03_Pb_driftadjust"]])


#STEP SIX: Plot Pb standard curve and use gradient to convert all sample Pb signals from cps to ppb

#Group (under Lead_standard curve) then rename all standard curve samples and convert to numeric
Lead_standardcurve <- filter(Lead_analysis, Sample_name == "0.0004 ppb metal standard" 
                             | Sample_name == "0.002 ppb metal standard" | Sample_name == "0.02 ppb metal standard"
                             | Sample_name == "0.1 ppb metal standard" | Sample_name == "0.5 ppb metal standard"
                             | Sample_name == "1 ppb metal standard" | Sample_name == "5 ppb metal standard" 
                             | Sample_name == "10 ppb metal standard")
Lead_stdconccharacters <- gsub(" ppb metal standard", "", Lead_standardcurve$Sample_name) 
Lead_standardcurve$Concentration <- as.numeric(Lead_stdconccharacters)
#Plot trendline and use gradient to convert signal to ppb for all samples
trendline(Lead_standardcurve$Concentration, Lead_standardcurve$Z04_Pb_blankadjust,model = "line2P", 
          Pvalue.corrected = TRUE, linecolor = "blue", lty = 1, lwd = 3, 
          show.equation = TRUE, show.Rpvalue = TRUE, Rname = 1, Pname = 0, xname = "x", yname = "y",
          yhat = FALSE, summary = TRUE, text.col = "black", main="Pb standard curve", 
          xlab = "Pb concentration (ppb)", ylab = "Signal", CI.fill = FALSE, CI.lty = "blank")
Lead_stdplot <- lm(Lead_standardcurve$Z04_Pb_blankadjust ~ Lead_standardcurve$Concentration)
#This gives new column (Z05_Pb_ppb) in Lead_analysis
Lead_analysis$Z05_Pb_ppb <- Lead_analysis$Z04_Pb_blankadjust / coef(Lead_stdplot)[2]


#STEP SEVEN: Calculate Pb concentration in original samples by accounting for 1/50 dilution during sample prep/digestion

#Sample_Pb_ppb (column Z07) - Adjust sample Pb values by dilution factor from digestion and sample prep
Lead_analysis$Z07_Sample_Pb_ppb <- Lead_analysis$Z05_Pb_ppb * 50


#STEP EIGHT: Apply RSDs and summarise results for Pb (original sample and SD in ppb, detection limit)

#Multiply all Sample_Pb_ppb values by numeric equivalent of (RSD/100) to get error
Lead_analysis$Z08Pb_RSD_ppb <- Lead_analysis$Z07_Sample_Pb_ppb * (as.numeric(Lead_analysis$RSD_percent) / 100)

#Use 3 * STD of Blanks to get detection limit (in ppb) to appear below results
Lead_detectionlimit_ppb <- (3 * sd(Lead_blanks[["Z03_Pb_driftadjust"]])) / coef(Lead_stdplot)[2]

#Arrange results summary in new data frame (Lead_results)
#Include Sample_name, Pb_ppb detected in measured samples, Pb_ppb concentration in original sample, RSD_percent, RSD_ppb
Lead_results <- select(Lead_analysis, Sample_name, Z07_Sample_Pb_ppb, RSD_percent, Z08Pb_RSD_ppb, Z05_Pb_ppb)

#Rename column titles to give output of results in console
Lead_results <- data.frame(Sample_name = Lead_results$Sample_name, 
                              Pb_ppb_original_sample = Lead_results$Z07_Sample_Pb_ppb,
                              RSD_percent = Lead_results$RSD_percent,
                              RSD_ppb = Lead_results$Z08Pb_RSD_ppb, 
                              Pb_ppb_detected_in_measured_sample = Lead_results$Z05_Pb_ppb, 
                              Detection_limit_ppb = Lead_detectionlimit_ppb)
View(Lead_results)


#STEP NINE: Save data to file

#Write results to file, first creating Results folder in working directory
if(!dir.exists("Results")) dir.create("Results")
#Make csv file of Lead_results dataframe
write.csv(Lead_results, file = file.path("Results", "Lead_results.csv"), row.names = FALSE)


#######################################
#END OF DATA ANALYSIS
#######################################


#View step by step data analysis
View(Lead_analysis)
#View raw data (csv files) to bring dataframes up in tabs and check them
View(Intensities)




#Clear plots
dev.off()
#Clear environment
rm(list = ls())
#Clear console
cat("\f")