#Simple script for ICP-MS data analysis of iron (single metal, with drift adjust) 
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

#Sample_prep csv file must be prepared to ensure that 'Sample_names' exactly match those in the ICPMS sequence
#Prepare Sample_run_intensities csv files by adding:
#...Sequence_position and Sequence_name column titles, underscores to metals (Fe 56 becomes Fe_56),
#...and removing columns for acquisition time, dataset file and method file
#...and deleting all worksheets except the one with both Intensities and RSDs
#Save this worksheet as a .csv file (Sample_run_XYZ)
#Files must be in same folder as working directory/project/script


#STEP ONE: Load relevant .csv files into RStudio

#Change file names and sample names (standards e.g. "1 ppb standard" or "10 ppb metal standard") in script below as required
#Create data frames for both .csv files (Sample_run_XYZ, Sample_digest_XYZ)
#Read as csv (read.csv) or excel file (read_excel from readxl package) 
Intensities <- read.csv("Sample_run_test.csv")
Sample_prep <- read.csv("Sample_digest_test.csv")


#STEP TWO: Rearrange 'Intensities' data frame to move RSD values from separate rows to new column
#Can always skip STEP TWO if this has already been done manually in Excel

#Retrieve sequences positions without 'NA's in new data frame (seq_pos)
seq_pos <- Intensities$Sequence_position
seq_pos <- seq_pos[!is.na(seq_pos)]

#Create empty dataframe (Intensities_Fe) and start loop command
Iron_analysis <- data.frame()
for (n in seq_pos) {
  #Get row number(idx) corresponding to sequence position
  idx <- which(Intensities$Sequence_position == n)
  
  #Retrieve RSD values from Fe_56 column in Intensities (one row below 'idx')
  df_rsd <- Intensities[idx+1, "Fe_56"]
  #Use function to remove % sign
  df_rsd <- lapply(df_rsd, function(x) { 
    gsub("%", "", x)
  })
  #Convert to data frame
  df_rsd <- as.data.frame(df_rsd, check.names = FALSE)
  colnames(df_rsd) <- "RSD_percent"
  
  #Use column bind (cbind) to assemble newly named columns into a new data frame (df)
  #combine Seq_pos, Sample name, Fe_56 intensity, RSD values
  df <- cbind(data.frame(Sequence_position = n, 
                         Sample_name = Intensities[idx,"Sample_name"]), 
              "Fe_56 (cps)" = Intensities$Fe_56[idx], 
              "RSD_percent" = df_rsd,
              "Sc_45 (cps)" = Intensities$Sc_45[idx])
  
  #At each iteration, stack df to new data frame (Intensities_Fe)
  Iron_analysis <- rbind(Iron_analysis, df)
}
#Look at resulting Fe dataframe with RSD as new column
View(Iron_analysis)


#STEP THREE: Account for machine sampling variation using Sc_45 internal standard

#Convert Fe and Sc values from characters to numeric using as.numeric function.
#Add new columns to Iron_analysis for: 
#Scfactor (all Sc values arbitrarily divided by 200000) 
#Scadjust (Fe values divided by Scfactor to adjust for amount of sample analyzed)
Iron_analysis$Z01_Scfactor <- as.numeric(Iron_analysis$Sc_45) / 200000
Iron_analysis$Z02_Scadjust <- as.numeric(Iron_analysis$Fe_56) / Iron_analysis$Z01_Scfactor


#STEP FOUR: Continue Fe analysis by adjusting samples for baseline drift using 1 ppb standards

#Create data frame with only 1 ppb standard samples (Iron_standards)
Iron_standards <- filter(Iron_analysis, Sample_name == "1 ppb standard")
#Use trendline function to create scatter plot for 1 ppb standards with trendline and equation
trendline(Iron_standards$Sequence_position, Iron_standards$Z02_Scadjust,model = "line2P", 
          Pvalue.corrected = TRUE, linecolor = "blue", lty = 1, lwd = 3, 
          show.equation = TRUE, show.Rpvalue = TRUE, Rname = 1, Pname = 0, xname = "x", yname = "y",
          yhat = FALSE, summary = TRUE, text.col = "black", main="Baseline Drift (1 ppb standards)", 
          xlab = "Sequence position", ylab = "Signal (Fe)", CI.fill = FALSE, CI.lty = "blank")
#Delete from here to remove drift adjust. Remove Z03 drift adjust appearances later and replace with Z02
#lm function allows use of m (gradient, coef[2]) in equation to apply drift adjustment
#This gives new column (Z03_Fe_driftadjust) in Iron_analysis
Iron_driftplot <- lm(Iron_standards$Z02_Scadjust ~ Iron_standards$Sequence_position)
Iron_analysis$Z03_Fe_driftadjust <- Iron_analysis$Z02_Scadjust - coef(Iron_driftplot)[2] * 
  (Iron_analysis$Sequence_position - Iron_standards[1,1])


#STEP FIVE: Continue Fe analysis by plotting blanks and deducting blank average from all values

#Make data frame with only Blank samples (Iron_blanks) and plot to check visually for outliers
Iron_blanks <- filter(Iron_analysis, Sample_name == "Blank")
plot(Iron_blanks$Sequence_position, Iron_blanks$Z03_Fe_driftadjust, 
     main="Trend of blank values", xlab = "Sequence position", ylab = "Signal (Fe)")
#Take average of all Blank values and use this to blank correct all other values
#This gives new column (Z04_Fe_blankadjust) in Iron_analysis
Iron_analysis$Z04_Fe_blankadjust <- Iron_analysis$Z03_Fe_driftadjust - mean(Iron_blanks[["Z03_Fe_driftadjust"]])


#STEP SIX: Plot Fe standard curve and use gradient to convert all sample Fe signals from cps to ppb

#Group (under Iron_standard curve) then rename all standard curve samples and convert to numeric
Iron_standardcurve <- filter(Iron_analysis, Sample_name == "0.0004 ppb metal standard" 
                             | Sample_name == "0.002 ppb metal standard" | Sample_name == "0.02 ppb metal standard"
                             | Sample_name == "0.1 ppb metal standard" | Sample_name == "0.5 ppb metal standard"
                             | Sample_name == "1 ppb metal standard" | Sample_name == "5 ppb metal standard" 
                             | Sample_name == "10 ppb metal standard")
Iron_stdconccharacters <- gsub(" ppb metal standard", "", Iron_standardcurve$Sample_name) 
Iron_standardcurve$Concentration <- as.numeric(Iron_stdconccharacters)
#Plot trendline and use gradient to convert signal to ppb for all samples
trendline(Iron_standardcurve$Concentration, Iron_standardcurve$Z04_Fe_blankadjust,model = "line2P", 
          Pvalue.corrected = TRUE, linecolor = "blue", lty = 1, lwd = 3, 
          show.equation = TRUE, show.Rpvalue = TRUE, Rname = 1, Pname = 0, xname = "x", yname = "y",
          yhat = FALSE, summary = TRUE, text.col = "black", main="Fe standard curve", 
          xlab = "Fe concentration (ppb)", ylab = "Signal", CI.fill = FALSE, CI.lty = "blank")
Iron_stdplot <- lm(Iron_standardcurve$Z04_Fe_blankadjust ~ Iron_standardcurve$Concentration)
#This gives new column (Z05_Fe_ppb) in Iron_analysis
Iron_analysis$Z05_Fe_ppb <- Iron_analysis$Z04_Fe_blankadjust / coef(Iron_stdplot)[2]


#STEP SEVEN: Calculate Fe concentration in original samples by accounting for dilution during sample prep/digestion

#Make Iron_analysis2 dataframe as subset of Iron_analysis dataframe, matching number of rows to Sample_prep.csv
Iron_analysis2 <- Iron_analysis %>% 
        filter(Sample_name != "0.0004 ppb metal standard" & Sample_name != "0.002 ppb metal standard"
           & Sample_name != "0.02 ppb metal standard" & Sample_name != "0.1 ppb metal standard" 
           & Sample_name != "0.5 ppb metal standard" & Sample_name != "1 ppb metal standard" 
           & Sample_name != "5 ppb metal standard" & Sample_name != "10 ppb metal standard" 
           & Sample_name != "Blank" & Sample_name != "1 ppb standard") %>%
        select(Sample_name, Z05_Fe_ppb, RSD_percent)

#Merge Sample_prep.csv data into Iron_analysis2, match rows by Sample_name
Iron_analysis2 <- merge(Iron_analysis2, 
                           Sample_prep[,c("Sample_name", "standard_factor_digestion", "Sample_dilution_factor")], 
                           by = "Sample_name")

#Wetashcorrect (column Z06) - Adjust metal concentration for standard dilution during Wet Ash digestion (Dry Ash = 1X)  
Iron_analysis2$Z06_Wetashcorrect <- Iron_analysis2$Z05_Fe_ppb / Iron_analysis2$standard_factor_digestion
#Sample_Fe_ppb (column Z07) - Adjust sample Fe values by dilution factor from digestion and sample prep
Iron_analysis2$Z07_Sample_Fe_ppb <- Iron_analysis2$Z06_Wetashcorrect / Iron_analysis2$Sample_dilution_factor


#STEP EIGHT: Apply RSDs and summarise results for Fe (original sample and SD in ppb, detection limit)

#Multiply all Sample_Fe_ppb values by numeric equivalent of (RSD/100) to get error
Iron_analysis2$Z08Fe_RSD_ppb <- Iron_analysis2$Z07_Sample_Fe_ppb * (as.numeric(Iron_analysis2$RSD_percent) / 100)

#Use 3 * STD of Blanks to get detection limit (in ppb) to appear below results
Iron_detectionlimit_ppb <- (3 * sd(Iron_blanks[["Z03_Fe_driftadjust"]])) / coef(Iron_stdplot)[2]

#Arrange results summary in new data frame (Iron_results)
#Include Sample_name, Fe_ppb detected in measured samples, Fe_ppb concentration in original sample, RSD_percent, RSD_ppb
Iron_results <- select(Iron_analysis2, Sample_name, Z07_Sample_Fe_ppb, RSD_percent, Z08Fe_RSD_ppb, Z05_Fe_ppb)

#Rename column titles to give output of results in console
Iron_results <- data.frame(Sample_name = Iron_results$Sample_name, 
                              Fe_ppb_original_sample = Iron_results$Z07_Sample_Fe_ppb,
                              RSD_percent = Iron_results$RSD_percent,
                              RSD_ppb = Iron_results$Z08Fe_RSD_ppb, 
                              Fe_ppb_detected_in_measured_sample = Iron_results$Z05_Fe_ppb, 
                              Detection_limit_ppb = Iron_detectionlimit_ppb)
View(Iron_results)


#STEP NINE: Save data to file

#Write results to file, first creating Results folder in working directory
if(!dir.exists("Results")) dir.create("Results")
#Make csv file of Iron_results dataframe
write.csv(Iron_results, file = file.path("Results", "Iron_results.csv"), row.names = FALSE)


#######################################
#END OF DATA ANALYSIS
#######################################


#View step by step data analysis
View(Iron_analysis)
View(Iron_analysis2)
#View raw data (csv files) to bring dataframes up in tabs and check them
View(Intensities)
View(Sample_prep)




#Clear plots
dev.off()
#Clear environment
rm(list = ls())
#Clear console
cat("\f")