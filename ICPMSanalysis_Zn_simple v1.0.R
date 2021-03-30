#Simple script for ICP-MS data analysis of zinc (single metal, no drift adjust) 
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
#...Sequence_position and Sequence_name column titles, underscores to metals (Zn 66 becomes Zn_66),
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

#Create empty dataframe (Intensities_Zn) and start loop command
Zinc_analysis <- data.frame()
for (n in seq_pos) {
  #Get row number(idx) corresponding to sequence position
  idx <- which(Intensities$Sequence_position == n)
  
  #Retrieve RSD values from Zn_66 column in Intensities (one row below 'idx')
  df_rsd <- Intensities[idx+1, "Zn_66"]
  #Use function to remove % sign
  df_rsd <- lapply(df_rsd, function(x) { 
    gsub("%", "", x)
  })
  #Convert to data frame
  df_rsd <- as.data.frame(df_rsd, check.names = FALSE)
  colnames(df_rsd) <- "RSD_percent"
  
  #Use column bind (cbind) to assemble newly named columns into a new data frame (df)
  #combine Seq_pos, Sample name, Zn_66 intensity, RSD values
  df <- cbind(data.frame(Sequence_position = n, 
                         Sample_name = Intensities[idx,"Sample_name"]), 
              "Zn_66 (cps)" = Intensities$Zn_66[idx], 
              "RSD_percent" = df_rsd,
              "Sc_45 (cps)" = Intensities$Sc_45[idx])
  
  #At each iteration, stack df to new data frame (Intensities_Zn)
  Zinc_analysis <- rbind(Zinc_analysis, df)
}
#Look at resulting Zn dataframe with RSD as new column
View(Zinc_analysis)


#STEP THREE: Account for machine sampling variation using Sc_45 internal standard

#Convert Zn and Sc values from characters to numeric using as.numeric function.
#Add new columns to Zinc_analysis for: 
#Scfactor (all In values arbitrarily divided by 200000) 
#Scadjust (Zn values divided by Scfactor to adjust for amount of sample analyzed)
Zinc_analysis$Z01_Scfactor <- as.numeric(Zinc_analysis$Sc_45) / 200000
Zinc_analysis$Z02_Scadjust <- as.numeric(Zinc_analysis$Zn_66) / Zinc_analysis$Z01_Scfactor


#STEP FOUR: Continue Zn analysis by adjusting samples for baseline drift using 1 ppb standards

#Create data frame with only 1 ppb standard samples (Zinc_standards)
Zinc_standards <- filter(Zinc_analysis, Sample_name == "1 ppb standard")
#Use trendline function to create scatter plot for 1 ppb standards with trendline and equation
trendline(Zinc_standards$Sequence_position, Zinc_standards$Z02_Scadjust,model = "line2P", 
          Pvalue.corrected = TRUE, linecolor = "blue", lty = 1, lwd = 3, 
          show.equation = TRUE, show.Rpvalue = TRUE, Rname = 1, Pname = 0, xname = "x", yname = "y",
          yhat = FALSE, summary = TRUE, text.col = "black", main="Baseline Drift (1 ppb standards)", 
          xlab = "Sequence position", ylab = "Signal (Zn)", CI.fill = FALSE, CI.lty = "blank")


#STEP FIVE: Continue Zn analysis by plotting blanks and deducting blank average from all values

#Make data frame with only Blank samples (Zinc_blanks) and plot to check visually for outliers
Zinc_blanks <- filter(Zinc_analysis, Sample_name == "Blank")
plot(Zinc_blanks$Sequence_position, Zinc_blanks$Z02_Scadjust, 
     main="Trend of blank values", xlab = "Sequence position", ylab = "Signal (Zn)")
#Take average of all Blank values and use this to blank correct all other values
#This gives new column (Z04_Zn_blankadjust) in Zinc_analysis
Zinc_analysis$Z04_Zn_blankadjust <- Zinc_analysis$Z02_Scadjust - mean(Zinc_blanks[["Z02_Scadjust"]])


#STEP SIX: Plot Zn standard curve and use gradient to convert all sample Zn signals from cps to ppb

#Group (under Zinc_standard curve) then rename all standard curve samples and convert to numeric
Zinc_standardcurve <- filter(Zinc_analysis, Sample_name == "0.0004 ppb metal standard" 
                             | Sample_name == "0.002 ppb metal standard" | Sample_name == "0.02 ppb metal standard"
                             | Sample_name == "0.1 ppb metal standard" | Sample_name == "0.5 ppb metal standard"
                             | Sample_name == "1 ppb metal standard" | Sample_name == "5 ppb metal standard" 
                             | Sample_name == "10 ppb metal standard")
Zinc_stdconccharacters <- gsub(" ppb metal standard", "", Zinc_standardcurve$Sample_name) 
Zinc_standardcurve$Concentration <- as.numeric(Zinc_stdconccharacters)
#Plot trendline and use gradient to convert signal to ppb for all samples
trendline(Zinc_standardcurve$Concentration, Zinc_standardcurve$Z04_Zn_blankadjust,model = "line2P", 
          Pvalue.corrected = TRUE, linecolor = "blue", lty = 1, lwd = 3, 
          show.equation = TRUE, show.Rpvalue = TRUE, Rname = 1, Pname = 0, xname = "x", yname = "y",
          yhat = FALSE, summary = TRUE, text.col = "black", main="Zn standard curve", 
          xlab = "Zn concentration (ppb)", ylab = "Signal", CI.fill = FALSE, CI.lty = "blank")
Zinc_stdplot <- lm(Zinc_standardcurve$Z04_Zn_blankadjust ~ Zinc_standardcurve$Concentration)
#This gives new column (Z05_Zn_ppb) in Zinc_analysis
Zinc_analysis$Z05_Zn_ppb <- Zinc_analysis$Z04_Zn_blankadjust / coef(Zinc_stdplot)[2]


#STEP SEVEN: Calculate Zn concentration in original samples by accounting for dilution during sample prep/digestion

#Make Zinc_analysis2 dataframe as subset of Zinc_analysis dataframe, matching number of rows to Sample_prep.csv
Zinc_analysis2 <- Zinc_analysis %>% 
        filter(Sample_name != "0.0004 ppb metal standard" & Sample_name != "0.002 ppb metal standard"
           & Sample_name != "0.02 ppb metal standard" & Sample_name != "0.1 ppb metal standard" 
           & Sample_name != "0.5 ppb metal standard" & Sample_name != "1 ppb metal standard" 
           & Sample_name != "5 ppb metal standard" & Sample_name != "10 ppb metal standard" 
           & Sample_name != "Blank" & Sample_name != "1 ppb standard") %>%
        select(Sample_name, Z05_Zn_ppb, RSD_percent)

#Merge Sample_prep.csv data into Zinc_analysis2, match rows by Sample_name
Zinc_analysis2 <- merge(Zinc_analysis2, 
                           Sample_prep[,c("Sample_name", "standard_factor_digestion", "Sample_dilution_factor")], 
                           by = "Sample_name")

#Wetashcorrect (column Z06) - Adjust metal concentration for standard dilution during Wet Ash digestion (Dry Ash = 1X)  
Zinc_analysis2$Z06_Wetashcorrect <- Zinc_analysis2$Z05_Zn_ppb / Zinc_analysis2$standard_factor_digestion
#Sample_Zn_ppb (column Z07) - Adjust sample Zn values by dilution factor from digestion and sample prep
Zinc_analysis2$Z07_Sample_Zn_ppb <- Zinc_analysis2$Z06_Wetashcorrect / Zinc_analysis2$Sample_dilution_factor


#STEP EIGHT: Apply RSDs and summarise results for Zn (original sample and SD in ppb, detection limit)

#Multiply all Sample_Zn_ppb values by numeric equivalent of (RSD/100) to get error
Zinc_analysis2$Z08Zn_RSD_ppb <- Zinc_analysis2$Z07_Sample_Zn_ppb * (as.numeric(Zinc_analysis2$RSD_percent) / 100)

#Use 3 * STD of Blanks to get detection limit (in ppb) to appear below results
Zinc_detectionlimit_ppb <- (3 * sd(Zinc_blanks[["Z02_Scadjust"]])) / coef(Zinc_stdplot)[2]

#Arrange results summary in new data frame (Zinc_results)
#Include Sample_name, Zn_ppb detected in measured samples, Zn_ppb concentration in original sample, RSD_percent, RSD_ppb
Zinc_results <- select(Zinc_analysis2, Sample_name, Z07_Sample_Zn_ppb, RSD_percent, Z08Zn_RSD_ppb, Z05_Zn_ppb)

#Rename column titles to give output of results in console
Zinc_results <- data.frame(Sample_name = Zinc_results$Sample_name, 
                              Zn_ppb_original_sample = Zinc_results$Z07_Sample_Zn_ppb,
                              RSD_percent = Zinc_results$RSD_percent,
                              RSD_ppb = Zinc_results$Z08Zn_RSD_ppb, 
                              Zn_ppb_detected_in_measured_sample = Zinc_results$Z05_Zn_ppb, 
                              Detection_limit_ppb = Zinc_detectionlimit_ppb)
View(Zinc_results)


#STEP NINE: Save data to file

#Write results to file, first creating Results folder in working directory
if(!dir.exists("Results")) dir.create("Results")
#Make csv file of Zinc_results dataframe
write.csv(Zinc_results, file = file.path("Results", "Zinc_results.csv"), row.names = FALSE)


#######################################
#END OF DATA ANALYSIS
#######################################


#View step by step data analysis
View(Zinc_analysis)
View(Zinc_analysis2)
#View raw data (csv files) to bring dataframes up in tabs and check them
View(Intensities)
View(Sample_prep)




#Clear plots
dev.off()
#Clear environment
rm(list = ls())
#Clear console
cat("\f")