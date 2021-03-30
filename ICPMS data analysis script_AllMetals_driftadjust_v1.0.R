#Script for ICP-MS data analysis of all metals (Latest version of RStudio recommended)

#Update working directory, file names, and sample names (change to "1 ppb standard" and correct metal standard concentrations)

#R SET UP

#install then load tidyverse (contains all packages needed such as dplyer and ggplot2) and basicTrendline
#Make sure these are loaded each time this is used
if(!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)
if(!require("basicTrendline")) install.packages("basicTrendline")
library(basicTrendline)

#Set working directory - change user name and file path as required
setwd("C:/Users/tomasrs/Desktop/R working directory ICPMS")



#LOAD .csv FILES CORRESPONDING TO DATA INTO RStudio

#Create data frames for all .csv files (Sample_run_intensities_XYZ, Sample_digest_XYZ)
#Change file names in the script below as required
#Files must be in same folder as working directory/project/script
#Read as csv (read.csv) or excel file (read_excel from readxl package) 
#View csv files to bring dataframes up in tabs

#Sample_prep csv file must be prepared to ensure that order of samples matches ICPMS sequence
#Prepare Sample_run_intensities csv files by adding:
#...Sequence_position and Sequence_name column titles, underscores to metals (Fe 56 becomes Fe_56),
#...and removing columns for acquisition time, dataset file and method file
#...and deleting all worksheets except the one with both Intensities and RSDs
#Save this worksheet as a .csv file (Sample_run_XYZ)

Intensities <- read.csv("Sample_run_test.csv", stringsAsFactors = FALSE)
View(Intensities)
Sample_prep <- read.csv("Sample_digest_test.csv", stringsAsFactors = FALSE)
view(Sample_prep)

#BEGIN DATA ANALYSIS

#STEP ONE: Reshape intensities data frame to 'long' format

#Sequence position,Sample names, metal, measure (cps or rsd) and value 
#Split intensities into csp measures and RSD values

#Set up some data frames in preparation for a 'for loop'
#Make new data frame (Intensities_names) for the column names in 'Intensities'
Intensities_names <- colnames(Intensities)
#Make another data frame (Metals_names) by excluding the column names 'Sequence_position' and 'Sample_name'
Metals_names <- Intensities_names[!Intensities_names %in% c("Sequence_position", "Sample_name")]

#Keep track of cps and rsd measures by making new data frames with cps and rsd appended to names
cps_names <- paste0("cps ", Metals_names)
RSD_names <- paste0("RSD ", Metals_names)
#Retrieve sequences positions without 'NA's in new data frame (seq_pos)
seq_pos <- Intensities$Sequence_position
seq_pos <- seq_pos[!is.na(seq_pos)]

#Create 'for loop' through each sequence position to create a new data.frame (df)
df <- data.frame()
for (n in seq_pos) {
  #Get row number(idx) corresponding to sequence 
  idx <- which(Intensities$Sequence_position == n)
  #Retrieve cps values from Intensities (character), transform them to numeric in a data frame
  df_cps <- Intensities[idx, Metals_names]
  df_cps <- lapply(df_cps, as.numeric)
  df_cps <- as.data.frame(df_cps, check.names = FALSE)
  #Make column names 'cps metal' (from cps_names data frame)
  colnames(df_cps) <- cps_names
 
  #Retrieve RSD values from Intensities (one row below 'idx')
  df_rsd <- Intensities[idx+1, Metals_names]
  #Make new function to remove % sign then transform them to numeric and divide by 100 (5% becomes 0.05)
  df_rsd <- lapply(df_rsd, function(x) { 
    rsd = gsub("%", "", x)
    as.numeric(rsd)/100})
  #Convert to data frame
  df_rsd <- as.data.frame(df_rsd, check.names = FALSE)
  #Make column names 'RSD metal' (from RSD_names data frame)
  colnames(df_rsd) <- RSD_names
  
  #Use column bind (cbind) to assemble cps and rsd all in one row of a new data frame (df_row)
  df_row <- cbind(data.frame(Sequence_position = n, 
                             Sample_name = Intensities[idx,"Sample_name"]), 
                  df_cps, 
                  df_rsd)
  
  #At each iteration, stack df_row to df
  df <- rbind(df, df_row)
}

#Transform df to long format
#Gather makes one long data frame with column "cps/RSD metal"
df_long <- gather(df, measure, value, c(cps_names, RSD_names)) %>%
  separate(measure, c("measure", "metal"), sep = " ")
#Separate splits "cps/RSD metal" column into "measure" (cps or RSD) and "metal" (metal) due to space

#Create two separate data.frame for cps and rsd
#This means each data frame can be accessed 
df_cps <- filter(df_long, measure == "cps")
df_rsd <- filter(df_long, measure == "RSD")

#Get all sample names using 'unique' then remove the blanks and standards
sample_names <- unique(df_cps$Sample_name)
sample_names <- sample_names[!sample_names %in% c("0.1 ppb metal standard", 
                                                  "1 ppb metal standard", 
                                                  "5 ppb metal standard", 
                                                  "10 ppb metal standard", 
                                                  "50 ppb metal standard", 
                                                  "100 ppb metal standard", 
                                                  "400 ppb metal standard", 
                                                  "Blank", 
                                                  "1 ppb standard")]


#STEP TWO: Account for machine sampling variation using Sc internal standard

#Divide all Sc intensities by 50000 (arbitrary number) 
#Divide all metals cps by Sc intensity to adjust for amount of sample analyzed (Sc factor)
#Add new columns to df_cps for: 
#Scfactor (all Sc values arbitrarily divided by 50000) 
#Scadjust (metal values divided by Scfactor to adjust for amount of sample that was analyzed)

Scfactor <- df_cps %>%
  group_by(Sequence_position) %>% #Sequence_position is unique to each sample
  filter(metal == "Sc_45") %>%
  mutate(Scfactor = value / 50000)

#We have one adjust factor per sample (seq position)
#We use merge to combine df_csp and Scfactor
df_cps <- merge(df_cps, Scfactor[,c("Sequence_position", "Scfactor")], by = "Sequence_position")
df_cps$Scadjust <- df_cps$value / df_cps$Scfactor


#Start a for loop going through each metal

for (analyzedMetal in Metals_names) {
  
  #Filter for analyzed metal
  df_cps_metal <- df_cps %>%
    filter(metal == analyzedMetal)
  
  #STEP THREE: Adjust samples for baseline drift using 1 ppb standards
  
  #Create data frame (df_standards) with only 1 ppb standard samples
  df_standards <- df_cps_metal %>%
    filter(Sample_name == "1 ppb standard")
 
  #Use trendline function to create scatter plot for 1 ppb standards with trendline and equation
  trendline(df_standards$Sequence_position, df_standards$Scadjust,model = "line2P", 
            Pvalue.corrected = TRUE, linecolor = "blue", lty = 1, lwd = 3, 
            show.equation = TRUE, show.Rpvalue = TRUE, Rname = 1, Pname = 0, xname = "x", yname = "y",
            yhat = FALSE, summary = TRUE, text.col = "black", main=paste0(analyzedMetal, " Baseline Drift (1 ppb standards)"), 
            xlab = "Sequence position", ylab = paste0("Signal (",analyzedMetal,")"), CI.fill = FALSE, CI.lty = "blank")
  
  #lm function allows use of m (gradient, coef[2]) in equation to apply drift adjustment (adjust_drift)
  lm_drift <- lm(df_standards$Scadjust ~ df_standards$Sequence_position)
  
  df_cps_metal$adjust_drift <- df_cps_metal$Scadjust - coef(lm_drift)[2] * (df_cps_metal$Sequence_position - df_standards$Sequence_position[1])
  
  #STEP FOUR: Plot blanks and deduct blank average from all values
  
  #Make data frame with only Blank samples (df_blanks) and plot to check visually for outliers
  df_blanks <- filter(df_cps_metal, Sample_name == "Blank")
  
  plot(df_blanks$Sequence_position, df_blanks$adjust_drift, 
       main=paste0(analyzedMetal, " Trend of blank values"), xlab = "Sequence position", ylab = paste0("Signal (",analyzedMetal,")"))
  
  #Take average of all Blank values and use this to blank correct all other values (adjust_blank)
  df_cps_metal$adjust_blank <- df_cps_metal$adjust_drift - mean(df_blanks$adjust_drift)
  

  
  #STEP FIVE: Plot metal standard curve and use gradient to convert all sample metal signals to ppb
  
  #Group then rename all standard curve samples and convert to numeric (0.1 ppb metal standard becomes 0.1)
  metal_standardcurve <- filter(df_cps_metal, Sample_name %in% c("0.1 ppb metal standard", 
                                                                 "1 ppb metal standard", 
                                                                 "5 ppb metal standard", 
                                                                 "10 ppb metal standard", 
                                                                 "50 ppb metal standard", 
                                                                 "100 ppb metal standard", 
                                                                 "400 ppb metal standard"))
  
  metal_stdconccharacters <- gsub(" ppb metal standard", "", metal_standardcurve$Sample_name) 
  metal_standardcurve$Concentration <- as.numeric(metal_stdconccharacters)
  
  #Plot trendline and use gradient to convert signal to ppb for all samples
  trendline(metal_standardcurve$Concentration, metal_standardcurve$adjust_blank,model = "line2P", 
            Pvalue.corrected = TRUE, linecolor = "blue", lty = 1, lwd = 3, 
            show.equation = TRUE, show.Rpvalue = TRUE, Rname = 1, Pname = 0, xname = "x", yname = "y",
            yhat = FALSE, summary = TRUE, text.col = "black", main=paste0(analyzedMetal, " standard curve"), 
            xlab = "Fe concentration (ppb)", ylab = "Signal", CI.fill = FALSE, CI.lty = "blank")
  
  metal_stdplot <- lm(adjust_blank ~ Concentration, data = metal_standardcurve)
  df_cps_metal$ppb <- df_cps_metal$adjust_blank / coef(metal_stdplot)[2]
  
  
  
  #STEP SIX: Calculate metal concentration in samples by accounting for dilution during sample prep/digestion
  
  #Make df_cps_metal_samples dataframe by merging Sample_prep.csv data into df_cps_metal, match rows by Sample_name

  df_cps_metal_samples <- merge(df_cps_metal, Sample_prep[,c("Sample_name","Sc_factor_digestion", "Sample_dilution_factor")], by = "Sample_name")
  
  #Wetashcorrect - Adjust metal concentration for Sc dilution during Wet Ash digestion (Dry Ash = 1X)  
  df_cps_metal_samples$Wetashcorrect <- df_cps_metal_samples$ppb / df_cps_metal_samples$Sc_factor_digestion
 
  #Sample_ppb - Adjust metal ppb values by dilution factor from digestion and sample prep
  df_cps_metal_samples$Sample_ppb <- df_cps_metal_samples$Wetashcorrect / df_cps_metal_samples$Sample_dilution_factor
  
  
  
  #STEP SEVEN: Apply RSDs and summarise results for metal (original sample and RSD in ppb, detection limit)
  
  
  #Make df_rsd_metal dataframe as subset of df_rsd data frame
  df_rsd_metal <- filter(df_rsd, metal == analyzedMetal)
  df_rsd_metal_samples <- df_rsd_metal %>% 
    filter(Sample_name %in% sample_names)
  
  #Need to make sure that samples are in the same order between df_rsd and df_cps
  #Calculate RSD for samples by multiplying % error (df_rsd_metal_samples) by sample value (Sample_ppb)
  df_cps_metal_samples$SD_ppb <- df_cps_metal_samples$Sample_ppb * (df_rsd_metal_samples$value)
  
  #Use 3 * STD of Blanks to get detection limit for metal
  metal_detectionlimit <- (3 * sd(df_blanks$adjust_drift)) / coef(metal_stdplot)[2]
  
  #Arrange results summary (Sample_name, Sample_ppb in original sample, SD_ppb) in new data frame (metal_results)
  metal_results <- select(df_cps_metal_samples, Sample_name, 
                          Sample_ppb, SD_ppb)
  
  #Rename column titles to give output of results in console
  metal_results <- data.frame(Sample_name = df_cps_metal_samples$Sample_name, 
                              metal = analyzedMetal, 
                              Concentration = df_cps_metal_samples$Sample_ppb, 
                              RSD = df_cps_metal_samples$SD_ppb, metal_detectionlimit)
  

  
  #STEP EIGHT: Make plots and save data to file
  
  plot <- ggplot(metal_results, aes(x = Sample_name, y = Concentration)) +  
    geom_point() + 
    geom_errorbar(aes(ymin = Concentration - RSD,
                      ymax = Concentration + RSD),
                  width = 0.5,
                  position=position_dodge(0.05)) +
    geom_abline(intercept = metal_detectionlimit, slope = 0, color = "red", linetype = "dashed") +
    ylab("Concentration (ppb)") + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.text.x = element_text(angle = - 90, hjust = 0), 
          axis.title.x = element_blank())
  
  plot
  
  plot_log10 <- ggplot(metal_results, aes(x = Sample_name, y = Concentration)) + 
    geom_point() + 
    geom_errorbar(aes(ymin = Concentration - RSD,
                      ymax = Concentration + RSD),
                  width = 0.5,
                  position=position_dodge(0.05)) +
    geom_abline(intercept = metal_detectionlimit, slope = 0, color = "red", linetype = "dashed") + 
    scale_y_log10() + 
    ylab("Concentration (ppb)") + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.text.x = element_text(angle = - 90, hjust = 0), 
          axis.title.x = element_blank())
  
  plot_log10
  
  # Write results and figure to file 
  
  # Change column names
  colnames(metal_results) <- c("Sample", 
                               "metal",
                               "Concentration (ppb)", 
                               "RSD (± ppb)")
  
  if(!dir.exists("Results")) dir.create("Results")
  
  write.csv(metal_results, file = file.path("Results", paste0(analyzedMetal,"_results.csv")), row.names = FALSE)
  write.csv(metal_detectionlimit, file = file.path("Results", paste0(analyzedMetal,"_detectionlimit.csv")), row.names = FALSE)
  ggsave(plot_log10, filename = file.path("Results", paste0(analyzedMetal,"_log10.png")), units = "in", 
         width = 5, height = 4, dpi = 300)
  ggsave(plot, filename = file.path("Results", paste0(analyzedMetal,"_linear.png")), units = "in", 
         width = 5, height = 4, dpi = 300)
}

#END OF ANALYSIS



#Clear plots
dev.off()
#Clear environment
rm(list = ls())
#Clear console
cat("\f")