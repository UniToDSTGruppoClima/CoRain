###############################################################################
############################################################################### 
########################         Co.Rain         ##############################
######################## Comparing series of Rain #############################
###############################################################################
###############################################################################
###                                                                         ###
### Copyright (C) 2015                                                      ###
### Fiorella Acquaotta, Diego Garzena, Diego Guenzi and Simona Fratianni    ###
###                                                                         ###
### This program is free software: you can redistribute it and/or modify    ###
### it under the terms of the GNU General Public License as published by    ###
### the Free Software Foundation, either version 3 of the License, or       ###
### (at your option) any later version.                                     ###
###                                                                         ###
### This program is distributed in the hope that it will be useful,         ###
### but WITHOUT ANY WARRANTY; without even the implied warranty of          ###
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           ###
### GNU General Public License for more details.                            ###
###                                                                         ###
### You should have received a copy of the GNU General Public License       ###
### along with this program.  If not, see <http://www.gnu.org/licenses/>    ###
###                                                                         ###
###############################################################################
###############################################################################


###############################  Co.Rain  #####################################
# Description of the software
#
# This program works in three steps: statistical analysis, comparison between 
# the cleaned series (with same missing values and daily values <= 1mm) and 
# comparison between precipitation classes. For major details, see the paper 
# by Acquaotta F., Fratianni S. and Venema V. (2016) - Assessment of parallel 
# precipitation measurements networks in Piedmont, Italy (International 
# Journal of Climatology - DOI 10.1002/joc.4606).
# 1] Statistical analysis
# A statistical analysis is carried out on the two series, the candidate 
# series (column 4 of the input file) and the reference one (column 5). The 
# analysis starts with the two raw series, removes all values with daily 
# rain < 1mm and then puts the same missing values on both series, cleaning 
# them. The results are reported in the file "0_statistics_input_file.txt". 
# The statistical analysis calculates the minimum value, the 1st quantile, 
# the median, the mean, the 3rd quantile, the maximum value, the number of 
# missing values, the number of values and the results (W and p_value) of 
# Shapiro-Wilk test. In the 2 files that end with "_daily_avalibillity.csv" 
# are written the monthly number of days used in this analysis; in the 2 
# files that end with "_monthly_rain.csv" you can find the monthly amount of 
# rain; at the end, in the 2 files with "_day_month_rain.png" are reported 
# the plots of daily and monthly precipitation series with same missing 
# values and daily values <= 1mm.
# 2] Comparison between the cleaned series (with same missing values and 
# daily values <= 1mm)
# A statistical comparison analysis between the two cleaned series is 
# carried out. We apply the Student test (t), the Kolmogorov-Smirmov test 
# (ks), the Wilcoxon rank sum test (W) and the Kruskal test (kru), in 
# addition to the correlation coefficient (cor) with Spearman's method 
# and the Root Mean Square Error (RMSE). The results of the comparison 
# analysis are reported in the file "4_statistics_between_daily_series.csv" 
# where, for every test, we report the statistical value and the p_value. 
# A scatter plot, "5_total_precipitation_plot_with_15_5.png", with daily 
# values is also created. The maximum and minimum error is calculated for 
# the reference series. The maximum error (15) is equal to ± 15% of daily 
# value while the minimum error (5) is equal to ± 5% of daily value. The 
# monthly values is calculated for the two cleaned series and, after that, 
# the percentage relative error is computed; the results are reported in 
# the file "6_percentage_relative_error.csv", with a box plot called 
# "7_percentage_relative_error_boxplot.png". On the percentage relative 
# series, the non-parametric trends are elaborated using the Theil-Sen 
# approach (TSA). The Mann-Kendall test for the trend is then run on the 
# resulting time series to compute the level of significance. The 
# results of the trend are written in the file 
# "8_stats_percentage_relative_error_trend.csv", where we report the 
# trend, the intercept, the trend over all period, the tau of Mann 
# Kendall test and the p_value. After this, if some values of percentage 
# relative error are outside a specific range (set at ±500 by default), 
# we start a re-computation of the only values inside the range, writing 
# results in a file called "6_filtered_percentage_relative_error.csv", 
# re-plotting data in the file 
# "7_filtered_percentage_relative_error_boxplot.png" and re-calculating 
# the trend; the results of the trend are written in the file 
# "8_filtered_stats_percentage_relative_error_trend.csv" and plotted in 
# "8_filtered_percentage_relative_error_trend.png".
# 3] Comparison between five classes of precipitation events
# In this final step we calculate the quantile of the reference series 
# to identify the thresholds of the class. Every daily rain event is 
# classified as weak, mean, heavy, very heavy (R95) or extreme (R99). 
# For each class we create a scatterplot: "9_events_weak.png", 
# "9_events_mean.png", "9_events_heavy.png", "9_events_R95.png" and 
# "9_events_R99.png". In the scatterplots are reported the maximum 
# error, ± 15% of daily value, and the minimum error, ± 5% of daily 
# value. In the file "10_class_events_and_RMSE.csv" for every class, 
# we report: the thresholds, the mean value of maximum and minimum 
# error, number of events for reference and candidate series, 
# precipitation amount for reference and candidate series, number 
# of events recorded in the same day, precipitation amount for 
# reference and candidate series recorded in the same day, number 
# of events outside the range for maximum and minimum error and RMSE.
#
#
# When you start the program, (if interactive mode is available) you will be 
# asked for an input text file and a folder where results will be stored. If
# you are running in batch mode, please modify the code to use the correct 
# input file and create the right output folder.
#
# Please, note that a step of cleaning is made on the input file. The file 
# used as starting point of the program and its statistics are available in 
# output number 0 and 1. If these doesn't suit your needs, please discard 
# the whole computation.
#
# INPUT
# The text file has to be formatted in five TAB-separated columns. The first
# row of the file has to contain the headers (column names) and the first
# three columns have to contain the year, the month and the day of the
# series. Column four is the candidate rain serie and column five is the
# reference rain serie. Missing values must be marked as NA. It is very
# important to start the file from the first of January (of any year) and
# end it at the 31st of December (of any year) with at least 5 years of
# data. See the attached file called example.txt
#
# OUTPUT
# 0_statistics_input_file.csv - Important statistics on input file and 
#                               cleaned input file
# 1_cleaned_input_file.txt - The file that is really used as input
# 2_can_Pidro_daily_availability.csv - Number of days of available data 
#                                      on candidate serie
# 2_can_Pidro_monthly_rain.csv - Total amount of rain by month on 
#                                candidate serie
# 2_ref_Parpa_daily_availability.csv - Number of days of available data 
#                                      on reference serie
# 2_ref_Parpa_monthly_rain.csv - Total amount of rain by month on 
#                                reference serie
# 3_can_Pidro_day_month_rain.png - Plot of daily and monthly data on 
#                                  candidate serie
# 3_ref_Parpa_day_month_rain.png - Plot of daily and monthly data on 
#                                  reference serie
# 4_statistics_between_daily_series.csv - RMSE, T, KS, Wilcox, 
#                                         Kruskal and Spearman between 
#                                         the two daily series
# 5_total_precipit_plot_with_15_5.png - Scatter plot of daily events
# 6_percentage_relative_error.csv - Percentage of relative errors 
#                                   grouped by month (filtered optional)
# 7_percentage_relative_error_boxplot.png - Box plot of previous data 
#                                           (filtered optional)
# 8_percentage_relative_error_trend.png - Trend of previous data 
#                                         (filtered optional)
# 8_stats_percentage_relative_error_trend.csv - Statistics on trend 
#                                               (filtered optional)
# 9_events_R99.png - Scatter plot of R99 precipitation
# 9_events_R95.png - Scatter plot of R95 precipitation
# 9_events_heavy.png - Scatter plot of heavy precipitation
# 9_events_mean.png - Scatter plot of mean precipitation
# 9_events_weak.png - Scatter plot of weak precipitation
# 10_class_event_and_RMSE.csv - Summary of information on classes
#
# Following R-packages have to be installed:
#    "class", "zoo", "hydroGOF", "xts", "hydroTSM", "zyp", "Kendall"
#
# This code has been written under R-Version 3.2.2; for older or newer
# versions problems might occur.
#
###############################################################################
# Versions:
# 
# v0.9 - 20151216: First creation (see publication from Acquaotta et al., 2015
#                  Assessment of parallel precipitation measurements networks
#                  in Piedmont, Italy)
# v1.0 - 20160128: First public release after code cleaning and review
# 
###############################################################################


###############################################################################
####                         INITIALIZATION                                ####
###############################################################################

start.time = Sys.time()
require(class)
require(zoo)
require(hydroGOF)
require(xts)
require(hydroTSM)
require(zyp)
require(Kendall)

###############################################################################
####                       FUNCTION DEFINITION                             ####
###############################################################################

serie_stat_analysis <- function(name, serie) {
  ##########################################################
  # Compute statistical analysis on the serie, plot results 
  # and write tables
  #
  # INPUT
  # name: the suffix of the file names generated by this function
  # serie: a single rain serie to analyze
  #
  # OUTPUT
  # Matrix containing monthly informations on rain
  # This function also saves three files: two csv and a png plot
  ##########################################################
  
  serie.zoo = zooreg(serie, start=as.Date(paste(y.start, "-01-01", sep="")))
  num.serie = dwi(serie.zoo, out.unit="mpy")
  
  # Save first CSV with the number of wet days
  Y.num.serie = cbind(years, num.serie)
  write.table(Y.num.serie, file=paste("2_",name,"_daily_availability.csv",sep=""), 
              sep=";", row.names=FALSE)
  
  # Save PNG plot
  png(paste("3_",name,"_day_month_rain.png",sep=""))
  hydroplot(serie.zoo, na.rm=TRUE, var.type="Precipitation", pfreq="dm", 
            ptype="ts", ylab=name)
  dev.off()
  
  # Monthly elaborations
  m.serie = daily2monthly(serie.zoo, FUN=sum, na.rm=TRUE) 
  monthly = matrix(m.serie, ncol=12, byrow=TRUE)
  
  # Save second CSV for monthly information
  Y.monthly = cbind(years, monthly)
  colnames(Y.monthly) = c("Year", months)
  write.table(Y.monthly, file=paste("2_",name,"_monthly_rain.csv",sep=""), 
              sep=";", row.names=FALSE)
  
  return(monthly)
} # end of function serie_stat_analysis


class_formatting <- function(my_class) {
  ##########################################################
  # Formats the class passed as input
  #
  # INPUT
  # my_class: class to be formatted
  #
  # OUTPUT
  # Dataframe containing the formatted class
  ##########################################################
  
  # Insert zero before single digit months and days
  my_class[,2] = sprintf("%02d", my_class[,2])
  my_class[,3] = sprintf("%02d", my_class[,3])
  
  # Create yyyymmdd column
  my_data = paste(my_class[,1], my_class[,2], my_class[,3])
  my_data = as.matrix(my_data, ncol=1)
  my_data_def = cbind(my_data, my_class[,4])
  
  # Transform into data frame
  df_data = as.data.frame(my_data_def)
  colnames(df_data) = c("data", "val")
  
  return(df_data)
} # end of function class_formatting


class_analysis <- function(join_c, join_r, name) {
  ##########################################################
  # Computes statistics on the common days of the two series 
  # in input and plots the result
  #
  # INPUT
  # join_c: candidate serie to be joined
  # join_r: reference serie to be joined
  # name: class used for computations
  #
  # OUTPUT
  # Matrix containing the values
  # This function also saves a png plot
  ##########################################################
  
  # Merge the series by data
  join = merge(join_r, join_c, by="data") 
  join = as.matrix(join)
  scd = join[,-1]
  scd = matrix(as.numeric(scd), nrow=nrow(scd), ncol=2)
  progressive = c(1:nrow(scd))
  scd = cbind(progressive,scd)
  
  # 15% and 5% error on daily data from the joined ref serie
  er_p_15 = ((scd[,2] * 15) / 100) + scd[,2]
  er_n_15 = (-1 * (scd[,2] * 15) / 100) + scd[,2]
  er_p_5 = ((scd[,2] * 5) / 100) + scd[,2]
  er_n_5 = (-1 * (scd[,2] * 5) / 100) + scd[,2]
  
  # Mean errors from the joined ref serie
  m_scd_15 = mean((scd[,2] * 15) / 100)
  m_scd_5 = mean((scd[,2] * 5) / 100)
  m_scd = mean(scd[,2])
  
  # Number of common events in the 15% range
  new_scd_15 = cbind(scd, er_p_15, er_n_15)
  max_15 = new_scd_15[,3]<=new_scd_15[,4]
  min_15 = new_scd_15[,3]>=new_scd_15[,5]
  eve_15 = max_15&min_15
  pm_15 = t(matrix(summary(eve_15)))
  pmn_scd_15 = as.numeric(pm_15[,2])
  
  # Number of common events in the 5% range
  new_scd_5 = cbind(scd, er_p_5, er_n_5)
  max_5 = new_scd_5[,3]<=new_scd_5[,4]
  min_5 = new_scd_5[,3]>=new_scd_5[,5]
  eve_5 = max_5&min_5
  pm_5 = t(matrix(summary(eve_5)))
  pmn_scd_5 = as.numeric(pm_5[,2])
  
  # Save plot
  png(paste("9_events_",name,".png",sep=""))
  plot(scd[,2], scd[,2], main=paste(name,"precipitation",sep=" "), 
       ylab="[mm]", xlab="[mm]")
  lines(scd[,2], er_p_15)
  lines(scd[,2], er_n_15)
  lines(scd[,2], er_p_5, col="green")
  lines(scd[,2], er_n_5, col="green")
  lines(scd[,2], scd[,3], type="p", col="red")
  legend('topleft', lty=1, col=c('red','black','green','black'), 
         legend=c("(circles) Candidate values","(circles) Reference values",
                  "(lines) 5% error","(lines) 15% error"), bty='n', cex=.75)
  dev.off()
  
  return(as.matrix(data.frame(scd, m_scd_15, m_scd_5, m_scd, pmn_scd_15, pmn_scd_5)))
} # end of function class_analysis

###############################################################################
####                              MAIN PROGRAM                             ####
###############################################################################

# Read input data and output folder creation
if (interactive()) { # User has to choose file in input and path of results
  tab = read.table(file.choose(), header=T, na.strings="NA")
  setwd(choose.dir(getwd(), "Choose a folder to save your results:"))
} else { # File in input and path of results are hardcoded
  setwd("/data/test") # Path where input file is located
  tab = read.table("example.txt", header=T, na.strings="NA") # Input file name
  dir.create("./results/", showWarnings = FALSE) # Results folder
  setwd("./results/") # Results folder
}

# Prepares series, years, months and other global variables
treshold = 500 # Set treshold for percentage relative error
can.serie = tab[,1:4]
ref.serie = tab[,-4]
y.start = as.numeric(head(tab[1], 1))
y.end = as.numeric(tail(tab[1], 1))
years = c(y.start:y.end)
months = c(1:12)
months = month.abb[months]

####    STATISTICAL ANALYSIS   ####

# Obtain stats about useful data in raw series
stats.can.raw = as.vector(summary(can.serie[,4]))
stats.can.raw[8] = length(which(!is.na(can.serie[,4])))
stats.ref.raw = as.vector(summary(ref.serie[,4]))
stats.ref.raw[8] = length(which(!is.na(ref.serie[,4])))
if (length(which(!is.na(can.serie[,4]))) > 5000) {
  # Shapiro-Wilk only on <5000 values
  stats.can.raw[9] = NA
  stats.can.raw[10] = NA
} else {
  test.shap.can = shapiro.test(can.serie[,4])
  stats.can.raw[9] = round(test.shap.can$statistic, digits=2)
  stats.can.raw[10] = round(test.shap.can$p.value, digits=5)
}
if (length(which(!is.na(ref.serie[,4]))) > 5000) {
  # Shapiro-Wilk only on <5000 values
  stats.ref.raw[9] = NA
  stats.ref.raw[10] = NA
} else {
  test.shap.ref = shapiro.test(ref.serie[,4])
  stats.ref.raw[9] = round(test.shap.ref$statistic, digits=2)
  stats.ref.raw[10] = round(test.shap.ref$p.value, digits=5)
}
stat.name = c("Min", "1st quantile", "Median", "Mean", "3rd quantile", 
              "Max", "Number of NA", "Number of values", 
              "Shapiro-Wilk normality test W", 
              "Shapiro-Wilk normality test p-value")
stats = cbind(stat.name, stats.can.raw, stats.ref.raw)

# Set as NA all precipitations < 1mm, ignoring them on both series
can.serie[can.serie<1] = NA
ref.serie[ref.serie<1] = NA

# Obtain stats about useful data in partially cleaned series
stats.can.partCleaned = as.vector(summary(can.serie[,4]))
stats.can.partCleaned[8] = length(which(!is.na(can.serie[,4])))
stats.ref.partCleaned = as.vector(summary(ref.serie[,4]))
stats.ref.partCleaned[8] = length(which(!is.na(ref.serie[,4])))
if (length(which(!is.na(can.serie[,4]))) > 5000) {
  # Shapiro-Wilk only on <5000 values
  stats.can.partCleaned[9] = NA
  stats.can.partCleaned[10] = NA
} else {
  test.shap.can = shapiro.test(can.serie[,4])
  stats.can.partCleaned[9] = round(test.shap.can$statistic, digits=2)
  stats.can.partCleaned[10] = round(test.shap.can$p.value, digits=5)
}
if (length(which(!is.na(ref.serie[,4]))) > 5000) {
  # Shapiro-Wilk only on <5000 values
  stats.ref.partCleaned[9] = NA
  stats.ref.partCleaned[10] = NA
} else {
  test.shap.ref = shapiro.test(ref.serie[,4])
  stats.ref.partCleaned[9] = round(test.shap.ref$statistic, digits=2)
  stats.ref.partCleaned[10] = round(test.shap.ref$p.value, digits=5)
}
stats = cbind(stats, stats.can.partCleaned, stats.ref.partCleaned)

# Transform the series having the same missing data 
can.serie[is.na(ref.serie)] = NA
ref.serie[is.na(can.serie)] = NA

# Obtain and save stats about useful data in cleaned series
stats.can.cleaned = as.vector(summary(can.serie[,4]))
stats.can.cleaned[8] = length(which(!is.na(can.serie[,4])))
stats.ref.cleaned = as.vector(summary(ref.serie[,4]))
stats.ref.cleaned[8] = length(which(!is.na(ref.serie[,4])))
if (length(which(!is.na(can.serie[,4]))) > 5000) {
  # Shapiro-Wilk only on <5000 values
  stats.can.cleaned[9] = NA
  stats.can.cleaned[10] = NA
} else {
  test.shap.can = shapiro.test(can.serie[,4])
  stats.can.cleaned[9] = round(test.shap.can$statistic, digits=2)
  stats.can.cleaned[10] = round(test.shap.can$p.value, digits=5)
}
if (length(which(!is.na(ref.serie[,4]))) > 5000) {
  # Shapiro-Wilk only on <5000 values
  stats.ref.cleaned[9] = NA
  stats.ref.cleaned[10] = NA
} else {
  test.shap.ref = shapiro.test(ref.serie[,4])
  stats.ref.cleaned[9] = round(test.shap.ref$statistic, digits=2)
  stats.ref.cleaned[10] = round(test.shap.ref$p.value, digits=5)
}
stats = cbind(stats, stats.can.cleaned, stats.ref.cleaned)
write.table(stats, file="0_statistics_input_file.csv", sep=";", row.names=FALSE)

# Save cleaned input file
input_file = cbind(can.serie,ref.serie[4])
write.table(input_file, file="1_cleaned_input_file.txt", sep="\t", row.names=FALSE)

# Ask to the user if he wants to proceed
quit = FALSE
if (interactive()) {
  cat("Here you can find statistics about your data", "\n")
  print(stats)
  max_num_val = max(as.numeric(stats[8,2]), as.numeric(stats[8,3]))
  cat("Pay attention especially to this: we are using", stats[8,6], 
      "values of yours", max_num_val,"\n")
  go_ahead = FALSE
  while (!go_ahead) {
    cat("\n", "Do you want to proceed with the analysis? (Y/N)")
    yn = scan(n=1, what=character())
    yn = toupper(yn)
    if (yn=="N" || yn=="NO") {
      quit = TRUE
      go_ahead = TRUE
    } else {
      if (yn!="Y" && yn!="YES") {
        cat("\n", "Please, reply only Y or N (yes or no)")
      } else {
        go_ahead = TRUE
      }
    }
  }
}

if (!quit) { # User wants to proceed
  
  # Single serie analysis
  monthly_can = serie_stat_analysis(paste("can",colnames(tab)[4],sep="_"), 
                                    can.serie[,4])
  monthly_ref = serie_stat_analysis(paste("ref",colnames(tab)[5],sep="_"), 
                                    ref.serie[,4])
  
  # Calculates RMSE, T, KS, Wilcox, Kruskal and Cor between the two daily series
  error = rmse(ref.serie[,4], can.serie[,4], na.rm=T)
  rmse.can.ref = round(error, digits=2)
  test.t.can.ref = t.test(can.serie[,4], ref.serie[,4])
  test.ks.can.ref = suppressWarnings(ks.test(can.serie[,4], ref.serie[,4]))
  test.wil.can.ref = wilcox.test(can.serie[,4], ref.serie[,4])
  test.kru.can.ref = kruskal.test(can.serie[,4], ref.serie[,4])
  test.cor.can.ref = cor.test(can.serie[,4], ref.serie[,4], method="spearman", 
                              exact=FALSE)
  test.stats.name = c("RMSE", "T test t", "T test df", "T test p-value", 
                      "Kolmogorov-Smirnov test D", "Kolmogorov-Smirnov test p-value", 
                      "Wilcoxon test W", "Wilcoxon test p-value", 
                      "Kruskal-Wallis test chi-squared", "Kruskal-Wallis test df", 
                      "Kruskal-Wallis test p-value", "Spearman correlation S", 
                      "Spearman correlation rho", "Spearman correlation p-value")
  test.stats.value = c(rmse.can.ref, round(test.t.can.ref$statistic,digits=2), 
                       round(test.t.can.ref$parameter,digits=2), 
                       round(test.t.can.ref$p.value,digits=2), 
                       round(test.ks.can.ref$statistic,digits=2), 
                       round(test.ks.can.ref$p.value,digits=2), 
                       round(test.wil.can.ref$statistic,digits=2), 
                       round(test.wil.can.ref$p.value,digits=2), 
                       round(test.kru.can.ref$statistic,digits=2), 
                       round(test.kru.can.ref$parameter,digits=2), 
                       round(test.kru.can.ref$p.value,digits=2), 
                       round(test.cor.can.ref$statistic,digits=2), 
                       round(test.cor.can.ref$estimate,digits=2), 
                       round(test.cor.can.ref$p.value,digits=2))
  test.stats = cbind(test.stats.name, test.stats.value)
  write.table(test.stats, file="4_statistics_between_daily_series.csv", sep=";", 
              row.names=FALSE)
  
  # Define 15% and 5% error on daily data from ref.serie
  err_pos_15 = ((ref.serie[,4] * 15) / 100) + ref.serie[,4]
  err_neg_15 = (-1 * (ref.serie[,4] * 15) / 100) + ref.serie[,4]
  err_pos_5 = ((ref.serie[,4] * 5) / 100) + ref.serie[,4]
  err_neg_5 = (-1 * (ref.serie[,4] * 5) / 100) + ref.serie[,4]
  
  # Plot of daily data
  png("5_total_precipit_plot_with_15_5.png")
  plot(ref.serie[,4], ref.serie[,4], main="Total precipitation (with p >= 1mm)", 
       ylab="[mm]", xlab="[mm]")
  lines(ref.serie[,4], err_pos_15)
  lines(ref.serie[,4], err_neg_15)
  lines(ref.serie[,4], err_pos_5, col="green")
  lines(ref.serie[,4], err_neg_5, col="green")
  lines(ref.serie[,4], can.serie[,4], type="p", col="red")
  legend('topleft', lty=1, col=c('red','black','green','black'), 
         legend=c("(circles) Candidate values","(circles) Reference values",
                  "(lines) 5% error","(lines) 15% error"), bty='n', cex=.75)
  dev.off()
  
  # Percentage relative error
  rat_month = ((monthly_can - monthly_ref) / monthly_ref) * 100
  rat_month[rat_month==Inf] = NA
  rat_month[rat_month==-Inf] = NA
  rat_month = round(rat_month, digits=2)
  
  # Check if the percentage relative error is good (every value <= treshold)
  M = abs(rat_month)
  M[M<=treshold] = 0
  M[is.na(M)] = 0
  M[M>treshold] = 1
  if (sum(M)>0) good=FALSE else good=TRUE
  
  # Save percentage relative error csv
  Y.rat_month = cbind(years, rat_month)
  colnames(Y.rat_month) = c("Year", months)
  write.table(Y.rat_month, file="6_percentage_relative_error.csv", sep=";", 
              row.names=FALSE) 
  
  # Save percentage relative error boxplot
  png("7_percentage_relative_error_boxplot.png")
  boxplot(rat_month, ylab="Percentage relative error (%)", xlab="Month", 
          main="Box plot of total percentage relative error")
  dev.off()
  
  # If percentage relative error is not good, re-compute previous outputs
  if (!good) {
    new_rat_month = rat_month
    new_rat_month[M==1] = NA
    Y.new_rat_month = cbind(years, new_rat_month)
    colnames(Y.new_rat_month) = c("Year", months)
    write.table(Y.new_rat_month, 
                file="6_filtered_percentage_relative_error.csv", sep=";", 
                row.names=FALSE) 
    
    png("7_filtered_percentage_relative_error_boxplot.png")
    boxplot(new_rat_month, ylab="Percentage relative error (%)", xlab="Month", 
            main="Box plot of total percentage relative error (FILTERED)")
    dev.off()
  }
  
  #### TREND ANALYSIS ####
  
  # Monthly trend analysis
  rat_month.mat = as.matrix(rat_month)
  t.rat_month = as.vector(t(rat_month.mat))
  rat_month.zoo = zooreg(t.rat_month, start=y.start, frequency=12)
  
  # Trend matrix
  trend = zyp.trend.vector(rat_month.zoo, method="yuepilon") 
  trend = as.vector(trend)
  trend.name = c("Trend", "Intercept", "Trend over total period", 
                 "Kendall tau", "Kendall p-value")
  trend = c(trend[2], trend[11], trend[3], trend[5:6])
  trend = round(trend, digits=2)
  tab_trend = cbind(trend.name, trend)
  write.table(tab_trend, file="8_stats_percentage_relative_error_trend.csv", 
              sep=";", row.names=FALSE, col.names=FALSE)
  
  # Plot of monthly trend
  png("8_percentage_relative_error_trend.png")
  x=c(1:length(rat_month.zoo))
  y=as.numeric(rat_month.zoo)
  plot(x, y, ylab="Percentage relative error (%)", 
       xlab=paste("Months (from ",y.start," - ",y.end,")",sep=""), 
       main="Percentage total relative error trend", type="l")
  abline(trend[2], trend[1], col="red")
  dev.off()
  
  # If percentage relative error is not good, re-compute previous outputs
  if (!good) {
    rat_month.mat = as.matrix(new_rat_month)
    t.rat_month = as.vector(t(rat_month.mat))
    rat_month.zoo = zooreg(t.rat_month, start=y.start, frequency=12)
    
    trend = zyp.trend.vector(rat_month.zoo, method="yuepilon") 
    trend = as.vector(trend)
    trend.name = c("Trend", "Intercept", "Trend over total period", 
                   "Kendall tau", "Kendall p-value")
    trend = c(trend[2], trend[11], trend[3], trend[5:6])
    trend = round(trend, digits=2)
    tab_trend = cbind(trend.name, trend)
    write.table(tab_trend, file="8_filtered_stats_percentage_relative_error_trend.csv", 
                sep=";", row.names=FALSE, col.names=FALSE)
    
    png("8_filtered_percentage_relative_error_trend.png")                                 
    x=c(1:length(rat_month.zoo))
    y=as.numeric(rat_month.zoo)
    plot(x, y, ylab="Percentage relative error (%)", 
         xlab=paste("Months (from ",y.start," - ",y.end,")",sep=""), 
         main="Percentage total relative error trend", type="l")
    abline(trend[2], trend[1], col="red")
    dev.off()
  }
  
  #### CLASSES CALCULATION ####
  
  # Quantile evaluated on the ref.serie for thresholds 
  qq_ref.serie = quantile(na.omit(ref.serie[,4]), probs=seq(0,1,0.01))
  
  # Definition of 5 classes 
  q_R_cc = c(1, qq_ref.serie[51], qq_ref.serie[81], qq_ref.serie[96], 
             qq_ref.serie[100])
  
  # Candidate Classes definition
  # Weak
  weak_C = can.serie[can.serie[,4]>=q_R_cc[1]&can.serie[,4]<q_R_cc[2],]
  weak_C = na.omit(weak_C)
  # Mean
  mean_C = can.serie[can.serie[,4]>=q_R_cc[2]&can.serie[,4]<q_R_cc[3],]
  mean_C = na.omit(mean_C)
  # Heavy
  heavy_C = can.serie[can.serie[,4]>=q_R_cc[3]&can.serie[,4]<=q_R_cc[4],]
  heavy_C = na.omit(heavy_C)
  # R95p
  R95_C = can.serie[can.serie[,4]>q_R_cc[4],]
  R95_C = na.omit(R95_C)
  # R99p
  R99_C = can.serie[can.serie[,4]>q_R_cc[5],]
  R99_C = na.omit(R99_C)
  
  # Reference Classes definition
  # Weak
  weak_R = ref.serie[ref.serie[,4]>=q_R_cc[1]&ref.serie[,4]<q_R_cc[2],]
  weak_R = na.omit(weak_R)
  # Mean
  mean_R = ref.serie[ref.serie[,4]>=q_R_cc[2]&ref.serie[,4]<q_R_cc[3],]
  mean_R = na.omit(mean_R)
  # Heavy
  heavy_R = ref.serie[ref.serie[,4]>=q_R_cc[3]&ref.serie[,4]<=q_R_cc[4],]
  heavy_R = na.omit(heavy_R)
  # R95p
  R95_R = ref.serie[ref.serie[,4]>q_R_cc[4],]
  R95_R = na.omit(R95_R)
  # R99p
  R99_R = ref.serie[ref.serie[,4]>q_R_cc[5],]
  R99_R = na.omit(R99_R)
  
  # Vector events and sums
  vC = c(length(weak_C[,4]), length(mean_C[,4]), length(heavy_C[,4]), 
         length(R95_C[,4]), length(R99_C[,4]))
  vR = c(length(weak_R[,4]), length(mean_R[,4]), length(heavy_R[,4]), 
         length(R95_R[,4]), length(R99_R[,4]))
  sum_vC = c(sum(weak_C[,4]), sum(mean_C[,4]), sum(heavy_C[,4]), sum(R95_C[,4]), 
             sum(R99_C[,4]))
  sum_vR = c(sum(weak_R[,4]), sum(mean_R[,4]), sum(heavy_R[,4]), sum(R95_R[,4]), 
             sum(R99_R[,4]))
  cl_q = cbind(vC, sum_vC, vR, sum_vR) 
  
  #### CLASS FORMATTING ####
  
  df_R99_C = class_formatting(R99_C)
  df_R95_C = class_formatting(R95_C)
  df_heavy_C = class_formatting(heavy_C)
  df_mean_C = class_formatting(mean_C)
  df_weak_C = class_formatting(weak_C)
  
  df_R99_R = class_formatting(R99_R)
  df_R95_R = class_formatting(R95_R)
  df_heavy_R = class_formatting(heavy_R)
  df_mean_R = class_formatting(mean_R)
  df_weak_R = class_formatting(weak_R)
  
  #### CLASS ANALYSIS ####
  
  R99_results = class_analysis(df_R99_C, df_R99_R, "R99")
  R95_results = class_analysis(df_R95_C, df_R95_R, "R95")
  heavy_results = class_analysis(df_heavy_C, df_heavy_R, "heavy")
  mean_results = class_analysis(df_mean_C, df_mean_R, "mean")
  weak_results = class_analysis(df_weak_C, df_weak_R, "weak")
  
  #### SUMMARY TABLES ####
  
  # Common events between classes
  l_w = nrow(weak_results)
  l_m = nrow(mean_results)
  l_h = nrow(heavy_results)
  l_R95 = nrow(R95_results)
  l_R99 = nrow(R99_results)
  l_cl = c(l_w, l_m, l_h, l_R95, l_R99)
  
  s_w = apply(weak_results[,2:3], 2, sum)
  s_m = apply(mean_results[,2:3], 2, sum)
  s_h = apply(heavy_results[,2:3], 2, sum)
  s_R95 = apply(R95_results[,2:3], 2, sum)
  s_R99 = apply(R99_results[,2:3], 2, sum)
  
  s_cl = c(s_w, s_m, s_h, s_R95, s_R99)
  s_cl = matrix(s_cl, nrow=5, ncol=2, byrow=T)
  m_cl = cbind(s_cl, l_cl)
  
  # Classes errors
  er_p_max = c(weak_results[1,4], mean_results[1,4], heavy_results[1,4], 
               R95_results[1,4], R99_results[1,4])
  er_p_min = c(weak_results[1,5], mean_results[1,5], heavy_results[1,5], 
               R95_results[1,5], R99_results[1,5])
  media = c(weak_results[1,6], mean_results[1,6], heavy_results[1,6], 
            R95_results[1,6], R99_results[1,6])
  num_15 = c(weak_results[1,7], mean_results[1,7], heavy_results[1,7], 
             R95_results[1,7], R99_results[1,7])
  num_5 = c(weak_results[1,8], mean_results[1,8], heavy_results[1,8], 
            R95_results[1,8], R99_results[1,8])
  errr = cbind(er_p_max, er_p_min, media, num_15, num_5)
  er_mat = matrix(errr, nrow=5, ncol=5)
  mcl = cbind(m_cl, er_mat, cl_q)
  
  # RMSE on common elements
  rmse99 = rmse(R99_results[,2], R99_results[,3])
  rmse95 = rmse(R95_results[,2], R95_results[,3])
  rmseh = rmse(heavy_results[,2], heavy_results[,3])
  rmsem = rmse(mean_results[,2], mean_results[,3])
  rmsew = rmse(weak_results[,2], weak_results[,3])
  rmse_mat = matrix(c(rmsew, rmsem, rmseh, rmse95, rmse99), nrow=5, ncol=1)
  rmse_mat = round(rmse_mat, digits=2)
  
  # Save classes and RMSE
  col.clas = c("weak", "mean", "heavy", "R95", "R99")
  last = cbind(col.clas, mcl, rmse_mat)
  last = cbind(last[,1], last[,7], last[,5], last[,6], last[,12], last[,10], 
               last[,13], last[,11], last[,4], last[,2], last[,3], last[,8], 
               last[,9], last[,14])
  last[,2:14]=round(as.numeric(last[,2:14]),digits=2)
  colnames(last) = c("Class", "Mean daily (mm)", "Err max mean daily (mm)", 
                     "Err min mean daily (mm)", "Num tot events ref", 
                     "Num tot events can", "Sum tot ref (mm)", "Sum tot can (mm)", 
                     "Num common events", "Sum common events ref (mm)", 
                     "Sum common events can (mm)", "Num common outside 15", 
                     "Num common outside 5", "Tot RMSE")
  write.table(last, file="10_class_event_and_RMSE.csv", sep=";", row.names=FALSE)
}

# Print runtime stats and clean all the environment
if (quit) {
  cat("Execution aborted by the user; in", getwd(), 
      "you can find additional information \n")
} else {
  cat("Elaboration completed. You can find your results in", getwd(), "\n")
}
end.time = Sys.time()
cat(capture.output(end.time - start.time))
rm(list=ls())
