# Co.Rain
Comparing series of Rain

This program works in three steps: statistical analysis, comparison between the cleaned series (with same missing values and daily values <= 1mm) and comparison between precipitation classes. For major details, see the paper by Acquaotta F., Fratianni S. and Venema V. (2016) - Assessment of parallel precipitation measurements networks in Piedmont, Italy (International Journal of Climatology - DOI 10.1002/joc.4606).

1] Statistical analysis
A statistical analysis is carried out on the two series, the candidate series (column 4 of the input file) and the reference one (column 5). The analysis starts with the two raw series, removes all values with daily rain < 1mm and then puts the same missing values on both series, cleaning them. The results are reported in the file "0_statistics_input_file.txt". The statistical analysis calculates the minimum value, the 1st quantile, the median, the mean, the 3rd quantile, the maximum value, the number of missing values, the number of values and the results (W and p_value) of Shapiro-Wilk test. In the 2 files that end with "_daily_avalibillity.csv" are written the monthly number of days used in this analysis; in the 2 files that end with "_monthly_rain.csv" you can find the monthly amount of rain; at the end, in the 2 files with "_day_month_rain.png" are reported the plots of daily and monthly precipitation series with same missing values and daily values <= 1mm.

2] Comparison between the cleaned series (with same missing values and daily values <= 1mm)
A statistical comparison analysis between the two cleaned series is carried out. We apply the Student test (t), the Kolmogorov-Smirmov test (ks), the Wilcoxon rank sum test (W) and the Kruskal test (kru), in addition to the correlation coefficient (cor) with Spearman's method and the Root Mean Square Error (RMSE). The results of the comparison analysis are reported in the file "4_statistics_between_daily_series.csv" where, for every test, we report the statistical value and the p_value. A scatter plot, "5_total_precipitation_plot_with_15_5.png", with daily values is also created. The maximum and minimum error is calculated for the reference series. The maximum error (15) is equal to ± 15% of daily value while the minimum error (5) is equal to ± 5% of daily value. The monthly values is calculated for the two cleaned series and, after that, the percentage relative error is computed; the results are reported in the file "6_percentage_relative_error.csv", with a box plot called "7_percentage_relative_error_boxplot.png". On the percentage relative series, the non-parametric trends are elaborated using the Theil-Sen approach (TSA). The Mann-Kendall test for the trend is then run on the resulting time series to compute the level of significance. The results of the trend are written in the file "8_stats_percentage_relative_error_trend.csv", where we report the trend, the intercept, the trend over all period, the tau of Mann Kendall test and the p_value. After this, if some values of percentage relative error are outside a specific range (set at ±500 by default), we start a re-computation of the only values inside the range, writing results in a file called "6_filtered_percentage_relative_error.csv", re-plotting data in the file "7_filtered_percentage_relative_error_boxplot.png" and re-calculating the trend; the results of the trend are written in the file "8_filtered_stats_percentage_relative_error_trend.csv" and plotted in "8_filtered_percentage_relative_error_trend.png".

3] Comparison between five classes of precipitation events
In this final step we calculate the quantile of the reference series to identify the thresholds of the class. Every daily rain event is classified as weak, mean, heavy, very heavy (R95) or extreme (R99). For each class we create a scatterplot: "9_events_weak.png", "9_events_mean.png", "9_events_heavy.png", "9_events_R95.png" and "9_events_R99.png". In the scatterplots are reported the maximum error, ± 15% of daily value, and the minimum error, ± 5% of daily value. In the file "10_class_events_and_RMSE.csv" for every class, we report: the thresholds, the mean value of maximum and minimum error, number of events for reference and candidate series, precipitation amount for reference and candidate series, number of events recorded in the same day, precipitation amount for reference and candidate series recorded in the same day, number of events outside the range for maximum and minimum error and RMSE.

When you start the program, (if interactive mode is available) you will be asked for an input text file and a folder where results will be stored. If you are running in batch mode, please modify the code to use the correct input file and create the right output folder.

Please, note that a step of cleaning is made on the input file. The file used as starting point of the program and its statistics are available in output number 0 and 1. If these doesn't suit your needs, please discard the whole computation.

# INPUT
The text file has to be formatted in five TAB-separated columns. The first row of the file has to contain the headers (column names) and the first three columns have to contain the year, the month and the day of the series. Column four is the candidate rain serie and column five is the reference rain serie. Missing values must be marked as NA. It is very important to start the file from the first of January (of any year) and end it at the 31st of December (of any year) with at least 5 years of data. See the attached files called exampleX.txt

# OUTPUT
0_statistics_input_file.csv - Important statistics on input file and cleaned input file
1_cleaned_input_file.txt - The file that is really used as input
2_can_Pidro_daily_availability.csv - Number of days of available data on candidate serie
2_can_Pidro_monthly_rain.csv - Total amount of rain by month on candidate serie
2_ref_Parpa_daily_availability.csv - Number of days of available data on reference serie
2_ref_Parpa_monthly_rain.csv - Total amount of rain by month on reference serie
3_can_Pidro_day_month_rain.png - Plot of daily and monthly data on candidate serie
3_ref_Parpa_day_month_rain.png - Plot of daily and monthly data on reference serie
4_statistics_between_daily_series.csv - RMSE, T, KS, Wilcox, Kruskal and Spearman between the two daily series
5_total_precipit_plot_with_15_5.png - Scatter plot of daily events
6_percentage_relative_error.csv - Percentage of relative errors grouped by month (filtered optional)
7_percentage_relative_error_boxplot.png - Box plot of previous data (filtered optional)
8_percentage_relative_error_trend.png - Trend of previous data (filtered optional)
8_stats_percentage_relative_error_trend.csv - Statistics on trend (filtered optional)
9_events_R99.png - Scatter plot of R99 precipitation
9_events_R95.png - Scatter plot of R95 precipitation
9_events_heavy.png - Scatter plot of heavy precipitation
9_events_mean.png - Scatter plot of mean precipitation
9_events_weak.png - Scatter plot of weak precipitation
10_class_event_and_RMSE.csv - Summary of information on classes

Following R-packages have to be installed: "class", "zoo", "hydroGOF", "xts", "hydroTSM", "zyp", "Kendall"

This code has been written under R-Version 3.2.2; for older or newer versions problems might occur.
