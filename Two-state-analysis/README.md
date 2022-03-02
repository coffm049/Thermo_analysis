# Two-state-analysis
# The analysis will be contained entirely within one .R file. The different files in this repository will therefore
# represent different troubleshooting processes in the development.

Complete analysis of thermodynamic transitions in DSC data with application of confidence regions and log-likeliness contours for model evaluation.

This analysis was made to be compatible with .csv file types with an arbitrary number of scans. The column names for the 
scans are arbitrary, but the first column must be the temperature data followed by the heat capacity data, and any 
subsequent experiment should follow the same ordering. (i.e. columns should be as follows: Temp, Cp, Temp, Cp, etc...)
Note: Data does not need to be baselined or converted to Kelvin. The program will take care of that if you haven't already.
Steps for the analysis.
1. Load in the .csv file with DSC data as df.dsc through the interactive file chooser.
2. Enter the lower and upper temperatures to trim the data. This should include any oddities that are resultant of the machine
at the beginning and end of a scan. This will also serve as the peripheral bounds for the upper and lower baselines, respectively.
3. Choose the remaining bounds for the baselines.
4. The Graphical.features function will automatically analyze the data and create estimates for the three thermeodynamic
parameters of interest (dH, dCp, Tm) as well as a graph annotating the features.
5. The two.state function will run a nonlinear regression (Gauss-Newton method) starting from the initial parameters estimated
by the Graphical.features function. It will also produce graphs displaying the evaluation of the confidence regions as well as
log-likeliness contours.
