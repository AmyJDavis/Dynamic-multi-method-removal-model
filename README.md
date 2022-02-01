# Guide to use a dynamic multi-method removal model 

This document provides guidance for the implementation of the dynamic multi-method removal model described in “An efficient method of evaluating multiple concurrent management actions on invasive populations” by Amy J. Davis, Randy Farrar, Brad Jump, Parker Hall, Travis Guerrant, Kim M. Pepin.  Three R scrips are provided. 

## Files included:
1.	RemDataProcessing.R – Code to process removal data and to run the MCMC model
2.	remfuncs.R – Code to determine the spatial area impacted by each removal event. This code uses latitude and longitude data if available, but will calculate area impacted based on effort generically if no location data are provided
3.	RemMultNomMultMethLambdaMCMC.R – The dynamic multi-method removal model implemented in a custom coded hierarchical Bayesian MCMC code

## Data needs:
The dynamic multi-method removal model is designed to estimate population abundance, population removal rate, and population growth rate using only management removal data. Since management removal often comes from a variety of removal methods this approach allows for and uses all removal methods and estimates removal rates for each method. The data needed for this analysis are often routinely collected as a part of management removal activities. The data needed are:

1.	Methods = The type of removal activity
2.	Date = The date of the removal activity
3.	The location of the removal activity (column for each Latitude/Longitude) if available
4.	Take = The number of individuals removed
5.	Effort = The effort employed in this removal event
6.	Site = A number designated the site of the removal event

The units for the effort will depend on the method. In our example, we used hours in the helicopter, number of nights trapping, and number of days when ground shooting occurred as our effort for those three removal methods. The Site number is to allow for the joint estimation across more than one area in which abundance is desired.  There may only be one Site, in which case a column of ones in the data is needed.  However, if the study area is very large, it may need to be sectioned into smaller areas to reduce bias. In that case, choose the Sites based on a priori knowledge of the study area or use a spatial clustering algorithm to group the removal locations.  
In addition to the data needs described above, to use the code there are two other a priori user inputs needed.  For each removal effort, there needs to be some a prior information about the maximum area of impact for a single unit of effort.  In the case described in the associated article, we used feral swine (Sus scrofa), and used three removal methods: aerial gunning, trapping, and ground shooting.  Therefore, using prior literature or expert opinion, we set the maximum area impacted by an hour in a helicopter, one trap night, or one day of ground shooting.  For the aerial gunning and ground shooting we used expert opinion to inform the area of impact.  For the trapping we used a buffer around the traps (assuming a circle) with a radius based on work by (McRae et al. 2020). The area information needs to be provided on Line 48 in the RemDataProcessing.R code.  The area of each Site needs to be provided (line 45 of RemDataProcessing.R code). This area should include a buffer accounting for the potential for some removal efforts to be on the edge of the property or region of interest.  We used the same buffer that we used around the trap.  However, the approach can be modified for many other species and the area of impact will depend on the species and removal methods used. 

## Steps to run the code
To use this code with your own data, you may need to make some changes (including the ones listed above in the data needs section).  The following steps take you through how to modify RemDataProcessing.R to work with your data. This code was created under R version 4.1.2. 

1.	Save all three scrips in the same working directory (or modify the location of the remfuncs.R and RemMultNomMultMethLambdaMCMC.R on Lines 14 and 15).
2.	Ensure you have all of the packages loaded on your machine. 
3.	Format your data as described above and load it on Line 26. 
4.	If your date format differs from Line 31, change the input so your data loads correctly. I do not use the times in this script. 
5.	Change Line 37 to match the removal methods in your data. Each different factor will have a separate removal rate estimate.  These should be ordered based on the order in which they would occur if events occurred on the same day. 
6.	Modify the ardf data frame on Line 45 to be the areas of each of your Sites.
7.	Modify the detdf data frame on line 48 to give the methods used and the maximum area impacted by one unit of effort for each method.  Keep the last option as Method = “0” and Area = 0, this will be when no removal occurs
8.	Change GPS to TRUE if latitude/longitude information is available and set it to FALSE if not on Line 81. If you want to see the area of impact buffers when you have lat/long data, set plotareas to TRUE.
9.	Set your own design matrix for the linear relationship with growth rate on Line 85.  Or use the default of allowing growth rate to vary across time. 
10.	Lines 88-98 are starting values, prior values, and tuning values for the model.  These may need to be changed to fit your model better.
    a) For example the number of starting values for 'pstart' needs to match the number of removal methods in your data. These values should represent your best guess at the removal rate for one unit of that method (e.g., removal rate for one hour in a helicopter or one night of trapping). It should not matter if your estimate is far from the posterior mean, but try a few value to make sure your starting value is not influencing the results. 
    b) If you are not familiar with Bayesian models and MCMC tuning I will not provide detailed information here, please see additional resources such as van Ravenzwaaij et al. (2018). Very simply I will say that you want to make sure the posterior distribution shows good mixing (the trace plot looks like a fuzzy caterpillar, see Fig. 1A). If the trace plot is a wandering line (Fig. 1B) you need in increase the tuning parameter. If the trace plot looks like it was drawn with an Etch A Sketch (Fig. 1C) you need to decrease the tuning parameter.
  
![alt text](https://github.com/AmyJDavis/Dynamic-multi-method-removal-model/blob/main/MCMC_Tuning.jpg?raw=true)
Figure 1. Example of trace plots of posterior distributions. A) This is an example of the type of trace plot you would like to see, this shows good mixing. B) Shows poor mixing and in particular the tuning parameter is too small, the distribution is not exploring the range efficiently. C) Shows an example of poor mixing where the value is getting stuck and the tuning parameter needs to be reduced to mix better. 

11.	Line 99 set the number of MCMC iterations you want to use.  Keep in mind this code takes a while to run. You might start with 1,000 just to make sure things run smoothly, then based on how long that took and the computing power of your machine, change this to fit your situation. I would recommend 10,000 to 20,000 at least.  
12.	Lines 102-103 run the MCMC code with your inputs

## Model outputs
This model outputs posterior results for abundance (by Site and Month; nt.save), growth rate (by Site and Month; lambda.save), removal rate (by removal method; p.save), and beta estimates for growth rate (depending on the number of covariates examined; beta.save).  The outputs include the entire chains.  To calculate posterior means and credible intervals a burn in should be used, and thinning can be used if desired.  
The model produces diagnostic plots for each abundance, growth rate, removal rate, and beta estimate.  These include trace plots (using all of the data) and posterior distribution plots assuming a burn-in of 50%. The model also produces a plot of the posterior means of the abundance estimates by Site and by Month. 
Code at the end of the RemDataProcessing.R script calculates the posterior means and 95% credible intervals for abundance, growth rate, removal rate, and beta estimates.  

![alt text](https://github.com/AmyJDavis/Dynamic-multi-method-removal-model/blob/main/ModelFlow.jpg?raw=true)
Figure 2. Diagram showing the model input elements (data, informed parameter, and model priors) and the model outputs.  The notation is described in the associated article. 

## Literature Cited
McRae, J. E., P. E. Schlichting, N. P. Snow, A. J. Davis, K. C. VerCauteren, J. C. Kilgo, D. A. Keiter, J. C. Beasley, and K. M. Pepin. 2020. Factors Affecting Bait Site Visitation: Area of Influence of Baits. Wildlife Society Bulletin 44:362-371.

van Ravenzwaaij, D., P. Cassey, and S. D. Brown. 2018. A simple introduction to Markov Chain Monte–Carlo sampling. Psychonomic Bulletin & Review 25:143-154.10.3758/s13423-016-1015-8

