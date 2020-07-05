# FitFunctionANDkCluster
Data is comprised of the spike activity of neurons recording in an awake behaving mouse. Mice were presented with a light stimulus that repeated at a regular interval. Averaging the neural responses of each cell over all intervals reveals a gaussian-like function for many cells, with each cell showing an activity peak at shifting timepoints over the interval. The goal of this project is to use gaussian fit and initial response parameters for each cell and determine whether classes of response types exist to represent the passage of time in the brain between the light flashes. 

To Run:
Load one dataset and then type 'curvefitANDcluster' (no quotes) in command window to run.
Two datasets are provided, one with neurons recorded during a 300 ms interval and another with neurons recorded during an 800 ms interval. The data suggest one or two additional classes of response may support longer intervals.
Changing the appropriate flag in the variable initiation of the code enables the ability to plot the activity of each cell and it's gaussian fit.

Skills: 
Demonstrates implentation of fitting functions to data and use of coefficients as additional variables. Use of functions in coding is implemented here, just to demonstrate my ability to implement them if desired, in addition to user queries. Demonstrates use of unsupervised learning algorithm (Kmeans Clustering), and data transformation and handling. Outliers are handled for free parameters using 'IQR outer fence' strategy, and variables for clustering rescaled as min-max (0-1) for each variable. Clustering algorithm was replicated with new random seeds x100 (as a built-in function) with lowest cluster-distance used to choose best result.

NOTE: Code to pull the data from the original files is not included (see SimpleSVM for example of that). Pre-processed structured sample data is included with this package, in addition to some key variables calculated when originally pulling the data. 

General Information:
The passage of interval time is defined as discrete 10 ms time bins.

400 interval samples of spike data were taken and averaged for each neuron. The averages of 2 intervals are used here for function fitting. The averages were smoothed with a 5-bin moving mean filter and min-max normalized (0-1) for these purposes. An additional smoothing option is included in the variable initiation flags at the beginning of the code, after which the data are re-normalized if invoked.

Two figures are included with the sample data. The first is a pseudocolored plot showing the averaged activity of each cell for a 2-interval average, sorted by time to peak activity and permitting visualization of the activity peak shifting for each cell over time during the repeating interval. The second plot shows the average of groups of neurons with activity peaks during the same 50 ms window, also allowing for visualization of the shifting activity over time in addition to notable changes in the shape of the function as cells represent later moments in time during the interval. This observation was the inspiration for defining classes of cell responses using unsupervised clustering. 

Results and Interpretation:
If plotted, it is easy to see that most cells show a good gaussian fit when the interval is shorter (300 ms), and supported by graphing the root of the mean squared error (RMSE). The gaussian functions become slightly wider over the interval, notable by plotting the gaussian fit coefficient 'sigma' and invoking a linear fit. The sigma variable, along with the peak of activity (gaussian coefficiant mu) and a "difference" measure taken from each neuron's initial response to stimulus presentation, were used as classification variables for Kmeans Clustering. Iterative clustering by increasing k by +1 was used, and the sum of squared distances from each neuron's variables to it's assigned cluster was used with the 'elbow method' to determine the number of classes of response types. When compared to the results from data in which the longer interval was presented (800 ms), an additional class of temporal responses are observed when the interval is lengthened. This class may be represented by a non-gaussian temporal function with the longer interval, resulting in poorer gaussian fits (higher RMSE) for some of these cells.
 
 
 
 
