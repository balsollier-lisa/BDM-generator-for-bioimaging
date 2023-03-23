# Code

In this folder, you find two Python files:

- ''process.py'' gathers all Python functions necessary to generate a BDM process. The main function, named ''proctotal’’, takes in argument all the characteristics of the process (final time of simulation, move, intensities of births, of deaths and of transformation, transition probabilities for the births, for the deaths and for the transformations) and returns the coordinates of the generated particles at each time of the simulation (over a fine temporal grid). See the folders ''A-simple-example-of-simulation'' and  ''A-realistic-example-of-simulation'' for examples.  

- ''draw.py'' contains several functions that take in argument the output of the function “proctotal” (see above), in order to plot the generated trajectories, or the histograms of their size,  or the boxplots of the number of particles per frame, or to show a movie of the simulated process.  See the folders ''A-simple-example-of-simulation'' and  ''A-realistic-example-of-simulation'' for examples.  
