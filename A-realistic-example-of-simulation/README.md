# Realistic-example-of-simulation

The Python file realistic-example contains programs that allow to define all the characteristics (displacement, intensities and kernels) of the example presented in the article where two types of particles move in a cell following three possible types of motion: a Brownian motion, a Brownian motion with linear drift and a confined motion.  
At the end of this program, you find the command to generate the process and then the commands necessary to generate the different plots presented in the article, namely the trajectories of the process, the boxplots of the number of particles at each moment, the histograms of the length of the trajectories... 

You have to put this program in the same folder as the files in the Code folder to run it. You will also have to put in this same folder the .pickle files of this folder so that the python programs can import them.
When you start this program, the current time is displayed. The program ends when the current time reaches the final time of the simulation (here 167.86).