# Real-dataset

In this folder, you find:

- "M10_Rab11mCherry_09_GaussianWeights.mp4" and "M10_LangerinYFp_09_GaussianWeights.mp4" are the raw sequence for each type of proteins.

- "Utrack_Langerin_GaussianWeights.mp4" and "Utrack_Rab11_GaussianWeights.mp4" are the result of proteins tracking (by the U-track algorithm) for the previous sequences.

- "real-data-video.mp4" gathers the two previous sequences by representing  Langerin vesicles as circles and Rab11 vesicles as triangles. Moreover, each particle is coloured with respect to its motion regime: Brownian in blue, superdiffusive in red and subdiffusive in green. 

- The two .csv files are the results of the tracking procedure for each type of protein. Each file contains the spatial coordinates of the particles, at each frame, the index of the trajectory to which they belong and their motion regime (1: Brownian, 2: superdiffusive, 3: subdiffusive).

- The two images "trajectories_Langerin.png" and "trajectories_Rab11.png" depicts all trajectories of Langerin and  Rab11 proteins, tracked over the image sequences, with the same color label as above (blue: Brownian, red: superdiffusive, green: subdiffusive).  

- “number_per-frame.png” contains the boxplots of the number of particles observed in each frame, with respect to their motion regime (with the same colour label as above). 

- “size_of_trajectories.png” depicts the size of all trajectories (in frames), according to their motion regime (with the same colour label as above). 
