# Real-dataset

In this folder, you find: 

- "M10_Rab11mCherry_09_GaussianWeights.mp4" and "M10_LangerinYFp_09_GaussianWeights.mp4" are the raw sequence for each type of proteins.

- "Utrack_Langerin_GaussianWeights.mp4" and "Utrack_Rab11_GaussianWeights.mp4" are the result of proteins tracking (by the U-track algorithm) for the previous sequences.

- "real-data-video.mp4" is also the result of the traceability of proteins but adding the result of the type of movement followed by each protein: we find in blue the Brownian particles, in red the superdiffusive ones, in green the subdiffusive ones and the Langerin ones are represented by a point while the Rab11 are represented by a triangle.

- The two .csv files containing the real data on which we rely, one for the Langerin particles and one for the Rab11. Each of the two files contains the spatial coordinates of the particles, at each frame, as well as the number of the trajectory to which they belong, as well as the type of motion they have (1 for Brownian, 2 for superdiffusive and 3 for subdiffusive).

- the two files "trajectories_Langerin.png" and "trajectories_Rab11.png" are plots of the Langerin and Rab11 proteins trajectories (with the same color code as above) then the other .png files are plots of the Langerin proteins descriptors: the number of particles of each type per frame, the histogram of the trajectory length, etc...
