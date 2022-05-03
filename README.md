# ArorualLK

Code of "Multi-resolution Data Assimilation for Auroral Energy Flux and Mean Energy"

1. Data Input from different sources: Satellite, Ground-base and Empirical Data.
3. Pre-processing 
   2.1 Interpolated data from satellite observations
   2.2 Downsampling interpolated and empirical data
   2.3 Combined pre-processing data sets
4. Apply Lattice Kriging to obtain the Intermediate results
5. Implement KNN based on combined pre-processing data to generate the weight coefficients
6. Simulate and Predict the Auroral Energy Flux and Mean Energy (two ways: empirical model as Input and Covariate in LatticeKrig function)

