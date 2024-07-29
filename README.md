# Temp Dependent QE Plotting Program
 This program plots the raw data from the University of Toledo's PV Measurements, Inc.'s EQE measurement device for tandem solar cell characterization. Given the wavelength and quantum efficiency (QE) of each subcell at different temperatures, the program will plot each curve and display the sample's Jsc in AM1.5G and AM0, Jsc mismatch, and band gap (Eg) as a function of temperature. 

![QE Example Plot](https://github.com/KhoaWeston/Temp-Dependent-QE-Plotting-Program/blob/master/Screenshot%202024-07-26%20145629.png)



## How to Build
1. Clone the repository
2. Open the project folder using your personal Python IDE
3. Run the following commands in the terminal:
```
pip install matplotlib
```
```
pip install numpy
```

## How to Use
1. Specify the directory where your data is stored in the `directory` variable.
2. Change the `sample_ID` variables to correspond with your personal file naming convention.
3. Run the program.

## Result
- A new folder will be created within the specified directory named `Plots`
- Within the new folder, you will get a couple of plots:
  * QE curves for each temperature stacked onto one plot
  * Jsc and Eg vs temperature with their corresponding temperature coefficients 

## References
