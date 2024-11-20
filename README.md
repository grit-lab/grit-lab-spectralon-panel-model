We provide the Spectralon BRF model as described in the corresponding IEEE TGRS publication (DOI: 10.1109/TGRS.2024.3361392): A Comprehensive BRF Model for Spectralon® and Application to Hyperspectral Field Imagery.

Currently, there are four main files for the model:

1.) BRF_normalization_factor.csv

2.) BRF-model.py

3.) Spectralon_Num4.txt

File 1 contains the values for the normalization coefficient, A(\theta_i, \lambda), which is described in our article. The first column contains the values for the incident zenith angle in degrees. For each incident zenith angle, its corresponding row contains the values for A(\theta_i, \lambda) for every wavelength between 350-2500 nm for a total of 2151 entries.

File 2 contains the Python script that computes the Spectralon BRF for the desired incident and viewing angles, and the wavelength. This script requires Numpy and Scipy. In its current state, the BRF_model() function is not vectorized; it can only take single values for the parameters due to the way it references files 1 and 3. Descriptions of the required inputs can be found within the script.

File 3 contains the 8deg/h calibration coefficients provided with the Spectralon panel from Labsphere, which is the panel that was used for our publication. The first column contains the wavelengths in nm, the second column contains the calibration coefficients, and the third column contains the corresponding uncertainties.



We also provide the diffuse-directional reflectance factor (DDRF) of our panel in the file:

Spectralon-panel-4-diffuse-reflectance_2023Jan18.txt

which was acquired using our laboratory integrating sphere specifically designed for diffuse reflectance measurements.



An example output of the panel model for a viewing zenith of 0 deg and illumination zenith of 10 deg is given in the file:

BRF_i10_e0.txt