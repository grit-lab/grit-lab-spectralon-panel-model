"""
Working BRF model as of 2023 January 10.
"""

import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import RegularGridInterpolator
import sys

"""
These are the optimized parameters for the Spectralon BRF model.
"""

params = np.array([
0.5696875921173422 , #0
0.27903122059147045 , #1
-1.7476750611608318 , #2
2.5855843488200034 , #3
8.263416395698357 , #4
1.6450590437852088 , #5
1.851701611373643 , #6
0.09100373143807905 , #7
0.5718024986324589 , #8
0.009912786942804066 , #9
0.04026383165673021 , #10
1.1856631736844976 , #11
1.4094448692363608 , #12
3.753748299350997 , #13
0.26592501709914296 , #14
0.25982186894746806 , #15
4.719404248241861 , #16
0.8025276470171108 , #17
0.29606885318487797 , #18
0.4965306551513434 , #19
0.0 , #20
5.866204408080527 , #21
1.212578167338859 , #22
0.4107185337489295 , #23
0.02091690159556928 , #24
0.9782260653794538 , #25
2.944780616597054 , #26
])


"""
These functions define the BRF components.
"""

## define slope-intercept line formula for input angle space
def slope_int(angle_space, slope, intercept):
    return(
    (slope * angle_space) + intercept
    )


## model forward subsurface scattering
def forward_ss(params, azim_space, angle_space):
    
    forward_ss_slope_height, forward_ss_intercept_height = params[0], params[1]
    forward_ss_slope_width, forward_ss_intercept_width = params[2], params[3]
    forward_ss_height_power, forward_ss_width_power = params[4], params[5]
    
    shift = np.pi
    
    small_mat = np.zeros_like(angle_space)
    small_mat.fill(1e-8)
    
    forward_ss_height = np.maximum(small_mat, slope_int(angle_space, forward_ss_slope_height, forward_ss_intercept_height))
    forward_ss_width = np.maximum(small_mat, slope_int(angle_space, forward_ss_slope_width, forward_ss_intercept_width))

    forward_ss_height = forward_ss_height**forward_ss_height_power
    forward_ss_width = forward_ss_width**forward_ss_width_power
    
    return(
    forward_ss_height * (1 + (azim_space-shift)**2 / forward_ss_width)**(-(forward_ss_width+1) / 2)
    )

def diffuse_power(params, azim_space, zen_space, light_angles):
    
    power = params[6]
    coeff = params[7]
    
    return(
           1 - (coeff*zen_space**power)
    )

def diffuse_forward(params, azim_space, angle_space):
    
    slope_diffuse_forward_height, intercept_diffuse_forward_height = params[8], params[9]
    slope_diffuse_forward_width, intercept_diffuse_forward_width = params[10], params[11]
    diffuse_forward_height_power, diffuse_forward_width_power = params[12], params[13]
    
    shift = np.pi
    
    small_mat = np.zeros_like(angle_space)
    small_mat.fill(1e-8)
    
    diffuse_forward_height = np.maximum(small_mat, slope_int(angle_space, slope_diffuse_forward_height, intercept_diffuse_forward_height))
    diffuse_forward_height = diffuse_forward_height**diffuse_forward_height_power
    
    diffuse_forward_width = np.maximum(small_mat,slope_int(angle_space, slope_diffuse_forward_width, intercept_diffuse_forward_width))
    diffuse_forward_width = diffuse_forward_width**diffuse_forward_width_power
    
    return(
        diffuse_forward_height * np.exp(-(azim_space - shift)**2 / diffuse_forward_width**2)
    )


def back_scattering(params, azim_shift_space, zen_space, light_angles):
    
    back_height_slope, back_height_int, back_height_power = params[14], params[15], params[16]
    back_width_azim, back_width_zen = params[17], params[18]
    
    # figure we need a light angle-dependent backscattering height
    small_mat = np.zeros_like(light_angles)
    small_mat.fill(1e-8)
    
    
    # truncate backscattering peak since we don't have that data.
    # to truncate, take the np.min axis=2 of the backscattering peak, or a constant function that is
    # defined on the angle limit of the backscattering peak.
    
    # +/- 10deg in azimuth, +/- 5deg in zenith FROM incidence angle.
    # this is a 20deg window in azimuth and a 10deg window in zenith.
    azim_const = np.zeros_like(azim_shift_space)
    azim_const.fill(np.deg2rad(10))
    
    zen_const = np.zeros_like(zen_space)
    zen_const.fill(np.deg2rad(5))
    
    # roll the shifted azim space to align the indices with the original azim space
    # changed zenith dependence from gaussian to lorentzian
    back_azim = np.exp(-(azim_shift_space)**2 / back_width_azim**2)
    back_zen = np.exp(-(zen_space - light_angles)**2 / back_width_zen**2)
#     back_zen = 1 / (back_width_zen**2 * (zen_space - light_angles)**2 + 1)

    back_azim_const = np.exp(-(azim_const)**2 / back_width_azim**2)
    back_zen_const = np.exp(-(zen_const)**2 / back_width_zen**2)
#     back_zen_const = 1 / (back_width_zen**2 * (zen_const)**2 + 1)
    
    back_height = slope_int(light_angles, back_height_slope, back_height_int)**back_height_power
    
    backscattering_trunc = np.minimum(back_azim*back_zen, back_azim_const*back_zen_const) 
    
    return(
        back_height * backscattering_trunc
    )

def specular(params, azim_space, zen_space, light_angles):
    
    specular_height_slope, specular_height_int, specular_height_power = params[19], params[20], params[21]
    specular_width_azim, specular_width_zen = params[22], params[23]
    
    shift = np.pi
    
    small_mat = np.zeros_like(light_angles)
    small_mat.fill(1e-8)
    
    specular_height = np.maximum(small_mat, slope_int(light_angles, specular_height_slope, specular_height_int))
    specular_height = specular_height**specular_height_power
    
    specular_azim = np.exp(-(azim_space - shift)**2 / specular_width_azim**2)
    specular_zen = np.exp(-(zen_space - light_angles)**2 / specular_width_zen**2)

    return(
        specular_height * specular_azim * specular_zen
    )

def reddening(params, wav):
    redden_slope, redden_int, redden_power = params[24], params[25], params[26]
    return(
        slope_int(wav, redden_slope, redden_int)**redden_power
    )

## define the scattering model
def BRF_model(params, azim_space, azim_shift_space, zen_space, light_angles, wav):
    
    # parse out first column, LUT currently saved 1st column as view zenith angles -- 2023 Jan 11.
    BRF_norm_data = np.loadtxt('BRF_normalization_factor.csv', delimiter=',')[:,1:]
    
    # assume top to down in column is 0-70deg view zenith in 1deg steps.
    # assume left to right in row in 350nm - 2500nm wavelength in 1nm steps.
    BRF_norm = RegularGridInterpolator(
        (np.deg2rad(np.arange(0, 70+1, 1)), np.arange(350, 2500+1, 1)),
        BRF_norm_data)([light_angles, wav])[0]
    
    diffuse_forward_wrt_zen = diffuse_forward(params, azim_space, zen_space)
    diffuse_forward_wrt_light = diffuse_forward(params, azim_space, light_angles)
    diffuse_forward_model = diffuse_forward_wrt_zen*diffuse_forward_wrt_light
    
    diffuse_power_model = diffuse_power(params, azim_space, zen_space, light_angles)
    
    forward_ss_wrt_light = forward_ss(params, azim_space, light_angles)
    forward_ss_wrt_zen = forward_ss(params, azim_space, zen_space)
    forward_ss_model = forward_ss_wrt_light*forward_ss_wrt_zen
    
    back_scattering_model = back_scattering(params, azim_shift_space, zen_space, light_angles)
    
    specular_model = specular(params, azim_space, zen_space, light_angles)
    
    # convert wav to um from nm
    reddening_model = reddening(params, wav*1e-3)
    
    return(
    1/BRF_norm * \
        (diffuse_power_model + \
         diffuse_forward_model*reddening_model + \
         forward_ss_model*reddening_model + \
         specular_model*reddening_model + \
         back_scattering_model)
    )
        
        
        
"""
Inputs into BRF_model():
    
    
azim_space = the relative azimuthal angle between sensor and sun azimuth angles defined between 0 and +2pi,
azim_shift_space = the same azimuthal angle as above but defined between -pi and +pi,
zen_space = the sensor view zenith angle,
light_angles = the solar zenith angle,
wav = wavelength of light in nanometers.

"""

"""
Inputs into BRF_model():
    
    
azim_space = the relative azimuthal angle between sensor and sun azimuth angles defined between 0 and +2pi,
azim_shift_space = the same azimuthal angle as above but defined between -pi and +pi,
zen_space = the sensor view zenith angle,
light_angles = the solar zenith angle,
wav = wavelength of light in nanometers.

"""

print(BRF_model(params, np.deg2rad([120]), np.deg2rad(120), np.deg2rad(48), np.deg2rad(30), 355))
