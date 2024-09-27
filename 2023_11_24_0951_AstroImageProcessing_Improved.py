# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 09:50:41 2023

@author: kyebchoo
"""
#%%

# NUMPY
import numpy as np

# SCIPY
import scipy as sci
# import scipy.signal as sig
import scipy.special as scisp
import scipy.stats as stats
from scipy.optimize import curve_fit

# MATPLOTLIB
import matplotlib.pyplot as plt

# PANDAS

# import os
import math
import pandas as pd
import pyvisa
# import time
from matplotlib import cbook
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
# import serial
# import re
# from pygame import mixer
# from time import sleep
from datetime import datetime
from tqdm import tqdm
# from tqdm.auto import tqdm
# import keyboard
# from datetime import date, datetime
from matplotlib.pyplot import figure

from astropy.io import fits
# from aspylib import astro

import matplotlib.colors as colors

pd.options.mode.chained_assignment = None  # default='warn'

from scipy import ndimage
# import imageio.v3 as iio
# import ipympl
import skimage as ski

from matplotlib.gridspec import GridSpec
from matplotlib.ticker import FormatStrFormatter

# Example figsize control, advise to use x width as 9cm to match with latex styles
# cm = 1/2.54
# figsize=(9*cm, 5*cm)




#%%
'''
PLOT CONTROLS
'''

# Example figsize control, advise to use width as 9cm to match with latex styles
cm = 1/2.54

# https://matplotlib.org/stable/api/matplotlib_configuration_api.html#matplotlib.rcParams

# setting default plot dimensions as 9cm by 9cm
plt.rcParams['figure.figsize'] = [9.0*cm, 9.0*cm]
plt.rcParams['figure.dpi'] = 1200

# set font sizes
SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 12

LATEX_SIZE_HEADER = 9
LATEX_SIZE_STANDARD = 8
LATEX_SIZE_SMALL = 7
LATEX_SIZE_SMALLER = 5

plt.rc('font', size = LATEX_SIZE_SMALL)          # controls default text sizes
plt.rc('axes', titlesize = LATEX_SIZE_SMALL)     # fontsize of the axes title
plt.rc('axes', labelsize = LATEX_SIZE_SMALL)     # fontsize of the x and y labels
plt.rc('xtick', labelsize = LATEX_SIZE_SMALLER)    # fontsize of the tick labels
plt.rc('ytick', labelsize = LATEX_SIZE_SMALLER)    # fontsize of the tick labels
plt.rc('legend', fontsize = LATEX_SIZE_SMALLER)    # legend fontsize
plt.rc('figure', titlesize = LATEX_SIZE_SMALL)   # fontsize of the figure title



# setting global line sizes
plt.rc('lines', linewidth = 0.5)
plt.rc('axes', linewidth = 0.5)
plt.rc('grid', linewidth = 0.5)
plt.rc('patch', linewidth = 0.5)
plt.rc('ytick.major', width = 0.5)
plt.rc('xtick.major', width = 0.5)


# setting global marker sizes
plt.rc('lines', markersize = 5)


# to be used to export and show images throughout
def export_plot(filename, item = plt, figsize = (9.0*cm, 9.0*cm)):
    item.tight_layout()
    plt.show()
    item.savefig(str('%s_%s_%4.2fcm_%4.2fcm.png' % (filename, datetime.now().strftime("%Y%m%d_%H_%M_%S"), figsize[0]/cm, figsize[1]/cm)))
    pass



#%%
def line(x, m, c):
    return m*x + c
    
    
def poisson(x, lamb, shift, A):
    return A*(lamb**(x - shift))*(np.exp(-1 * lamb))/(scisp.factorial((x - shift)))
    
def gaussian(x, A, b, sig):
    out = A * np.exp(-((x-b)**2)/(2*(sig**2)))
    return out

def gaussian_yshift(x, A, b, sig, y_shift):
    out = y_shift + A * np.exp(-((x-b)**2)/(2*(sig**2)))
    return out

def in_circle(x, y, x_center, y_center, radius):
    
    dx = abs(x - x_center)
    dy = abs(y - y_center)
    
    if (dx**2 + dy**2) < (radius)**2:
        return True
    else:
        return False
    pass

def in_annulus(x, y, x_center, y_center, radius_min, radius_max):
    
    dx = abs(x - x_center)
    dy = abs(y - y_center)
    
    if (dx**2 + dy**2) < (radius_max)**2 and (dx**2 + dy**2) > (radius_min)**2:
        return True
    else:
        return False
    pass

#%%


# define global variables:
bin_count = 10000
# catalogue = pd.DataFrame()
# catalogue_initialised = False
# catalogue = pd.DataFrame()
# hdulist = pd.DataFrame()
# hdulistdf = pd.DataFrame()
# hdulistdfflat = pd.DataFrame()
# hdulist_header = pd.DataFrame()
# Dataframe = pd.DataFrame()
# hdulistdf_mask = pd.DataFrame()
# MAGZPT = 0
# MAGZRR = 0
# popt_background = 0
# pcov_background = 0
# hdulistdf2 = pd.DataFrame()
# remove_count = 0
# signal_flags = pd.DataFrame()
# local_background_array = pd.DataFrame()




#%%

def initialise_catalogue():
    global catalogue, catalogue_initialised
    catalogue_initialised = True
    catalogue = pd.DataFrame(columns = ['X-Location', 'Y-Location', 'Brightness', 'Sigma', 'Object Sigma Limit', 'Background Sigma Lowest', 'Background Sigma Highest'])
    pass

def preview_catalogue(whole = False):
    if whole == False:
        return catalogue[['X-Location', 'Y-Location', 'Brightness', 'Sigma']]
    elif whole == True:
        return catalogue.to_string()
    pass

   
    
def load(filename = "Astro\Fits_Data\mosaic.fits"):
    global hdulist, hdulistdf, hdulistdfflat, hdulist_header, Dataframe, hdulistdf_mask
    hdulist = fits.open(filename)
    hdulist_header = hdulist[0].header
    hdulistdf = pd.DataFrame(hdulist[0].data) # Pandas array 2d
    # hdulistdf_mask = hdulistdf                # Duplicate mask
    
    hdulistdfflat = hdulistdf.values.flatten() # numpy array 1d format
    counts, bins = np.histogram(hdulistdfflat, bins = bin_count)
    Dataframe = pd.DataFrame({'Counts': counts, 'Brightness': bins[0:-1]})
    # print(Dataframe)
    get_MAGZPT()
    get_MAGZRR()
    initialise_catalogue()
    initialise_exception_list()
    pass

def load_excel(filename):
    global hdulist, hdulistdf, hdulistdfflat, Dataframe
    hdulist = pd.read_excel(filename)
    hdulistdf = pd.DataFrame(hdulist) # Pandas array 2d
    hdulistdfflat = hdulistdf.values.flatten() # numpy array 1d format
    counts, bins = np.histogram(hdulistdfflat, bins = bin_count)
    Dataframe = pd.DataFrame({'Counts': counts, 'Brightness': bins[0:-1]})
    # print(Dataframe)
    initialise_catalogue()
    initialise_exception_list()
    pass

def HDUList():
    return hdulist

def HDUList_header():
    return hdulist_header

def get_MAGZPT():
    global MAGZPT
    MAGZPT = HDUList_header()[157]
    pass
    # return MAGZPT

def get_MAGZRR():
    global MAGZRR
    MAGZRR = HDUList_header()[158]
    pass
    # return MAGZRR
    

def DataFrame2D():
    return hdulistdf

def DataFrameLinear():
    return Dataframe

def preview_linear():
    plt.figure()
    plt.hist(hdulistdfflat, bins = bin_count)
    plt.ylim(0, 0.3E6)
    export_plot('preview_linear', figsize = (9, 9), item = plt)
    pass

def fit_background(peak = 3400):
    global popt_background, pcov_background 
    popt_background, pcov_background = curve_fit(gaussian, DataFrameLinear()['Brightness'], DataFrameLinear()['Counts'], p0 = [3000, peak, 500])
    # popt_background, pcov_background = curve_fit(poisson, DataFrame()['Brightness'], DataFrame()['Counts'], p0 = [3, peak, 50000])
    
    plt.figure()
    plt.hist(hdulistdfflat, bins = bin_count)
    plt.plot(np.arange(0, 50000), gaussian(np.arange(0, 50000), popt_background[0], popt_background[1], popt_background[2]))
    # plt.plot(np.arange(0, 50000), poisson(np.arange(0, 50000), popt_background[0], popt_background[1], popt_background[2]))
    plt.ylim(0, 0.3E6)
    export_plot('fit_background1', item = plt)
    
    fig, ax = plt.subplots(figsize=(8, 5), dpi=120)
    ax.imshow(hdulistdf)
    ax.invert_yaxis()
    export_plot('fit_background2', item = plt)
    
    return popt_background, pcov_background


def convert_to_dms(deg, display = False):
    decimals, number = math.modf(deg)
    d = int(number)
    m = int(decimals * 60)
    s = (deg - d - m / 60) * 3600.00
    # return '{}º{}\'{:.2f}"{}'.format(abs(d), abs(m), abs(s))
    if display == False:
        return ((d, m, s))
    elif display == True:
        return str('%.0f:%02.0f:%05.2f' % (d, m, s))
    pass


    
    
def convert_to_hms(deg, display = False):
    # deg, arcmin, arcsec = convert_to_deg_arcmin_arcsec(deg, display = False)
    hours = deg * 24 / 360
    decimal, h = math.modf(hours)
    minutes = decimal * 60
    decimal, m = math.modf(minutes)
    s = decimal * 60
    if display == False:
        return ((h, m, s))
    elif display == True:
        return str('%.0f:%02.0f:%05.2f' % (h, m, s))
    pass



def convert_to_deg(input_array, form):
    if form == 'hms':
        hrsdeg = input_array[0] * 360 / 24
        mindeg = input_array[1] * 15 / 60
        secdeg = input_array[2] * 15 / 3600
        return hrsdeg + mindeg + secdeg
    elif form == 'dms':
        degdeg = input_array[0]
        mindeg = input_array[1] / 60
        secdeg = input_array[2] / 3600
        return degdeg + mindeg + secdeg
    else:
        raise Exception('Invalid form input (\'mhs\'\\\'dms\')')
        pass
    pass


def generate_total_window_area(unit = 'sr'):
    RA  = convert_to_deg((0, 0, 2570 * 0.258), form = 'dms')
    DEC = convert_to_deg((0, 0, 4611 * 0.258), form = 'dms')
    # note this is only an approximation, the surface elecment is not cartesian in celestial cooredinates
    if unit == 'sqdeg':
        return RA * DEC
    elif unit == 'sr':
        return RA * DEC *(np.pi/180)**2
    else:
        raise Exception('Invalid unit input.')




def map_pixel(form, initial_offset, pixel, arcsec_per_pixel, display = True):
    initial_offset = convert_to_deg(initial_offset, form = form)
    arcsec = arcsec_per_pixel * pixel
    delta = convert_to_deg([0, 0, arcsec], form = 'dms')
    position = initial_offset + delta
    
    if form == 'dms':
        position = convert_to_dms(position)
        pass
    elif form == 'hms':
        position = convert_to_hms(position)
        pass
    else:
        raise Exception('Invalid form input (\'mhs\'\\\'dms\')')
    
    if display == False:
        return position
    elif display == True and form == 'dms':
        return str('%.0f$^{\circ}$%02.0f$^{\prime}$%05.2f$^{\prime\prime}$' % (position[0], position[1], position[2]))
    elif display == True and form == 'hms':
        return str('%.0f$^{h}$%02.0f$^{m}$%05.2f$^{s}$' % (position[0], position[1], position[2]))
    else:
        raise Exception('Invalid display input.')
        pass
    pass



def preview_2D():
    
    from astropy.visualization import imshow_norm, MinMaxInterval, SqrtStretch
    
    fig, ax = plt.subplots(nrows = 1, ncols = 2, sharey = True)
    
    load()
    # im, norm = imshow_norm(hdulistdf, ax, origin='lower', interval=MinMaxInterval(), stretch=SqrtStretch())
    im, norm = imshow_norm(hdulistdf, ax[0], origin='lower', vmin = 3400, vmax = 36000, stretch=SqrtStretch())
    # im = ax.imshow(hdulistdf, cmap = 'jet', vmin = 3400, vmax = 36000)
    
    load_excel('HDUlist_replaced_with_3416.xlsx')
    im, norm = imshow_norm(hdulistdf, ax[1], origin='lower', vmin = 3400, vmax = 36000, stretch=SqrtStretch())
    
    fig.colorbar(im, label = 'Brightness (ADU)')
    
    
    x_ticks = np.linspace(0, 2570, 3)
    x_ticklabels = []
    for i in x_ticks:
        x_ticklabels.append(map_pixel('hms', [10, 46, 00.00], i, 0.258, display = True))
        pass
    ax[0].xaxis.set(ticks = x_ticks, ticklabels = x_ticklabels) 
    ax[1].xaxis.set(ticks = x_ticks, ticklabels = x_ticklabels) 
    ax[0].tick_params(axis = "x", rotation = 45)
    ax[1].tick_params(axis = "x", rotation = 45)
    
    y_ticks = np.linspace(0, 4610, 5)
    y_ticklabels = []
    for i in y_ticks:
        y_ticklabels.append(map_pixel('dms', [59, 1, 59.99], i, 0.258, display = True))
        pass
    ax[0].yaxis.set(ticks = y_ticks, ticklabels = y_ticklabels)
    ax[0].tick_params(axis = "y", rotation = 45)
    
    ax[1].set_xlim(0, 2570)
    ax[1].set_ylim(0, 4610)
    ax[0].set_xlim(0, 2570)
    ax[0].set_ylim(0, 4610)
    
    ax[0].set_title('Before Removal')
    ax[1].set_title('After Removal')

    fig.supxlabel('Declination (deg:arcmin:arcsec)')
    fig.supylabel('Right Ascension (hrs:min:sec)')
    
    plt.tight_layout()
    # plt.grid()
    export_plot('preview_linear', item = fig, figsize = (9, 9))
    pass

load_removed_loaded = False
def load_removed():
    global hdulistdf2, load_removed_loaded
    hdulist2 = pd.read_excel('HDUlist_replaced_with_3416.xlsx')
    hdulistdf2 = pd.DataFrame(hdulist2)
    load_removed_loaded = True
    pass

#%%
def preview_2D_alt():
    
    from astropy.visualization import imshow_norm, MinMaxInterval, SqrtStretch
    
    set_figsize = (9*cm, 24*cm)

    fig = plt.figure(figsize = set_figsize, layout = 'tight')
    gs = fig.add_gridspec(2, 2,  width_ratios = (10, 1), height_ratios = (1, 1), hspace = 0, wspace = 0)

    colorbar = fig.add_subplot(gs[:, 1])
    upper = fig.add_subplot(gs[0, 0])
    lower = fig.add_subplot(gs[1, 0], sharex = upper)
    
    global hdulistdf
    load()
    print('Grid 1 Loaded: ', datetime.now())
    
    # im, norm = imshow_norm(hdulistdf, ax, origin='lower', interval=MinMaxInterval(), stretch=SqrtStretch())
    im, norm = imshow_norm(hdulistdf, upper, origin='lower', vmin = 3400, vmax = 36000, stretch=SqrtStretch(), cmap = 'jet', aspect = 'auto')
    # im = ax.imshow(hdulistdf, cmap = 'jet', vmin = 3400, vmax = 36000)
    
    global hdulistdf2, load_removed_loaded
    if load_removed_loaded == False:
        load_removed()
        pass
    print('Grid 2 Loaded', datetime.now())
    
    
    # load_excel('HDUlist_replaced_with_3416.xlsx')
    im, norm = imshow_norm(hdulistdf2, lower, origin='lower', vmin = 3400, vmax = 36000, stretch=SqrtStretch(), cmap = 'jet', aspect = 'auto')
    
    cbar = fig.colorbar(im, cax = colorbar)
    cbar.formatter.set_powerlimits((0, 0))
    colorbar.tick_params(axis = "y", direction = 'in', rotation = 0, pad = - 13, width = 0.5)
    # colorbar.yaxis.set_label_position("left")
    colorbar.yaxis.tick_left()
    colorbar.set_ylabel('Luminosity (ADU)', fontsize = LATEX_SIZE_SMALL, rotation = - 90, labelpad = 8)
    
    
    x_ticks = np.linspace(400, 2200, 3)
    x_ticklabels = []
    for i in x_ticks:
        x_ticklabels.append(map_pixel('hms', [10, 46, 00.00], i, 0.258, display = True))
        pass
    lower.xaxis.set(ticks = x_ticks, ticklabels = x_ticklabels)
    lower.tick_params(axis = "x", rotation = 0)
    upper.tick_params(axis = "x", labelbottom = False, bottom = False)
    

    y_ticks = np.linspace(600, 4500, 4)
    y_ticklabels = []
    for i in y_ticks:
        y_ticklabels.append(map_pixel('dms', [59, 1, 59.99], i, 0.258, display = True))
        pass
    upper.yaxis.set(ticks = y_ticks, ticklabels = y_ticklabels)
    upper.tick_params(axis = "y", rotation = 90)
    lower.yaxis.set(ticks = y_ticks, ticklabels = y_ticklabels)
    lower.tick_params(axis = "y", rotation = 90)
    
    lower.set_xlim(0, 2570)
    lower.set_ylim(0, 4610)
    upper.set_xlim(0, 2570)
    upper.set_ylim(0, 4610)
    
    # upper.set_title('Before Removal', pad = -0.1)
    # lower.set_title('After Removal', pad = -0.4)

    lower.set_xlabel('RA (J2000)', loc = 'center')
    upper.set_ylabel('DEC (J2000)')
    lower.set_ylabel('DEC (J2000)')
    # fig.supylabel('RA (J2000)')
    

    export_plot('preview_2D_alt', item = fig, figsize = set_figsize)
    pass
#%%


def get_value(x, y):
    return DataFrame2D()[x][y]


def replace(x, y, replace_value):
    DataFrame2D()[x][y] = replace_value
    pass


def remove_pixels(filename = 'Removals.xlsx', replace_value = -1):
    
    removelist = pd.read_excel(filename)
    global remove_count
    remove_count = 0
    
    entries = len(removelist['x left'])
    print()
    for i in tqdm(range(0, entries), position = 0, leave = True):
        for j in range(int(removelist['x left'][i]), int(removelist['x right'][i] + 1)):
            for k in range(int(removelist['y bottom'][i]), int(removelist['y top'][i] + 1)):
                if get_value(j, k) != -1:
                    replace(j, k, replace_value)
                    remove_count += 1
                    pass
                pass
            pass
        pass
    pass


def remove_above_threshold(upper_threshold = 50000, replace_value = 3416):
    for i in tqdm(range(0, len(len(DataFrame2D().columns))), position = 0, leave = True):
        for j in range(0, len(DataFrame2D().index)):
            if get_value(i, j) > upper_threshold:
                replace(i, j, replace_value)
                pass
            pass
        pass
    

def object_search_max(xmin, xmax, ymin, ymax, background_threshold = 3500):
    
    max_count = 0
    max_location = [-1, -1]
    
    for i in tqdm(range(xmin, xmax + 1), position = 0, leave = True):
        for j in range(ymin, ymax + 1):
            if get_value(i, j) > max_count and get_value(i, j) > background_threshold:
                max_count = get_value(i, j)
                max_location = [i, j]
                pass
            pass
        pass
    
    return max_location

def get_x_cut(x_mean, y_mean, width = 15, plot = False):
    
    x_cut = []
    
    for i in range (x_mean - width, x_mean + width):
        x_cut.append(get_value(i, y_mean))
        pass
    
    if plot == True:
        plt.figure()
        plt.scatter(np.arange(0, len(x_cut)), x_cut)
        export_plot('get_x_cut', item = plt)
        pass
    
    return x_cut

def get_y_cut(x_mean, y_mean, width = 15, plot = False):
    
    y_cut = []
    
    for i in range (y_mean - width, y_mean + width):
        y_cut.append(get_value(x_mean, i))
        pass
    
    if plot == True:
        plt.figure()
        plt.scatter(np.arange(0, len(y_cut)), y_cut)
        export_plot('get_y_cut', item = plt)
        pass
    
    return y_cut

def get_object_radius(x_mean, y_mean, sigma, width = 15, plot = False):
    
    x_array = get_x_cut(x_mean, y_mean, width, plot = False)
    y_array = get_y_cut(x_mean, y_mean, width, plot = False)
    
    poptx, pcovx = curve_fit(gaussian_yshift, np.arange(len(x_array)), x_array)
    popty, pcovy = curve_fit(gaussian_yshift, np.arange(len(y_array)), y_array)
    
    one_sigma_x = poptx[2]
    one_sigma_y = popty[2]
    
    if one_sigma_x > one_sigma_y:
        one_sigma = one_sigma_x
        pass
    else:
        one_sigma = one_sigma_y
        pass
    
    if plot == True:
        plt.figure()
        plt.scatter(np.arange(len(x_array)), x_array)
        plt.axvline(poptx[1] - sigma * one_sigma)
        plt.axvline(poptx[1] + sigma * one_sigma)
        export_plot('get_object_radius1', item = plt)
        
        plt.figure()
        plt.scatter(np.arange(len(y_array)), y_array)
        plt.axvline(popty[1] - sigma * one_sigma)
        plt.axvline(popty[1] + sigma * one_sigma)
        export_plot('get_object_radius2', item = plt)
        pass
    
    return abs(sigma * one_sigma)

def get_object_count_sum(x, y, sigma_signal = 3, sigma_background_min = 3, sigma_background_max = 4, replace_value = 3416):

    # background
    background_sum_count = 0
    background_pixel_count = 0
    
    radius_annulus_min = get_object_radius(x, y, sigma_background_min, width = 20, plot = False)
    radius_annulus_max = get_object_radius(x, y, sigma_background_max, width = 20, plot = False)
    
    
    for i in range(x - int(radius_annulus_min) - 1, x + int(radius_annulus_min) + 1 + 1):
        for j in range(y - int(radius_annulus_min) - 1, y + int(radius_annulus_min) + 1 + 1):
            
            condition = in_annulus(i, j, x, y, radius_annulus_min, radius_annulus_max)
            
            if condition == True:
                background_sum_count += get_value(i, j)
                background_pixel_count += 1
                replace(i, j, replace_value)
                pass
            pass
        pass
    
    average_background = background_sum_count/background_pixel_count
    
    # object
    radius = get_object_radius(x, y, sigma_signal, width = 20, plot = False)
    
    sum_count = 0
    pixel_count = 0
    
    for i in range(x - int(radius) + 1, x + int(radius) + 1 + 1):
        for j in range(y - int(radius) + 1, y + int(radius) + 1 + 1):
            
            condition = in_circle(i, j, x, y, radius)
            
            if condition == True:
                sum_count += get_value(i, j)
                pixel_count += 1
                replace(i, j, replace_value)
                pass
            pass
        pass
    
    
    total_background = average_background * pixel_count
    
    object_count = sum_count - total_background
    
    return object_count
    
def catalogue_window(xmin, xmax, ymin, ymax, sigma_signal_input = 3, sigma_background_min_input = 3, sigma_background_max_input = 4, background_threshold = 3500, upper_threshold = 50000, exception_value = 3416):
    
    global catalogue_initialised
    if catalogue_initialised != True:
        raise Exception("No catalogue detected, initialise using : initialise_catalogue()")
    
    run_count = 0
    
    global catalogue
    continue_scanning = True
    
    while continue_scanning == True:
        
        max_coordinate = object_search_max(xmin, xmax, ymin, ymax, background_threshold)
        
        if max_coordinate == [-1, -1]:
            continue_scanning = False
            pass
        elif max_coordinate != [-1, -1]:
            
            max_pixel = get_value(max_coordinate[0], max_coordinate[1])
            
            try:
                print(max_coordinate)
                sigma = get_object_radius(max_coordinate[0], max_coordinate[1], 1, width = 20, plot = False)
                
                if sigma >= 1.5:
                    brightness = get_object_count_sum(max_coordinate[0], max_coordinate[1], sigma_signal_input, sigma_background_min_input, sigma_background_max_input, replace_value = exception_value)
                    
                    if max_pixel <= upper_threshold:                    
                        catalogue.loc[len(catalogue)] = [max_coordinate[0], max_coordinate[1], brightness, sigma, sigma_signal_input, sigma_background_min_input, sigma_background_max_input]
                        run_count += 1
                        print(run_count)
                        print('')
                        pass
                    pass
                
                elif sigma < 1.5:
                    replace(max_coordinate[0], max_coordinate[1], exception_value)
                    pass
                
            except RuntimeError:
                print('RuntimeError: Iteration not converging, trial passed')
                replace(max_coordinate[0], max_coordinate[1], exception_value)
                pass
            except KeyError:
                print('KeyError: Iteration not converging, trial passed')
                replace(max_coordinate[0], max_coordinate[1], exception_value)
                pass
            pass
        pass
    pass



def export():
    hdulistdf.to_excel('HDUlist_edited.xlsx')
    pass

def export_catalogue():
    catalogue.to_excel('Catalogue.xlsx')
    pass

def get_removed_count():
    return remove_count

def get_instrumental_mag(counts):
    return -2.5 * np.log10(counts)

def get_calibrated_mag(counts):
    return get_instrumental_mag(counts) + MAGZPT

def calibrate_image():
    global hdulistdf
    hdulistdf = get_calibrated_mag(hdulistdf)
    pass



def initialise_background_scanning(set_radius_input = 15):
    global first_run_query, set_radius, x_columns_counts, previous_y, background_scanning_initialised
    
    first_run_query = True
    set_radius = set_radius_input
    x_columns_counts = np.zeros((2 * set_radius + 1, 2 * set_radius + 1))
    previous_y = - 1
    background_scanning_initialised = True
    pass

initialise_background_scanning(set_radius_input = 15)



def get_local_pixel_array(x, y):
    radius = set_radius
    
    global x_columns_counts, previous_y, first_run_query
    
    # summing columns to create new array for first run or when y axis changes
    if first_run_query == True or previous_y != y:
        previous_y = y
        first_run_query = False
        for i in range(0, 2 * radius + 1):
            for j in range(0, 2 * radius + 1):
                x_columns_counts[i][j] = get_value(x + i - radius, y + j - radius)
                pass
            pass
        pass
    
    # if there is a previous matrice, remove last row and add in new row
    else:
        x_columns_counts = np.delete(x_columns_counts, 0, 0)
        total = np.zeros((1, 2 * radius + 1))
        for j in range(0, 2 * radius + 1):
            total[0][j] = get_value(x + radius, y + j - radius)
            pass
        x_columns_counts = np.append(x_columns_counts, total, axis = 0)
        pass
    
    return x_columns_counts

def get_local_background(x, y, scan_range = (3350, 3450)):
    
    pixel_matrice = get_local_pixel_array(x = x, y = y)
    
    if sum(i < scan_range[1] for i in pixel_matrice.flatten()) < 0.5 * 4 * (set_radius**2):
        count_max = scan_range[1]
        sigma_max = 0
        pass
    else:
        try:
            counts, bins = np.histogram(pixel_matrice.flatten(), bins = 10, range = scan_range)
    
            index_max = np.argmax(counts)
            count_max = bins[index_max]
        
            if count_max >= scan_range[1]:
                count_max = scan_range[1]
                sigma_max = 0
                pass

            else:
                popt_background, pcov_background = curve_fit(gaussian, bins[0: -1], counts, p0 = [150, 3420, 2])
                # popt_background, pcov_background = curve_fit(gaussian, bins[0: -1], counts)
        
                # A, b, sigma
                # sigma_max = np.sqrt(pcov_background[2, 2])
                sigma_max = np.abs(popt_background[2])
                # print(popt_background[2])
                # print(pcov_background[1, 1])
                pass
            pass
        except RuntimeError:
            count_max = scan_range[1]
            sigma_max = 0
            pass
        except 'OptimizeWarning':
            count_max = scan_range[1]
            sigma_max = 0
            pass
        pass
    
    return count_max, sigma_max



def flag_signals(xmin = 200, xmax = 2400, ymin = 300, ymax = 4300, sigmas = 5, scan_range = (3350, 3450)):
    
    global local_background_array, signal_flags, background_scanning_initialised, xorigin, yorigin
    
    if background_scanning_initialised != True:
        raise Exception("Initialise programme with: initialise_background_scanning()")
        pass
    
    xorigin = xmin
    yorigin = ymin
    
    local_background_array = np.zeros((xmax - xmin + 1, ymax - ymin + 1))
    signal_flags = np.zeros((xmax - xmin + 1, ymax - ymin + 1))
    
    for i in tqdm(range(ymin, ymax + 1), position = 0, leave = True):
        for j in range(xmin, xmax + 1):
            # print(i, j)
            local_background, fluctuation = get_local_background(j, i, scan_range = scan_range)
            reading = get_value(j, i)
            local_background_array[j - xmin][i - ymin] = local_background
            
            if reading >= (local_background + fluctuation * sigmas):
                signal_flags[j - xmin][i - ymin] = 1
                pass
            
            pass
        pass
    background_scanning_initialised = False
    pass            



def illustrate_flag_signals(set_radius_input = 3, xmin = 200, xmax = 2400, ymin = 300, ymax = 4300):
    load()
    initialise_background_scanning(set_radius_input = set_radius_input)
    # flag_signals(xmin = 200, xmax = 2400, ymin = 300, ymax = 4300, sigmas = 5, scan_range = (3350, 3450))
    flag_signals(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, sigmas = 5, scan_range = (3350, 3450))
    
    global signal_flags_imported, local_background_array_imported
    signal_flags_imported = True
    local_background_array_imported = True
    export_flagged_regions()
    export_background()
    
    preview_flagged_regions()
    preview_background()

    pass

#%%

def export_flagged_regions():
    global signal_flags
    export_signal_flags = pd.DataFrame(signal_flags)
    export_signal_flags.to_excel(str('flagged_region %s.xlsx' % datetime.now().strftime("%Y%m%d_%H_%M_%S")))
    pass

signal_flags_imported = False
def import_flagged_regions( filename = 'flagged_region.xlsx'):
    global signal_flags, signal_flags_imported
    signal_flags = pd.read_excel(filename)
    signal_flags = np.transpose(signal_flags)
    signal_flags_imported = True
    pass



def export_background():
    global local_background_array
    export_local_background_array = pd.DataFrame(local_background_array)
    export_local_background_array.to_excel(str('local_background_array_%s.xlsx' % datetime.now().strftime("%Y%m%d_%H_%M_%S")))
    pass

local_background_array_imported = False
def import_background(filename = 'background.xlsx'):
    global local_background_array, local_background_array_imported
    local_background_array = pd.read_excel(filename)
    local_background_array = np.transpose(local_background_array)
    local_background_array_imported = True
    pass


def export_mapping():
    global mapping
    export_mapping_array = pd.DataFrame(mapping)
    export_mapping_array.to_excel(str('mapping_%s.xlsx' % datetime.now().strftime("%Y%m%d_%H_%M_%S")))
    pass

#%%
def plot_mapping():
    
    SMALL_SIZE = 8
    MEDIUM_SIZE = 10
    BIGGER_SIZE = 12
    
    global signal_flags_imported, local_background_array_imported, signal_flags, local_background_array
    
    if signal_flags_imported == False:
        import_flagged_regions('flagged_region 20240108_08_37_51.xlsx')
        pass
    print('Flagged regions imported.')
    if local_background_array_imported == False:
        import_background('local_background_array_20240108_08_39_07.xlsx')
        pass
    print('Background imported.')
    
    from astropy.visualization import imshow_norm, MinMaxInterval, SqrtStretch, SquaredStretch
    
    set_figsize = (18*cm, 16*cm)

    fig = plt.figure(layout = 'tight', figsize = set_figsize)
    gs = fig.add_gridspec(1, 3,  width_ratios = (14, 14, 0.75), wspace = 0)

    left = fig.add_subplot(gs[0])
    right = fig.add_subplot(gs[1], sharey = left)
    colorbar = fig.add_subplot(gs[2])
        
    # im, norm = imshow_norm(signal_flags, left, vmin = 0, vmax = 1)
    imshow_norm(signal_flags, ax = left, cmap = 'binary', vmin = 0, vmax = 1, stretch=SquaredStretch())
    print('Flagged regions plotted.')
    left.set_title('Binary Mask')
    left.set_ylabel('RA (J2000)')
    
    # im = right.imshow(local_background_array, vmin = 3350, vmax = 3450, cmap = 'jet', )
    im, norm = imshow_norm(local_background_array, ax = right, origin='lower', vmin = 3350, vmax = 3450, stretch=SquaredStretch(), cmap = 'viridis')
    print('Background counts plotted.')
    right.set_title('Background Count Mapping')
    
    # ax.tick_params(axis = "x", labelbottom = False)
    right.tick_params(axis = "y", labelleft = False, left = False)
    import matplotlib
    fig.colorbar(mappable = matplotlib.cm.ScalarMappable(norm = colors.Normalize(vmin = 3350, vmax = 3450), cmap = 'viridis'), cax = colorbar)
    colorbar.text(x = 0.20, y = 3430, s = 'Luminosity (ADU)', rotation = -90)
    
    x_ticks = np.linspace(400, 1800, 3)
    x_ticklabels = []
    for i in x_ticks:
        x_ticklabels.append(map_pixel('hms', [10, 46, 00.00], i + 200, 0.258, display = True))
        pass
    left.xaxis.set(ticks = x_ticks, ticklabels = x_ticklabels)
    left.tick_params(axis = "x", rotation = 0)
    right.xaxis.set(ticks = x_ticks, ticklabels = x_ticklabels)
    right.tick_params(axis = "x", rotation = 0)
    

    y_ticks = np.linspace(600, 4000, 4)
    y_ticklabels = []
    for i in y_ticks:
        y_ticklabels.append(map_pixel('dms', [59, 1, 59.99], i + 300, 0.258, display = True))
        pass
    left.yaxis.set(ticks = y_ticks, ticklabels = y_ticklabels)
    left.tick_params(axis = "y", rotation = 90)
    right.tick_params(axis = "y", labelbottom = False, bottom = False)
    
    fig.supxlabel('DEC (J2000)')
    
    export_plot('plot_mapping', figsize = set_figsize, item = fig)
    pass
    
    
    
    
    
    
    
    
    
    
    
    
#%%
    

def extract_local_background(x, y):
    return local_background_array[x - xorigin][y - yorigin]
    


def preview_flagged_regions():
    fig, ax = plt.subplots(nrows = 1, ncols = 1,  gridspec_kw={}, layout = 'tight')
    ax.set_title("Flagged Regions")
    ax.imshow(signal_flags)
    # plt.invert_yaxis()
    export_plot('preview_flagged_regions', item = fig)
    pass

#%%

def set_flagged_zero(x, y):
    global xorigin, yorigin
    signal_flags[x - xorigin][y - yorigin] = 0
    pass



def is_flagged(x, y):
    global xorigin, yorigin
    if signal_flags[x - xorigin][y - yorigin] == 1:
        return True
    if signal_flags[x - xorigin][y - yorigin] == 0:
        return False
    
    

def get_x_cut_flagged(x_mean, y_mean, width = 20):
    
    x_cut = []
    
    for i in range (x_mean - width, x_mean + width):
        x_cut.append(is_flagged(i, y_mean))
        pass
    
    return x_cut

def get_y_cut_flagged(x_mean, y_mean, width = 20):
    
    y_cut = []
    
    for i in range (y_mean - width, y_mean + width):
        y_cut.append(is_flagged(x_mean, i))
        pass
    
    return y_cut



def preview_background():
    plt.figure()
    plt.title("Background Count Mapping")
    plt.imshow(np.flipud(np.transpose(local_background_array)))
    # plt.gca().invert_yaxis()
    # plt.invert_yaxis()
    plt.colorbar()
    export_plot('preview_background', item = plt)
    pass



def get_object_radius_improved(x_mean, y_mean, sigma, width = 20, plot = False):
    
    # x_array = []
    
    # for i in range(x_mean - width, x_mean + width + 1):
    #     x_array.append(sum(get_y_cut(i, y_mean, width = width, plot = False)))
    #     pass
    
    # poptx, pcovx = curve_fit(gaussian_yshift, np.arange(len(x_array)), x_array)
    
    x_array = []
    
    for i in range(x_mean - width, x_mean + width + 1):
        x_array.append(sum(get_y_cut_flagged(i, y_mean, width = width, plot = False)))
        pass
    
    poptx, pcovx = curve_fit(gaussian_yshift, np.arange(len(x_array)), x_array)
    
    # y_array = []
    
    # for i in range(y_mean - width, y_mean + width + 1):
    #     y_array.append(sum(get_x_cut(i, x_mean, width = width, plot = False)))
    #     pass
    
    # popty, pcovy = curve_fit(gaussian_yshift, np.arange(len(y_array)), y_array)
    
    y_array = []
    
    for i in range(y_mean - width, y_mean + width + 1):
        y_array.append(sum(get_x_cut_flagged(i, x_mean, width = width, plot = False)))
        pass
    
    popty, pcovy = curve_fit(gaussian_yshift, np.arange(len(y_array)), y_array)
    
    if abs(poptx[2]) >= abs(popty[2]):
        one_sigma = abs(poptx[2])
        pass
    else:
        one_sigma = abs(popty[2])
    
    if plot == True:
        plt.figure()
        plt.scatter(np.arange(len(x_array)), x_array)
        plt.axvline(poptx[1] - sigma * one_sigma)
        plt.axvline(poptx[1] + sigma * one_sigma)
        export_plot('get_object_radius_improved', item = plt)
        pass
    
    return sigma * one_sigma
        


def get_object_count_sum_improved(x, y, sigma = 0, sigma_signal = 4, sigma_background_min = 4, sigma_background_max = 6):

    # background
    # background_sum_count = 0
    # background_pixel_count = 0
    
    # sigma = get_object_radius_improved(x, y, 1, width = 20, plot = False)
    # radius_annulus_min = sigma * sigma_background_min
    # radius_annulus_max = sigma * sigma_background_max
    
    
    # for i in range(x - int(radius_annulus_min), x + int(radius_annulus_min) + 1):
    #     for j in range(y - int(radius_annulus_min), y + int(radius_annulus_min) + 1):
            
    #         condition = in_annulus(i, j, x, y, radius_annulus_min, radius_annulus_max)
            
    #         if condition == True:
    #             background_sum_count += get_value(i, j)
    #             background_pixel_count += 1
    #             # replace(i, j, extract_local_background(i, j))
    #             set_flagged_zero(i, j)
    #             pass
    #         pass
    #     pass
    
    # average_background = background_sum_count/background_pixel_count
    
    # object
    if sigma == 0:
        # radius = get_object_radius_improved(x, y, sigma_signal, width = 20, plot = False)
        radius = get_object_radius(x, y, sigma_signal, width = 20, plot = False)
        pass
    else:
        radius = sigma * sigma_signal
        pass
    
    sum_count = 0
    pixel_count = 0
    
    for i in range(x - int(radius), x + int(radius) + 1):
        for j in range(y - int(radius), y + int(radius) + 1):
            
            condition = in_circle(i, j, x, y, radius)
            
            if condition == True:
                sum_count += get_value(i, j)
                sum_count -= extract_local_background(i, j)
                pixel_count += 1
                # replace(i, j, extract_local_background(i, j))
                # print(i, j)
                # print(radius)
                set_flagged_zero(i, j)
                pass
            pass
        pass
    
    
    # total_background = average_background * pixel_count
    
    # object_count = sum_count - total_background
    
    object_count = sum_count
    
    return object_count



def object_search_max_improved(xmin, xmax, ymin, ymax):
    
    max_count = 0
    max_location = [-1, -1]
    
    for i in tqdm(range(xmin, xmax + 1), position = 0, leave = True):
        for j in range(ymin, ymax + 1):
            if get_value(i, j) > max_count and is_flagged(i, j) == True:
                max_count = get_value(i, j)
                max_location = [i, j]
                pass
            pass
        pass
    
    return max_location



def initialise_exception_list():
    global exception_list, exception_list_initialised
    exception_list = pd.DataFrame(columns = ['X-Location', 'Y-Location', 'Brightness', 'Sigma', 'Condition'])
    exception_list_initialised = True
    pass



def preview_exception_list():
    return exception_list.to_string()
    
    
    
def catalogue_window_improved(xmin, xmax, ymin, ymax, sigma_signal_input = 4, sigma_background_min_input = 4, sigma_background_max_input = 6, upper_threshold = 36000, sigmas = 5, scan_range = (3350, 3450), plot = False):
    
    global catalogue_initialised
    if catalogue_initialised != True:
        raise Exception("No catalogue detected, initialise using : initialise_catalogue()")
    
    global exception_list_initialised
    if exception_list_initialised != True:
        raise Exception("No catalogue detected, initialise using : initialise_exception_list()")
    
    initialise_background_scanning(set_radius_input = 15)
    flag_signals(xmin - 120, xmax + 120, ymin - 120, ymax + 120, sigmas, scan_range)
    
    if plot == True:
        preview_flagged_regions()
        preview_background()
        pass
        
    run_count = 0
    
    global catalogue
    continue_scanning = True
    
    global exception_list
    
    while continue_scanning == True:
        
        # stopped here
        max_coordinate = object_search_max_improved(xmin, xmax, ymin, ymax)
        
        if max_coordinate == [-1, -1]:
            continue_scanning = False
            pass
        elif max_coordinate != [-1, -1] and is_flagged(max_coordinate[0], max_coordinate[1]) == True:
            
            max_pixel = get_value(max_coordinate[0], max_coordinate[1])
            
            try:
                print(max_coordinate)
                # sigma = get_object_radius_improved(max_coordinate[0], max_coordinate[1], 1, width = 20, plot = False)
                sigma = get_object_radius(max_coordinate[0], max_coordinate[1], 1, width = 20, plot = False)
                
                if sigma >= 1:
                    brightness = get_object_count_sum_improved(max_coordinate[0], max_coordinate[1], sigma, sigma_signal_input, sigma_background_min_input, sigma_background_max_input)
                    
                    if max_pixel < upper_threshold:                    
                        catalogue.loc[len(catalogue)] = [max_coordinate[0], max_coordinate[1], brightness, sigma, sigma_signal_input, sigma_background_min_input, sigma_background_max_input]
                        run_count += 1
                        print(run_count)
                        print(sigma)
                        print('')
                        pass
                    elif max_pixel >= upper_threshold:
                        exception_list.loc[len(exception_list)] = [max_coordinate[0], max_coordinate[1], brightness, sigma, 'Saturated']
                        pass
                    pass
                
                elif sigma < 1:
                    set_flagged_zero(max_coordinate[0], max_coordinate[1])
                    # pd.DataFramecolumns = ['X-Location', 'Y-Location', 'Brightness', 'Sigma', 'Condition']
                    exception_list.loc[len(exception_list)] = [max_coordinate[0], max_coordinate[1], 'N/A', sigma, 'Sigma < 1']
                    pass
                
            except RuntimeError:
                print('RuntimeError: Iteration not converging, trial passed')
                set_flagged_zero(max_coordinate[0], max_coordinate[1])
                exception_list.loc[len(exception_list)] = [max_coordinate[0], max_coordinate[1], 'N/A', 'N/A', 'RuntimeError: Iteration not converging']
                pass
            except KeyError:
                print('KeyError: Iteration not converging, trial passed')
                set_flagged_zero(max_coordinate[0], max_coordinate[1])
                exception_list.loc[len(exception_list)] = [max_coordinate[0], max_coordinate[1], 'N/A', 'N/A', 'KeyError']
                pass
            pass
        pass
    pass



def categorise():
    global mapping
    binary_mask = signal_flags
    # perform connected component analysis
    mapping, count = ski.measure.label(binary_mask, connectivity = 2, return_num = True)
    return count



def preview_mapping():
    global mapping
    plt.figure()
    plt.title('Grouping of Objects')
    plt.imshow(np.flipud(np.transpose(mapping)))
    export_plot('preview_mapping', item = plt)
    pass



def catalogue_window_improved_improved(xmin, xmax, ymin, ymax, upper_threshold = 36000, sigmas = 5, scan_range = (3350, 3450), plot = False):
    
    global catalogue_initialised
    if catalogue_initialised != True:
        raise Exception("No catalogue detected, initialise using : initialise_catalogue()")
    
    global exception_list_initialised
    if exception_list_initialised != True:
        raise Exception("No catalogue detected, initialise using : initialise_exception_list()")
    
    initialise_background_scanning(set_radius_input = 10)
    flag_signals(xmin, xmax, ymin, ymax, sigmas, scan_range)
    
    export_flagged_regions()
    export_background()
    
    group_number = categorise()
    
    if plot == True:
        preview_flagged_regions()
        preview_background()
        preview_mapping()
        pass
    
    global catalogue, xorigin, yorigin
    catalogue = pd.DataFrame(np.zeros((group_number + 1, 4)), columns = ['X-Location', 'Y-Location', 'Brightness', 'Number of Pixels'])
    
    for i in range(0, len(mapping[:, 0])):
        for j in range(0, len(mapping[0, :])):
            group = mapping[i][j]
            luminosity = get_value(i + xorigin, j + yorigin) - extract_local_background(i + xorigin, j + yorigin)
            catalogue['Brightness'][group] += luminosity
            catalogue['Number of Pixels'][group] += 1
            pass
        pass
    return catalogue
            
            

def load_catalogue(filename = 'Catalogue_x200_2400_y200_4400.xlsx'):
    global mags, luminosities, catalogue, mags_lower, mags_upper
   
    catalogue = pd.read_excel(filename)
    luminosities = catalogue['Brightness']
    mags = []
    mags_lower = []
    mags_upper = []
    removed = 0
    for i in range(len(luminosities)):
        if luminosities[i] > 0:
            mag = get_calibrated_mag(luminosities[i])
            mag_err = np.sqrt((MAGZRR**2) + (2.5 * (1/np.log(10)) * (luminosities[i]**(-0.5)))**2)
            mag_upper = mag + mag_err
            mag_lower = mag - mag_err
            if mag != 0 and mag <= 9e9:
                mags.append(mag)
                mags_lower.append(mag_lower)
                mags_upper.append(mag_upper)
                pass
            pass
        else:
            removed += 1
            pass
        pass
    print('total removed = ', removed)
    pass

    
def magnitudes():
    return mags



def bin_magnitudes(set_range = (10, 15), bins = 10):
    hist, bin_edges = np.histogram(magnitudes(), range = (10, 15), bins = 10)
    return hist, bin_edges
   
   
def N_m(m):
    N = 0
    for i in range(0, int(len(mags))):
        if mags[i] < m:
            N += 1
    return N


def N_m_upper(m):
    N = 0
    for i in range(0, int(len(mags_upper))):
        if mags_upper[i] < m:
            N += 1
    return N


def N_m_lower(m):
    N = 0
    for i in range(0, int(len(mags_lower))):
        if mags_lower[i] < m:
            N += 1
    return N

       
def log_plot(start, stop, number):
    global m_list, log_N_m_list, N_m_list
    m_list = np.linspace(start, stop, number)
    
    # median
    log_N_m_list = []
    N_m_list = []
    for i in range(len(m_list)):
        log_N_m = np.log10(N_m(m_list[i]))
        N_m_list.append(N_m(m_list[i]))
        log_N_m_list.append(log_N_m)
        pass
    
    # upper
    log_N_m_upper_list = []
    for i in range(len(m_list)):
        log_N_m = np.log10(N_m_upper(m_list[i]))
        log_N_m_upper_list.append(log_N_m)
        pass
    
    # lower
    log_N_m_lower_list = []
    for i in range(len(m_list)):
        log_N_m = np.log10(N_m_lower(m_list[i]))
        log_N_m_lower_list.append(log_N_m)
        pass
    
    plt.figure()
    plt.scatter(m_list, log_N_m_list, color = 'black', marker = '.')
    # plt.errorbar(m_list, log_N_m_list, error_final, fmt = 'none', color = 'red')
    error_final = (1/np.log(10)) * (np.sqrt(N_m_list)/N_m_list)
    plt.fill_between(m_list, y1 = log_N_m_list - error_final, y2 = log_N_m_list + error_final, color = 'gray', alpha = 0.3)
    plt.plot(m_list, log_N_m_lower_list, color = 'blue')
    plt.plot(m_list, log_N_m_upper_list, color = 'red')
    plt.grid()
    
    export_plot('log_plot', item = plt)
    pass

def get_squared_angle(add_back = 1598843):
    # remember to check values for different windows
    total_pixel = 2200 * 4200
    removed_pixels = get_removed_count() - add_back
    each_pixel = 0.258**2
    total = (total_pixel - removed_pixels) * each_pixel / (3600)**2
    return total


def get_N_per_half_mag_per_deg2():
    N, bin_edges = bin_magnitudes(set_range = (10, 15), bins = 10)
    norm = N / get_squared_angle(add_back = 1598843)
    return norm, bin_edges


def get_uncertainty():
    values, bin_edges = get_N_per_half_mag_per_deg2()
    out = (1/np.log(10)) * ((values + (values**2))**0.5) / values
    return out

def fit_elliptical_gaussian(xcenter, ycenter, width, print_output = False):
    
    param = astro.fit_gauss_elliptical([xcenter,ycenter], np.array(hdulistdf)[xcenter - width: xcenter + width + 1, ycenter - width : ycenter + width + 1])
    
    if print_output == True:
        print("floor (in ADU)               =", param[1])
        print("height (in ADU)              =", param[2])
        print("x0 (in pixels)               =", param[3])
        print("y0 (in pixels)               =", param[4])
        print("fwhm, smallest (in pixels)   =", param[5])
        print("fwhm, largest (in pixels)    =", param[6])
        print("Degrees, that gives the direction of the largest FWHM, measured starting from x (vertical direction) in the clockwise direction. This angle is between -90 deg and +90 deg. =", param[7])
    
    return param



def fit_circular_gaussian(xcenter, ycenter, width, print_output = False):
    param = astro.fit_gauss_circular([xcenter - width, ycenter - width],np.array(hdulistdf)[xcenter - width: xcenter + width + 1, ycenter - width : ycenter + width + 1])
    
    if print_output == True:
        print("floor (in ADU)   =", param[1])
        print("height (in ADU)  =", param[2])
        print("x0 (in pixels)   =", param[3])
        print("y0 (in pixels)   =", param[4])
        print("fwhm (in pixels) =", param[5])
        pass
    
    return param



def geometric_center(xs, ys):
    
    x_number = 0
    x_sum    = 0
    y_number = 0
    y_sum    = 0
    
    # input check
    if len(xs) != len(ys):
        raise Exception('Length of X and Y is not the same.')
        pass
    
    for i in range(0, len(xs)):
        x_sum    += xs[i]
        x_number += 1
        y_sum    += ys[i]
        y_number += 1
        pass
    
    if y_number != x_number:
        raise Exception('NUmber of X\'s is not the same as number of Y\'s ')
        pass 
    
    x_geometric_center = x_sum/x_number
    y_geometric_center = y_sum/y_number
    
    return x_geometric_center, y_geometric_center



# def catalogue_window_improved_improved_improved(xmin, xmax, ymin, ymax, upper_threshold = 36000, sigmas = 5, scan_range = (3350, 3450), plot = False):
    
#     global catalogue_initialised
#     if catalogue_initialised != True:
#         raise Exception("No catalogue detected, initialise using : initialise_catalogue()")
    
#     global exception_list_initialised
#     if exception_list_initialised != True:
#         raise Exception("No catalogue detected, initialise using : initialise_exception_list()")
    
#     initialise_background_scanning(set_radius_input = 15)
#     flag_signals(xmin, xmax, ymin, ymax, sigmas, scan_range)
    
#     group_number = categorise()
    
#     if plot == True:
#         preview_flagged_regions()
#         preview_background()
#         preview_mapping()
#         pass
    
#     global catalogue, xorigin, yorigin
#     catalogue = pd.DataFrame(np.zeros((group_number + 1, 4)), columns = ['X-Location', 'Y-Location', 'Brightness', 'Number of Pixels'])
    
#     table = pd.DataFrame(columns = ['Object', 'Array of Coordinates [[x1, y1], [x2, y2], ...]', 'radius', 'Luminosity'])
#     for i in range(0, len(mapping[:, 0])):
#         for j in range(0, len(mapping[0, :])):
#             group = mapping[i][j]
#             luminosity = get_value(i + xorigin, j + yorigin) - extract_local_background(i + xorigin, j + yorigin)
#             catalogue['Brightness'][group] += luminosity
#             catalogue['Number of Pixels'][group] += 1
#             pass
#         pass
#     return catalogue


#%%
def methodology():
    # coordinates (935, 1620), (1000, 1690)
    # xmin = 935
    # xmax = 1000
    # ymin = 1620
    # ymax = 1690
    
    x_mean = 968
    y_mean = 1655
    width = 22
    xmin  = x_mean - width
    xmax  = x_mean + width
    ymin  = y_mean - width
    ymax  = y_mean + width
    
    load()
    maxcoord = object_search_max(xmin, xmax, ymin, ymax, background_threshold = 3416)
    
    set_figsize = (9*cm, 8.5*cm)
    
    fig = plt.figure(figsize = set_figsize, layout = 'tight')
    gs = fig.add_gridspec(2, 3,  width_ratios = (3, 10, 0.9), height_ratios = (10, 3), wspace = 0, hspace = 0)
    ax = fig.add_subplot(gs[0, 1], box_aspect = 1)
    colorbar = fig.add_subplot(gs[0, 2])
    ax_gausy = fig.add_subplot(gs[0, 0], sharey = ax)
    ax_gausx = fig.add_subplot(gs[1, 1], sharex = ax)
    
    im = ax.imshow(hdulist[0].data[:, :], cmap = 'jet', vmin = 3416, vmax = get_value(maxcoord[0], maxcoord[1]))
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])
    ax.axvline(maxcoord[0], linestyle = 'dashed', color = 'black', label = 'Cross-section')
    ax.axhline(maxcoord[1], linestyle = 'dashed', color = 'black')
    ax.tick_params(axis = "x", labelbottom = False, bottom = False)
    ax.tick_params(axis = "y", labelleft = False, left = False)
    ax.grid(alpha = 0.2)

    
    cbar = fig.colorbar(im, cax = colorbar)
    cbar.formatter.set_powerlimits((0, 0))
    colorbar.tick_params(axis = "y", direction = 'in', rotation = 0, pad = -12, width = 0.5)
    # colorbar.yaxis.set_label_position("left")
    colorbar.yaxis.tick_left()
    colorbar.set_ylabel('Luminosity (ADU)', fontsize = LATEX_SIZE_SMALL, rotation = - 90, labelpad = 8)
    
    
    ax.scatter(maxcoord[0], maxcoord[1], marker = 'x', color = 'black', s = 40, label = 'Max. luminosity pixel')
    
    radius1 = get_object_radius(maxcoord[0], maxcoord[1], 4, width = width, plot = False)
    radius2 = get_object_radius(maxcoord[0], maxcoord[1], 6, width = width, plot = False)
    circle1 = plt.Circle((maxcoord[0], maxcoord[1]), radius1, facecolor = 'None', edgecolor = 'limegreen', label = 'Object ($4 \sigma$)', linewidth = 0.8)
    circle2 = plt.Circle((maxcoord[0], maxcoord[1]), radius2, facecolor = 'None', edgecolor = 'red', linestyle = 'dashed', label = 'Background ($6 \sigma$)', linewidth = 0.8)
    ax.add_patch(circle1)
    ax.add_patch(circle2)
    ax.legend(loc = 'upper right')
    
    ax_gausx.plot(np.arange(start = xmin, stop = xmax, step = 1), get_x_cut(x_mean, y_mean, width, plot = False), zorder = 10)
    ax_gausx.axvline(maxcoord[0] - radius2, color = 'red', linestyle = 'dashed')
    ax_gausx.axvline(maxcoord[0] - radius1, color = 'limegreen')
    ax_gausx.axvline(maxcoord[0], color = 'black')
    ax_gausx.axvline(maxcoord[0] + radius1, color = 'limegreen')
    ax_gausx.axvline(maxcoord[0] + radius2, color = 'red', linestyle = 'dashed')
    ax_gausx.set_xlim([xmin, xmax - 1])
    ax_gausx.set_ylim([0, 35000])
    ax_gausx.yaxis.set_label_position("right")
    # ax_gausx.yaxis.tick_right()
    ax_gausx.tick_params(axis = "y", direction = 'in', pad = -12, rotation = 0, width = 0.5)
    ax_gausx.tick_params(axis = "x", width = 0.5)
    # ax_gausx.set_ylabel('ADU/$10^4$', rotation = -90, fontsize = LATEX_SIZE_SMALL)
    # y_labels = ax_gausx.get_yticks()
    # ax_gausx.set_yticklabels([y/10000 for y in y_labels])
    ax_gausx.grid()
    
    ax_gausy.plot(get_y_cut(x_mean, y_mean, width, plot = False), np.arange(start = ymin, stop = ymax, step = 1), zorder = 10)
    ax_gausy.axhline(maxcoord[1] - radius2, color = 'red', linestyle = 'dashed')
    ax_gausy.axhline(maxcoord[1] - radius1, color = 'limegreen')
    ax_gausy.axhline(maxcoord[1], color = 'black')
    ax_gausy.axhline(maxcoord[1] + radius1, color = 'limegreen')
    ax_gausy.axhline(maxcoord[1] + radius2, color = 'red', linestyle = 'dashed')
    ax_gausy.set_ylim([ymin, ymax - 1])
    ax_gausy.set_xlim([0, 35000])
    ax_gausy.tick_params(axis = "x", direction = 'in', pad = -10, width = 0.5)
    ax_gausy.tick_params(axis = "y", width = 0.5)
    # ax_gausy.xaxis.set_label_position("top")
    # ax_gausy.xaxis.tick_top()
    ax_gausy.set_xlabel('ADU/$10^4$', fontsize = LATEX_SIZE_SMALL)
    # ax_gausy.xaxis.set_major_formatter(FormatStrFormatter('%2.0e'))
    # x_labels = ax_gausy.get_xticks()
    # ax_gausy.set_xticklabels([x/10000 for x in x_labels])
    ax_gausy.grid()
    
    
    x_ticks = np.linspace(xmin, xmax, 3)
    x_ticklabels = []
    for i in x_ticks:
        x_ticklabels.append(map_pixel('hms', [10, 46, 00.00], i, 0.258, display = True))
        pass
    ax_gausx.xaxis.set(ticks = x_ticks, ticklabels = x_ticklabels)
    ax_gausx.tick_params(axis = "x", rotation = 0)
    
    y_ticks = np.linspace(10000, 30000, 3)
    y_ticklabels = []
    for i in y_ticks:
        y_ticklabels.append(i/10000)
        pass
    ax_gausx.yaxis.set(ticks = y_ticks, ticklabels = y_ticklabels)
    # ax_gausy.tick_params(axis = "y", rotation = 90)
    
    x_ticks = np.linspace(10000, 30000, 3)
    x_ticklabels = []
    for i in x_ticks:
        x_ticklabels.append(i/10000)
        pass
    ax_gausy.xaxis.set(ticks = x_ticks, ticklabels = x_ticklabels)
    # ax_gausy.tick_params(axis = "x", rotation = 0)
    
    y_ticks = np.linspace(ymin, ymax, 4)
    y_ticklabels = []
    for i in y_ticks:
        y_ticklabels.append(map_pixel('dms', [59, 1, 59.99], i, 0.258, display = True))
        pass
    # y_ticklabels[0] = ''
    ax_gausy.yaxis.set(ticks = y_ticks, ticklabels = y_ticklabels)
    ax_gausy.tick_params(axis = "y", rotation = 90)
    
    
    ax_gausx.set_xlabel('RA (J2000)')
    ax_gausy.set_ylabel('DEC (J2000)')
    
    
    
    export_plot('methodology', item = fig, figsize = set_figsize)
    pass
#%%


def SDSS_filter():
    
    set_figsize = (8*cm, 4*cm)
    
    # https://noirlab.edu/science/filters/kp1018
    data = np.loadtxt('C:\\Users\\kyebchoo\\OneDrive - Imperial College London\\Desktop\\Physics\\Year 3\\PHYS60004 - Third Year Physics Laboratory 2023-2024\\Astro Image Processing\\Astro\\Filter Data - Lambda 9 Spectrophotometer.txt', skiprows = 14)

    fig = plt.figure(figsize = set_figsize, layout = 'tight')
    plt.plot(data[:, 0], data[:, 1], color = 'blue')
    plt.xlabel('Wavelength ($\AA$)')
    plt.ylabel('Transmission (%)')
    plt.xlim(min(data[:, 0]), max(data[:, 0]))
    plt.text
    plt.grid()
    export_plot('SDSS_filter', item = fig, figsize = set_figsize)
    pass


def running(catalogue_name):
    load()
    get_MAGZPT()
    load_catalogue(catalogue_name)
    log_plot(9, 20, 100)
    pass
    


def illustration_on_summing_error(pixel_size = 1, plot_range = [0.01, 5], steps = 0.01):
    
    set_figsize = (9*cm, 7*cm)
    
    r_array = []
    circle_area_array = []
    pixels_array = []
    fractional_difference_array = []
    
    mesh = np.ones(shape = (plot_range[1]*2 + 1, plot_range[1]*2 + 1))
    
    lower_range_scaled = int(plot_range[0]/steps)
    upper_range_scaled = int(plot_range[1]/steps + 1)
    for i in range(lower_range_scaled, upper_range_scaled):
        r = i * steps
        r_array.append(r)
        actual_area = np.pi * (r)**2
        circle_area_array.append(actual_area)
        
        pixel_sum = 0
        for x in range(-plot_range[1], plot_range[1] + 1):
            for y in range(-plot_range[1], plot_range[1] + 1):
                if in_circle(x = x, y = y, x_center = 0, y_center = 0, radius = r) == True:
                    pixel_sum += mesh[x + plot_range[1], y + plot_range[1]]
                    pass
                pass
            pass
        pixels_array.append(pixel_sum)
        fractional_difference_array.append((pixel_sum - actual_area) / actual_area)
        pass
    
    fig, ax = plt.subplots(figsize = set_figsize, nrows = 3, ncols = 1, sharex = True, gridspec_kw={'height_ratios': [2, 1.5, 1.5], 'hspace' : 0}, layout = 'tight')
    ax[0].plot(r_array, np.sqrt(circle_area_array), color = 'blue', label = 'Area within circle')
    ax[0].plot(r_array, np.sqrt(pixels_array), color = 'red', label = 'Pixels within circle')
    # ax[0].plot(r_array, circle_area_array, color = 'blue', label = 'Area within circle')
    # ax[0].plot(r_array, pixels_array, color = 'black', label = 'Pixels within circle')
    ax[0].set_ylabel('$\sqrt{Area}$')
    ax[0].grid(alpha = 0.5)
    ax[0].legend()
    ax[0].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax[0].tick_params(axis = 'x', labelbottom = False, bottom = False)
    ax[0].tick_params(axis = 'y', width = 0.5)
    ax[0].set_ylim(-0.95, 10)
    
    ax[1].plot(r_array, fractional_difference_array, color = 'black', label = '$ERR_{frac} = \dfrac{A_{pixels} - A_{actual}}{A_{actual}}$')
    ax[1].set_ylabel('ERR$_{Frac}$')
    # ax[1].text(x = 7, y = -0.5, s = '$\dfrac{A_{pixels} - A_{actual}}{A_{actual}}$')
    ax[1].grid(alpha = 0.5)
    ax[1].legend()
    ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax[1].set_ylim(-1, 4.5)
    ax[1].tick_params(axis = 'x', labelbottom = False, bottom = False)
    ax[1].tick_params(axis = 'y', width = 0.5)
    
    load()
    load_catalogue()
    ax[2].hist(preview_catalogue()['Sigma'], range = (0, plot_range[1]))
    ax[2].set_ylabel('N$_{KPNO}$')
    ax[2].grid(alpha = 0.5)
    ax[2].tick_params(axis = 'x', width = 0.5)
    ax[2].tick_params(axis = 'y', width = 0.5)
    ax[2].set_ylim(0, 680)
    
    plt.xlabel('Radius of Object (Arbitrary Units)')
    plt.xlim(0, plot_range[1])
    export_plot('illustration_on_summing_error', item = fig, figsize = set_figsize)
    pass
    





def load_file(name = 'Catalogue_x200_2400_y200_4400.xlsx'):
    file = pd.read_excel(name)
    file = file.sort_values(['Brightness'], axis = 'index')
    file = file[file['Brightness'] > 0]
    file.loc[:, "Luminosities"] = get_calibrated_mag(file['Brightness'])
    return file



def log_plot_fit(method = 'annulus', start = 10, stop = 19, number = 18, display = True):
    
    set_figsize = (9*cm, 6*cm)
    
    load()
    
    catalogue_35 = load_file('Catalogue_x1300_2400_y200_2300_3_5_sigma.xlsx')
    catalogue_46 = load_file('Catalogue_x1300_2400_y200_2300.xlsx')
    catalogue_57 = load_file('Catalogue_x1300_2400_y200_2300_5_7_sigma.xlsx')
    
    if method == 'annulus':
        catalogue_46_whole = load_file('Catalogue_x200_2400_y200_4400.xlsx')
        pass
    elif method == 'background':
        catalogue_46_whole = load_file('catalogue manual run 2.xlsx')
        pass
    else:
        raise Exception('Invalid method input (annulus/background).')
        pass

    
    number_35 = []
    number_46 = []
    number_57 = []
    number_46_whole = []
    x_axes    = []
    for i in range(0, number):
        magnitude_current = start + (stop - start) * i/number
        x_axes.append(magnitude_current)
        number_35.append(len(catalogue_35[catalogue_35['Luminosities'] < magnitude_current]))
        number_46.append(len(catalogue_46[catalogue_46['Luminosities'] < magnitude_current]))
        number_57.append(len(catalogue_57[catalogue_57['Luminosities'] < magnitude_current]))
        number_46_whole.append(len(catalogue_46_whole[catalogue_46_whole['Luminosities'] < magnitude_current]))
        pass
    
    max_number_number_35 = max(number_35)
    max_number_number_46 = max(number_46)
    max_number_number_57 = max(number_57)
    
    saturation_error = max(max_number_number_35/max_number_number_46 - 1, max_number_number_46/max_number_number_57 - 1)
    
    quarter_catalogue1 = load_file('Catalogue_x200_1300_y200_2300.xlsx')
    quarter_catalogue2 = load_file('Catalogue_x200_1300_y2300_4400.xlsx')
    quarter_catalogue3 = load_file('Catalogue_x1300_2400_y200_2300.xlsx')
    quarter_catalogue4 = load_file('Catalogue_x1300_2400_y2300_4400.xlsx')
    
    quarter_catalogue1_num = len(quarter_catalogue1)
    quarter_catalogue2_num = len(quarter_catalogue2)
    quarter_catalogue3_num = len(quarter_catalogue3)
    quarter_catalogue4_num = len(quarter_catalogue4)
    
    random_fluctuation = np.std([quarter_catalogue1_num, quarter_catalogue2_num, quarter_catalogue3_num, quarter_catalogue4_num])
    
    random_fluctuation = 4 * random_fluctuation / (quarter_catalogue1_num + quarter_catalogue2_num + quarter_catalogue3_num + quarter_catalogue4_num)
    
    error_final = (1/np.log(10)) * ((np.sqrt(number_46_whole) + np.multiply(number_46_whole, (random_fluctuation + saturation_error))) / number_46_whole)
    
    fit_x_axes = []
    fit_number_46 = []
    fit_lower = 12
    fit_upper = 16
    intervals = (fit_upper - fit_lower) * 2 + 1
    for i in range(0, intervals):
        magnitude_current = fit_lower + (fit_upper - fit_lower) * i/intervals
        fit_x_axes.append(magnitude_current)
        fit_number_46.append(len(catalogue_46_whole[catalogue_46_whole['Luminosities'] < magnitude_current]))
        pass
    
    fit_number_46 = np.log10(fit_number_46)
    
    popt, pcov = curve_fit(line, fit_x_axes, fit_number_46)
    
    # fitarray = []
    # for i in range(0, len(x_axes)):
    #     fitarray.append(line(x_axes[i], popt[0], popt[1]))
    #     pass
    
    max_gradient = (np.log10(number_46_whole[x_axes.index(16)]) + error_final[x_axes.index(16)] - np.log10(number_46_whole[x_axes.index(12)]) + error_final[x_axes.index(12)])/(16 - 12)
    min_gradient = (np.log10(number_46_whole[x_axes.index(16)]) - error_final[x_axes.index(16)] - np.log10(number_46_whole[x_axes.index(12)]) - error_final[x_axes.index(12)])/(16 - 12)
    
    print('Mean Gradient = ', popt[0])
    print('Max Gradient = ', max_gradient)
    print('Min Gradient = ', min_gradient)
    
    if display == True:
    
        fig, ax = plt.subplots(figsize = set_figsize)
        
        scaling_to_unit_deg2 = np.log10(generate_total_window_area(unit = 'sqdeg'))
        
        plt.scatter(x_axes, np.log10(number_46_whole) + scaling_to_unit_deg2, marker = 'x', color = 'blue', label = 'Analysis of KPNO data', zorder = 5)
    
        plt.errorbar(x_axes, np.log10(number_46_whole) + scaling_to_unit_deg2, yerr = error_final, fmt = 'none', color = 'red')
        
        xs = np.linspace(8, 20, 12)
        fitarray = []
        for i in range(0, len(xs)):
            fitarray.append(line(xs[i], popt[0], popt[1]))
            pass
        
        
        plt.plot(xs, fitarray + scaling_to_unit_deg2, color = 'black', linestyle = 'dashed', label = 'Linear fit between %s - %s' %(fit_lower, fit_upper), zorder = 10)
    
        plt.fill_between(np.linspace(fit_lower, fit_upper, 3), -1 * np.ones(shape = (3)), 5 * np.ones(shape = (3)), alpha = 0.15, label = 'Fit region')
        
        plt.ylim(0.5 + scaling_to_unit_deg2, 4.0 + scaling_to_unit_deg2)
        plt.xlim(9.5, 20.5)
    
        plt.legend(framealpha = 1)
        plt.grid()
        plt.xlabel('Magnitude, $r^{\prime}$')
        plt.ylabel('Number of Sources, $log_{10}N(<r^{\prime}) /deg^{2}$')
        plt.text(x = 14.5, y = 1.15 + scaling_to_unit_deg2, s = 'Gradient = $%4.2f^{+%4.2f}_{-%4.2f}$' % (popt[0], max_gradient - popt[0], popt[0] - min_gradient))
        
        export_plot('log_plot_fit', item = fig, figsize = set_figsize)
        pass
    return popt, min_gradient, max_gradient

def luminosity_vs_sigma():
    load()    
    catalogue_object = load_file('Catalogue_x200_2400_y200_4400.xlsx')
    
    plt.figure()
    plt.scatter(catalogue_object['Sigma'], catalogue_object['Luminosities'], color = 'black', s = 1)
    plt.xlabel('Radius (arbitrary units)')
    plt.ylabel('Luminosity (arbitrary units)')
    plt.grid(alpha = 0.5)
    export_plot('luminosity_vs_sigma', item = plt)
    pass

def luminosity_vs_temperature():
    load()    
    catalogue_object = load_file('Catalogue_x200_2400_y200_4400.xlsx')
    catalogue_object.loc[:, "Temperature"] = (catalogue_object['Luminosities']/(4 * np.pi * sci.constants.sigma))**0.25
    
    plt.figure()
    plt.scatter(catalogue_object["Temperature"], catalogue_object['Luminosities'], color = 'black', s = 1)
    plt.xlabel('Temperature (arbitrary units)')
    plt.ylabel('Luminosity (arbitrary units)')
    export_plot('luminosity_vs_temperature', item = plt)
    pass



def luminosity_vs_sigma_temperature(scale = False):
    load()    
    catalogue_object = load_file('Catalogue_x200_2400_y200_4400.xlsx')
    catalogue_object.loc[:, "Temperature"] = (catalogue_object['Luminosities']/(4 * np.pi * sci.constants.sigma * catalogue_object['Sigma']**2))**0.25
    
    
    set_figsize = (9*cm, 9*cm)
    
    fig, ax = plt.subplots(figsize = set_figsize, nrows = 1, ncols = 2, sharey = True, layout = 'tight', gridspec_kw={'wspace' : 0})
    
    # ax[0].scatter(np.log10(catalogue_object['Sigma']*0.258), catalogue_object['Luminosities'], color = 'black', s = 1)
    ax[0].scatter(np.log10(catalogue_object['Sigma']*0.258), catalogue_object['Luminosities'], color = 'black', s = 1)
    ax[0].set_ylim(min(catalogue_object['Luminosities'])*0.95, max(catalogue_object['Luminosities'])*1.02)
    ax[0].set_xlabel('Radius, $log_{10}(\sigma)$ ($arcsec^{\prime\prime}$)')
    ax[0].set_ylabel('Magnitude, $r^{\prime}$')
    ax[0].grid(zorder = 1, alpha = 0.5)
    ax[0].set_xlim(-0.65, 0.85)
    
    
    
    ax[1].scatter(catalogue_object["Temperature"], catalogue_object['Luminosities'], color = 'black', s = (catalogue_object['Sigma']**2)/2.5, facecolor = 'none', edgecolor = 'black')
    ax[1].set_xlabel('$T^{\star}$ (arbitrary units)')
    ax[1].invert_xaxis()
    ax[1].grid(zorder = 1, alpha = 0.5)
    ax[1].tick_params(axis = "y", left = False)
    ax[1].set_xlim(75, 6)
    ax[1].text(x = 70, y = 10.5, s = '$r\propto\sigma$')
    
    export_plot('luminosity_vs_sigma_temperature1', figsize = set_figsize, item = fig)
    
    
    
    
    set_figsize = (9*cm, 9*cm)
    
    fig = plt.figure(figsize = set_figsize)
    if scale == True:
        factor = 6
        plt.scatter(np.log10(catalogue_object['Sigma']*0.258), catalogue_object["Temperature"], marker = 'o', facecolor = 'none', edgecolor = 'black', s = ((catalogue_object['Luminosities']**(factor))/(13**factor))**2)
        plt.text(x = 0.6, y = 55, s = '$r\propto\sigma^{%s}$' % factor)
        pass
    else:
        plt.scatter(np.log10(catalogue_object['Sigma']*0.258), catalogue_object["Temperature"], marker = 'o', facecolor = 'none', edgecolor = 'black', s = 1)
        pass
    plt.xlabel('Radius, $log_{10}(\sigma)$ ($arcsec^{\prime\prime}$)')
    plt.ylabel('$T^{\star}$ (arbitrary units)')
    plt.grid()

    export_plot('luminosity_vs_sigma_temperature2', figsize = set_figsize, item = fig)

    
    pass

#%%
def histogram_results():
    load()    
    catalogue_object = load_file('Catalogue_x200_2400_y200_4400.xlsx')
    catalogue_object.loc[:, "Temperature"] = (catalogue_object['Luminosities']/(4 * np.pi * sci.constants.sigma * catalogue_object['Sigma']**2))**0.25
    
    set_figsize = (18*cm, 5*cm)
    
    fig, ax = plt.subplots(figsize = set_figsize, nrows = 1, ncols = 3, gridspec_kw = {'wspace' : 0}, layout = 'tight', sharey = True)
    
    ax[0].hist(catalogue_object['Luminosities'], bins = 10, range = (10, 20))
    ax[0].set_xlabel('Magnitude, $r^{\prime}$')
    ax[0].set_ylabel('$N_{KPNO}$')
    ax[0].grid(alpha = 0.4)
    ax[0].set_ylim(0, 1400)
    
    
    ax[1].hist(catalogue_object['Sigma']*0.258, bins = 16, range = (0, 8))
    ax[1].set_xlabel('Radius, $\sigma$ ($arcsec^{\prime\prime}$)')
    ax[1].grid(alpha = 0.2)
    ax[1].tick_params(axis = "y", left = False)
    
    
    ax[2].hist(catalogue_object['Temperature'], bins = 12, range = (10, 70))
    ax[2].set_xlabel('$T^{\star}$ (arbitrary units)')
    ax[2].grid(alpha = 0.2)
    ax[2].tick_params(axis = "y", left = False)
    
    
    export_plot('histogram_results', item = fig, figsize = set_figsize)
    pass


#%%
def object_mapping(method = 'annulus'):
    load()
    
    if method == 'annulus':
        catalogue_object = load_file('Catalogue_x200_2400_y200_4400.xlsx')
        pass
    elif method == 'background':
        raise Exception('Method \'background\' not supported at the moment.')
        pass
    else:
        raise Exception('Invalid method input (annulus/background).')
        pass
    
    set_figsize = (18*cm, 12*cm)
    
    fig, ax = plt.subplots(figsize = set_figsize, nrows = 1, ncols = 3, sharey = True, gridspec_kw={'wspace' : 0}, layout = 'tight')
    
    ax[0].scatter(catalogue_object["X-Location"], catalogue_object['Y-Location'], color = 'black', s = 1)
    # ax[0].set_xlabel('RA (J2000)')
    ax[0].set_ylabel('DEC (J2000)')
    
    control_x_ticks = [400, 2200, 3]
    
    x_ticks = np.linspace(control_x_ticks[0], control_x_ticks[1], control_x_ticks[2])
    x_ticklabels = []
    for i in x_ticks:
        x_ticklabels.append(map_pixel('hms', [10, 46, 00.00], i, 0.258, display = True))
        pass
    ax[0].xaxis.set(ticks = x_ticks, ticklabels = x_ticklabels)
    ax[0].tick_params(axis = "x", rotation = 0)
    
    y_ticks = np.linspace(0, 4600, 6)
    y_ticklabels = []
    for i in y_ticks:
        y_ticklabels.append(map_pixel('dms', [59, 1, 59.99], i, 0.258, display = True))
        pass
    y_ticklabels[0] = ''
    # ax.set_yticklabels(ticks = y_ticks, ticklabels = y_ticklabels)
    ax[0].yaxis.set(ticks = y_ticks, ticklabels = y_ticklabels)
    ax[0].tick_params(axis = "y", rotation = 90)  
    ax[0].set_xlim(0, 2570)
    ax[0].set_ylim(0, 4611)
    ax[0].grid(alpha = 0.5, zorder = 1)
    ax[0].set_title('Location Mapping')
    
    
    ax[1].scatter(catalogue_object["X-Location"], catalogue_object['Y-Location'], edgecolor = 'black', s = ((catalogue_object['Luminosities']**(4 + 2))/2500**2)**2, marker = 'o', facecolor = 'none')
    ax[1].xaxis.set(ticks = x_ticks, ticklabels = x_ticklabels)
    ax[1].tick_params(axis = "x", rotation = 0)
    ax[1].tick_params(axis = "y", left = False)
    ax[1].set_xlim(0, 2570)
    ax[1].grid(alpha = 0.5, zorder = 1)
    ax[1].set_title('Luminosity Mapping ($r{\propto}{(r^{\prime})^{6}}$)')
    ax[1].set_xlabel('RA (J2000)')
    
    
    ax[2].scatter(catalogue_object["X-Location"], catalogue_object['Y-Location'], edgecolor = 'black', s = catalogue_object['Sigma']**2, marker = 'o', facecolor = 'none')
    ax[2].xaxis.set(ticks = x_ticks, ticklabels = x_ticklabels)
    ax[2].tick_params(axis = "x", rotation = 0)
    ax[2].tick_params(axis = "y", left = False)
    ax[2].set_xlim(0, 2570)
    ax[2].grid(alpha = 0.5, zorder = 1)
    ax[2].set_title('Radius Mapping ($r\propto\sigma$)')
    
    
    removal = pd.read_excel('Removals.xlsx')
    for i in range(0, len(removal)):
        xmin = removal['x left'][i]
        ymin = removal['y bottom'][i]
        width = removal['x right'][i] - removal['x left'][i]
        height = removal['y top'][i] - removal['y bottom'][i]
        rectangle = plt.Rectangle(xy = (xmin, ymin), width = width, height = height, edgecolor = 'none', fill = True, facecolor = 'blue', alpha = 0.05)
        ax[0].add_patch(rectangle)
        rectangle = plt.Rectangle(xy = (xmin, ymin), width = width, height = height, edgecolor = 'none', fill = True, facecolor = 'blue', alpha = 0.05)
        ax[1].add_patch(rectangle)
        rectangle = plt.Rectangle(xy = (xmin, ymin), width = width, height = height, edgecolor = 'none', fill = True, facecolor = 'blue', alpha = 0.05)
        ax[2].add_patch(rectangle)
        rectangle = plt.Rectangle(xy = (xmin, ymin), width = width, height = height, edgecolor = 'blue', fill = False)
        ax[0].add_patch(rectangle)
        rectangle = plt.Rectangle(xy = (xmin, ymin), width = width, height = height, edgecolor = 'blue', fill = False)
        ax[1].add_patch(rectangle)
        rectangle = plt.Rectangle(xy = (xmin, ymin), width = width, height = height, edgecolor = 'blue', fill = False)
        ax[2].add_patch(rectangle)
        pass  
    
    
    export_plot('object_mapping', item = fig, figsize = set_figsize)
    pass
#%%

def load_Durham():
    
    # ref: http://star-www.dur.ac.uk/~nm/pubhtml/counts/counts.html
    
    global durham
    durham = pd.read_csv('DurhamRCounts.txt', sep=" ")
    cites = pd.DataFrame(columns = ['cite'])
    
    for i in range(0, len(durham)):
        if durham['cite'][i] in cites['cite'].values:
            pass
        else:
            cites.loc[len(cites)] = durham['cite'][i]
            pass
        pass
    return cites



run = np.random.randint(low = 0, high = 100)
def get_marker():
    global run
    run += 1
    markers = ['.', 'o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X']
    colors = plt.cm.tab10.colors
    
    return str(markers[run%len(markers)]), colors[run%len(colors)]




def plot_Durham(method = 'annulus'):
    
    set_figsize = (18*cm, 16*cm)
    
    global durham
    cites = load_Durham()
    # print(cites)
    # print(durham)
    
    
    scaling_to_unit_deg2 = np.log10(generate_total_window_area(unit = 'sqdeg'))
    
    fig, ax = plt.subplots(nrows = 1, ncols = 2, figsize = set_figsize, gridspec_kw={'width_ratios': [4, 1], 'wspace' : 0}, layout = 'tight')
    
    for i in range(0, len(cites)):
        marker, color = get_marker()
        ax[0].scatter(durham.loc[durham['cite'] == cites['cite'][i]]['R'], durham.loc[durham['cite'] == cites['cite'][i]]['N'], label = str(cites['cite'][i]).replace("_", " "), s = 25, marker = marker, color = color)
        pass   
    
    
    if method == 'annulus':
        catalogue_46_whole = load_file('Catalogue_x200_2400_y200_4400.xlsx')
        pass
    elif method == 'background':
        catalogue_46_whole = load_file('catalogue manual run.xlsx')
        pass
    else:
        raise Exception('Invalid method input (annulus/background).')
        pass
    number_46_whole = []
    x_axes    = []
    start = 10
    stop = 30
    x_fit = []
    y_fit = []
    for i in range(start, stop):
        magnitude_current = i
        x_axes.append(magnitude_current)
        number_46_whole.append(len(catalogue_46_whole[catalogue_46_whole['Luminosities'] < magnitude_current]))
        if i > 12 and i < 16:
            x_fit.append(i)
            y_fit.append(len(catalogue_46_whole[catalogue_46_whole['Luminosities'] < magnitude_current]))
        pass
    
    
    marker, color = get_marker()
    ax[0].scatter(x_axes, np.log10(number_46_whole) + scaling_to_unit_deg2, label = 'Imperial College (2023)', s = 25, marker = marker, color = color)
    
    
    handles, labels = ax[0].get_legend_handles_labels()

    
    popt, pcov = curve_fit(line, durham['R'], durham['N'])
    mean_gradient = popt[0]
    max_gradient = popt[0] + np.sqrt(pcov[0, 0])
    min_gradient = popt[0] - np.sqrt(pcov[0, 0])
    plus = max_gradient - mean_gradient
    minus = mean_gradient - min_gradient
    
    
    
    popt2, pcov2 = curve_fit(line, x_fit, np.log10(y_fit))
    # mean_gradient2 = popt2[0]
    # max_gradient2 = popt2[0] + np.sqrt(pcov2[0, 0])
    # min_gradient2 = popt2[0] - np.sqrt(pcov2[0, 0])
    
    popt2, min_gradient2, max_gradient2 = log_plot_fit(method = method, display = False)
    
    mean_gradient2 = popt2[0]
    
    plus2 = max_gradient2 - mean_gradient2
    minus2 = mean_gradient2 - min_gradient2
    
    
    
    xs = np.linspace(6, 32, 10)
    ax[0].plot(xs, line(xs, popt[0], popt[1]), linestyle = 'dashed', color = 'black', linewidth = 1, zorder = 1, label = 'Durham Physics Cosmology Research Galaxy Counts')
    ax[0].plot(xs, line(xs, popt2[0], popt2[1]) + scaling_to_unit_deg2, linestyle = 'dashed', color = 'red', linewidth = 1, zorder = 1, label = 'Imperial College Undergraduate Laboratory Counts')
    ax[0].grid(alpha = 0.5)
    ax[0].set_xlabel('Magnitude, $r^{\prime}$')
    ax[0].set_ylabel('log(N) /$deg^{2}$')
    ax[0].set_xlim(8, 30)
    ax[0].set_ylim(-1, 6)
    ax[0].text(x = 17.7, y = 1.1, s = '$Gradient_{(Durham)} = %5.3f^{+%5.3f}_{-%5.3f}$' % (mean_gradient, plus, minus))
    ax[0].text(x = 10.1, y = 3.1 + scaling_to_unit_deg2, s = '$Gradient_{(Imperial)} = %5.3f^{+%5.3f}_{-%5.3f}$' % (mean_gradient2, plus2, minus2))

    
    
    ax[1].legend(handles, labels, framealpha = 1, loc = 'center')
    ax[1].tick_params(axis = "x", labelbottom = False, bottom = False)
    ax[1].tick_params(axis = "y", labelleft = False, left = False)
    ax[1].spines['top'].set_visible(False)
    ax[1].spines['right'].set_visible(False)
    ax[1].spines['bottom'].set_visible(False)
    # ax[1].spines['left'].set_visible(False)
    ax[1].set_xlim(0, 1)
    
    
    ax[1].plot(xs, line(xs, popt[0], popt[1]), linestyle = 'dashed', color = 'black', linewidth = 1, zorder = 1, label = 'Durham Physics Cosmology Research Galaxy Counts')
    ax[1].plot(xs, line(xs, popt2[0] + scaling_to_unit_deg2, popt2[1]), linestyle = 'dashed', color = 'red', linewidth = 1, zorder = 1, label = 'Imperial College Undergraduate Laboratory Counts')
    handles, labels = ax[1].get_legend_handles_labels()
    
    ax[0].legend(handles, labels, framealpha = 1)
    
    # ax[1].legend(ax = ax[0])
    export_plot('plot_Durham', figsize = set_figsize, item = fig)
    pass

#%%%


def hist_background_method():
    bg_catalogue = load_file('Catalogue manual run 2.xlsx')
    annulus_catalogue = load_file()
    
    set_figsize = (9*cm, 9*cm)
    
    fig, ax = plt.subplots(figsize = set_figsize, nrows = 1, ncols = 1, layout = 'tight')
    
    # ax.hist(bg_catalogue['Luminosities'], bins = 16, range = (10, 18), log = True, histtype = 'step')
    n, bins = np.histogram(bg_catalogue['Luminosities'], bins = 16, range = (10, 16))
    n2, bins2 = np.histogram(annulus_catalogue['Luminosities'], bins = 16, range = (10, 16))
    xs = []
    for i in range(0, len(bins) - 1):
        xs.append(np.average([bins[i], bins[i + 1]]))
        pass
    
    ys = []
    ys2 = []
    for i in range(0, len(n)):
        ys.append(np.log10(n[i]) + np.log10(generate_total_window_area(unit = 'sqdeg')))
        ys2.append(np.log10(n2[i]) + np.log10(generate_total_window_area(unit = 'sqdeg')))
        pass
    
    ax.step(xs, ys, where = 'mid', color = 'black', label = 'Background Method')
    ax.step(xs, ys2, where = 'mid', color = 'red', label = 'Annulus Method')
    
    ax.set_xlabel('Magnitude, $r^{\prime}$')
    ax.set_ylabel('$log_{10}(N) /(0.5 mag deg^{2})$')
    
    popt, pcov = curve_fit(line, xs, ys)
    popt2, pcov2 = curve_fit(line, xs, ys2)
    
    mean_gradient = popt[0]
    max_gradient = popt[0] + np.sqrt(pcov[0, 0])
    min_gradient = popt[0] - np.sqrt(pcov[0, 0])
    
    mean_gradient2 = popt2[0]
    max_gradient2 = popt2[0] + np.sqrt(pcov2[0, 0])
    min_gradient2 = popt2[0] - np.sqrt(pcov2[0, 0])
    
    ax.plot(np.linspace(5, 25, 20), line(np.linspace(5, 25, 20), popt[0], popt[1]), linestyle = 'dashed', label = 'Line fit for background method', color ='blue')
    ax.plot(np.linspace(5, 25, 20), line(np.linspace(5, 25, 20), popt2[0], popt2[1]), linestyle = 'dashed', label = 'Line fit for annulus method', color ='green')
    
    ax.text(x = 13.5, y = -0.60, s = 'Gradient = $%4.2f^{+%4.2f}_{-%4.2f}$' % (mean_gradient, max_gradient - mean_gradient, mean_gradient - min_gradient))
    ax.text(x = 9.5, y = 0.5, s = 'Gradient = $%4.2f^{+%4.2f}_{-%4.2f}$' % (mean_gradient2, max_gradient2 - mean_gradient2, mean_gradient2 - min_gradient2))
    
    
    ax.grid(alpha = 0.4)
    ax.set_xlim(9, 17)
    ax.set_ylim(-1, 1.75)
    ax.legend()
    
    export_plot('hist_background_method', figsize = set_figsize, item = fig)
    pass



def plot_Durham2():
    
    set_figsize = (18*cm, 16*cm)
    
    global durham
    cites = load_Durham()
    # print(cites)
    # print(durham)
    
    
    scaling_to_unit_deg2 = np.log10(generate_total_window_area(unit = 'sqdeg'))
    
    fig, ax = plt.subplots(nrows = 1, ncols = 2, figsize = set_figsize, gridspec_kw={'width_ratios': [4, 1], 'wspace' : 0}, layout = 'tight')
    
    for i in range(0, len(cites)):
        marker, color = get_marker()
        ax[0].scatter(durham.loc[durham['cite'] == cites['cite'][i]]['R'], durham.loc[durham['cite'] == cites['cite'][i]]['N'], label = str(cites['cite'][i]).replace("_", " "), s = 25, marker = marker, color = color)
        pass   
    
    bg_catalogue = load_file('Catalogue manual run 2.xlsx')
    annulus_catalogue = load_file()
    
    # ax.hist(bg_catalogue['Luminosities'], bins = 16, range = (10, 18), log = True, histtype = 'step')
    n, bins = np.histogram(bg_catalogue['Luminosities'], bins = 16, range = (10, 16))
    n2, bins2 = np.histogram(annulus_catalogue['Luminosities'], bins = 16, range = (10, 16))
    xs = []
    for i in range(0, len(bins) - 1):
        xs.append(np.average([bins[i], bins[i + 1]]))
        pass
    
    ys = []
    ys2 = []
    for i in range(0, len(n)):
        ys.append(np.log10(n[i]) + np.log10(generate_total_window_area(unit = 'sqdeg')))
        ys2.append(np.log10(n2[i]) + np.log10(generate_total_window_area(unit = 'sqdeg')))
        pass
    
    ax[0].step(xs, ys, where = 'mid', color = 'blue', label = 'Imperial (2023) (BG)')
    ax[0].step(xs, ys2, where = 'mid', color = 'red', label = 'Imperial (2023) (AN)')
    
    
    
    
    handles, labels = ax[0].get_legend_handles_labels()
    
    
    
    popt1, pcov1 = curve_fit(line, xs, ys)
    popt2, pcov2 = curve_fit(line, xs, ys2)
    
    mean_gradient1 = popt1[0]
    max_gradient1 = popt1[0] + np.sqrt(pcov1[0, 0])
    min_gradient1 = popt1[0] - np.sqrt(pcov1[0, 0])
    
    mean_gradient2 = popt2[0]
    max_gradient2 = popt2[0] + np.sqrt(pcov2[0, 0])
    min_gradient2 = popt2[0] - np.sqrt(pcov2[0, 0])

    
    popt, pcov = curve_fit(line, durham['R'], durham['N'])
    mean_gradient = popt[0]
    max_gradient = popt[0] + np.sqrt(pcov[0, 0])
    min_gradient = popt[0] - np.sqrt(pcov[0, 0])
    plus = max_gradient - mean_gradient
    minus = mean_gradient - min_gradient
    
    
    

    
    
    
    
    
    ax[0].grid(alpha = 0.5)
    ax[0].set_xlabel('Magnitude, $r^{\prime}$')
    ax[0].set_ylabel('$log_{10}(N) /\;(0.5\; mag^{-2}\; deg^{-2})$')
    
    
    ax[0].set_xlim(8, 30)
    ax[0].set_ylim(-1, 6)
    ax[0].text(x = 11.5, y = 2.5, s = '$Gradient_{Durham} = %5.3f^{+%5.3f}_{-%5.3f}$' % (mean_gradient, plus, minus))
    ax[0].text(x = 19.0, y = 1.5, s = '$Gradient_{Imperial(BG)}$ = $%4.2f^{+%4.2f}_{-%4.2f}$' % (mean_gradient1, max_gradient1 - mean_gradient1, mean_gradient1 - min_gradient1), color = 'blue')
    ax[0].text(x = 8.5, y = 1.0, s = '$Gradient_{Imperial(AN)}$ = $%4.2f^{+%4.2f}_{-%4.2f}$' % (mean_gradient2, max_gradient2 - mean_gradient2, mean_gradient2 - min_gradient2), color = 'red')
    
    
    handles, labels = ax[0].get_legend_handles_labels()
    
    
    ax[1].legend(handles, labels, framealpha = 1, loc = 'center')
    ax[1].tick_params(axis = "x", labelbottom = False, bottom = False)
    ax[1].tick_params(axis = "y", labelleft = False, left = False)
    ax[1].spines['top'].set_visible(False)
    ax[1].spines['right'].set_visible(False)
    ax[1].spines['bottom'].set_visible(False)
    # ax[1].spines['left'].set_visible(False)
    ax[1].set_xlim(0, 1)
    

    xs = np.linspace(6, 32, 10)
    ax[0].plot(xs, line(xs, popt[0], popt[1]), linestyle = 'dashed', color = 'black', linewidth = 1, zorder = 1, label = 'Durham Physics Cosmology Research Galaxy Counts')
    ax[0].plot(np.linspace(5, 30, 20), line(np.linspace(5, 30, 20), popt1[0], popt1[1]), linestyle = 'dashed', color ='blue', label = 'Imperial College Undergraduate Laboratory Counts')
    ax[0].plot(np.linspace(5, 30, 20), line(np.linspace(5, 30, 20), popt2[0], popt2[1]), linestyle = 'dashed', color ='red', label = 'Imperial College Undergraduate Laboratory Counts')
    ax[1].plot(xs, line(xs, popt[0], popt[1]), linestyle = 'dashed', color = 'black', linewidth = 1, zorder = 1, label = 'Durham Physics Cosmology Research Galaxy Counts')
    ax[1].plot(np.linspace(5, 30, 20), line(np.linspace(5, 30, 20), popt1[0], popt1[1]), linestyle = 'dashed', color ='blue', label = 'Imperial College Undergraduate Laboratory Counts')
    ax[1].plot(np.linspace(5, 30, 20), line(np.linspace(5, 30, 20), popt2[0], popt2[1]), linestyle = 'dashed', color ='red', label = 'Imperial College Undergraduate Laboratory Counts')
    handles, labels = ax[1].get_legend_handles_labels()
    

    
    ax[0].legend(handles, labels, framealpha = 1)
    export_plot('plot_Durham2', figsize = set_figsize, item = fig)
    pass


def convert_to_alpha(gradient):
    alpha = gradient/0.4
    return alpha
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
#%%

# # load()
# # # fit_background()
# # remove_pixels()
# # preview_2D()

# # export()

# # calibrate_image()
# # preview_2D()


# initialise_catalogue()
# load(filename = "Astro\Fits_Data\mosaic.fits")
# # load_excel('HDUlist_replaced_with_3416.xlsx')
# # calibrate_image()
# remove_pixels(filename = 'Removals.xlsx', replace_value = 3416)
# catalogue_window(200, 2400, 1301, 2300, sigma_signal_input = 4, sigma_background_min_input = 4, sigma_background_max_input = 6, background_threshold = 3500)
# # preview_2D()
# # export()
# # export_catalogue()

# # plt.hist(preview_catalogue()['Brightness'])






# plt.figure()
# plt.hist(preview_catalogue()['Brightness'])
# plt.show()




# plt.figure()
# plt.scatter(preview_catalogue()['X-Location'], preview_catalogue()['Y-Location'])
# plt.show()


#%%

# load()
# remove_pixels()
# initialise_catalogue()
# initialise_background_scanning()
# initialise_exception_list()
# catalogue_window_improved_improved(200, 2400, 200, 4400, plot = True)


# from astropy.stats import SigmaClip
# from photutils.background import Background2D, MedianBackground

# sigma_clip = SigmaClip(sigma=3.0)
# bkg_estimator = MedianBackground()
# bkg = Background2D(hdulistdf, (30, 30), filter_size=(5, 5),
#                    sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
# print(bkg.background_median)
# print(bkg.background_rms_median)
# plt.figure()
# plt.imshow(bkg.background, origin='lower', cmap='jet',
#            interpolation='nearest')
# plt.show()
