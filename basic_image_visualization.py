# Script to load, perform basic postprocessing, and plot data analysis for the 
# GCaMP mice data.

import sys
import os
sys.path.append("./utilities")
from pylab import *
from scipy.stats import pearsonr
from matplotlib.animation import FuncAnimation
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

from chungotools import *
from fluotools import ImStack 

def basic_processing(file_in, figure_out, video_out, stim_period, us_start_id, selected_point, exclude_stim, cycle_offset_id, midline_id, hemisphere_sign, brain_outline_file, unattended_flag = 0):
    """
    This function performs basic processing and plotting. Cycle averaging is 
    performed and cycle exclusion is supported. A video of the cycle averaged
    sequence is written. The unattended_flag turns on and off the pause dialogs
    to inspect the time traces and annotate the cycles to be excluded.
    
    file_in: full path and name of the file to be opened.
    figure_out: full path and name of the saved figures. It is automatically modified
        if there are cycles to be excluded.
    video_out: full path and name of the video (must include '.mp4'). It is 
        automatically modified if there are cycles to be excluded.
    stim_period: Period of the stimulation (s).
    us_start_id: Manual adjustment for the exact begin of the stimulation.
    selected_point: Array with the coordinates of the points to extract the time
        traces. array([[x1, y1], [x2, xy2]) 
    exclude_stim: List containing the cycles to be excluded. If empty, all cycles
        are included.
    cycle_offset_id: Number of samples to be considered as baseline before the
        stimulation begins.
    brain_outline_file
    unattended_flag: 0 (default) for interactive dialog, 1 for unattended processing.    
    """
    
    if len(exclude_stim) != 0:
        figure_out = figure_out + '_exclude'
    
    # (2.0) Normal analysis
    fluo = ImStack(file_in)
    fluo.brain_outline_file = brain_outline_file
    fluo.load_image_stack()
    print('The image stack size is : '+ repr(shape(fluo.image_stack)))
    fluo.stimulation_start_id = int(fluo.stim_frames[0] - us_start_id)
    fluo.spatial_filter_gauss(gauss_sigma = [3, 3])
    fluo.spot_size = [1, 1]
    fluo.mask_threshold = 0.0
    # Correcting for bias
    #fluo.image_stack = fluo.image_stack - 2*median(fluo.image_stack.flatten())
    # Bilateral points
    n_points = shape(selected_point)[0]
    bilateral_points = zeros((2*n_points, 2))
    for nn in range(n_points):
        bilateral_points[2*nn, :] = array([selected_point[nn, 0], midline_id - (selected_point[nn, 1] - midline_id)])
        bilateral_points[2*nn + 1, :] = selected_point[nn,:]
    
    fluo.selected_point = bilateral_points 
    fluo.exclude_stim = exclude_stim 
    
    ## (3.0) Inter-hemisphere 
    inter = ImStack('lol.lol')
    inter.sampling_freq = fluo.sampling_freq
    inter.stimulation_start_id = fluo.stimulation_start_id
    #inter.emission_pulse_duration = fluo.emission_pulse_duration # (s)
    inter.spot_size = fluo.spot_size 
    inter.mask_threshold = fluo.mask_threshold**2
    inter.stim_frames = fluo.stim_frames
    inter.selected_point = selected_point - [0, midline_id]
    inter.exclude_stim = exclude_stim
    
    ## (3.1) Robust hemisphere substraction
    inter_size = array([midline_id, shape(fluo.image_stack)[1] - midline_id])
    if inter_size[0] < inter_size[1]:
        left = fliplr(fluo.image_stack[:,:midline_id])
        right = fluo.image_stack[:, midline_id:inter_size[0] + midline_id]
        
        inter.image_stack = hemisphere_sign*(left - right)
        inter.ref_image = fluo.ref_image[:, midline_id:inter_size[0] + midline_id] \
                        *fliplr(fluo.ref_image[:,:midline_id])
        inter.brain_outline_mask = fluo.brain_outline_mask[:, midline_id:inter_size[0] \
                            + midline_id]*fliplr(fluo.brain_outline_mask[:,:midline_id])
        inter.brain_segment_mask = fluo.brain_outline_mask[:, midline_id:inter_size[0] \
                            + midline_id]*fliplr(fluo.brain_outline_mask[:,:midline_id])
    else:
        left = fliplr(fluo.image_stack[:,midline_id - inter_size[1]:midline_id])
        right = fluo.image_stack[:, midline_id:inter_size[1] + midline_id]
        
        inter.image_stack = hemisphere_sign*(left - right)
        inter.ref_image = fluo.ref_image[:, midline_id:inter_size[1] + midline_id] \
                        *fliplr(fluo.ref_image[:,midline_id - inter_size[1]:midline_id])
        inter.brain_outline_mask = fluo.brain_outline_mask[:, midline_id:inter_size[1] \
                            + midline_id]*fliplr(fluo.brain_outline_mask[:,midline_id - inter_size[1]:midline_id])
        inter.brain_segment_mask = fluo.brain_outline_mask[:, midline_id:inter_size[1] \
                            + midline_id]*fliplr(fluo.brain_outline_mask[:,midline_id - inter_size[1]:midline_id])
    
    ## (4.0) Common processing for visualization
    fluo.get_time_traces()
    inter.get_time_traces()
    
    # (4.1) Reference figure
    pageSetup(250, 0.5)
    figure()
    subplot(1, 2, 1)
    imshow(fluo.ref_image, cmap='gray')
    plot([midline_id, midline_id], [0, fluo.ref_image.shape[0]], '--r')
    subplot(1, 2, 2)
    imshow(inter.ref_image, cmap='gray')
    fig_name = figure_out  + 'ref_both.svg'
    savefig(fig_name)
    frame_id = argmax(inter.time_traces[-3, :])
                       
    fluo.polarity = 'monopolar'
    fluo.plot_visibility = True
    fluo.overlay_transparency = 'logistic'
    fluo.show_time_traces(frame_id = frame_id)
    fig_name = figure_out + 'raw' + repr(frame_id)  + '.svg'
    savefig(fig_name)
    
    # if unattended_flag == 0:
    #     lol = input("Check the plot, if it is fine, press ENTER")
    # close('all')
        
    ## (4.2) Cycle averaging and basic plotting
    fluo.mean_cycle(offset = fluo.stimulation_start_id/fluo.sampling_freq, 
                prf = 1/stim_period,
                cycle_offset_id = cycle_offset_id)

    inter.mean_cycle(offset = fluo.stimulation_start_id/fluo.sampling_freq, 
                prf = 1/(stim_period),
                pulse_duration = stim_duration,
                cycle_offset_id = cycle_offset_id)
       
    ## (4.3) Recalculating the base line after cycle averaging
    new_base = mean(fluo.Cycle.mean[:,:,0:cycle_offset_id], axis = 2)
    fluo.Cycle.mean = fluo.Cycle.mean - repeat(new_base[:,:,newaxis], fluo.Cycle.mean.shape[2], axis = 2)
    
    new_base = mean(inter.Cycle.mean[:,:,0:cycle_offset_id], axis = 2)
    inter.Cycle.mean = inter.Cycle.mean - repeat(new_base[:,:,newaxis], inter.Cycle.mean.shape[2], axis = 2)
    
    fluo.get_time_traces()
    inter.get_time_traces()
    
    ## (4.4) Plotting cycle averaged traces
    frame_id1 = argmax(inter.Cycle.time_traces[-1,:])
    fluo.polarity = 'bipolar'
    fluo.plot_visibility = True
    fluo.overlay_transparency = 'logistic'
    # Minimum position
    frame_id0 = cycle_offset_id + 3
    fluo.show_time_traces(frame_id = frame_id0)
    fig_name = figure_out + 'cycle' + repr(frame_id0)  + '.svg'
    savefig(fig_name)
    # Maximum position
    #frame_id1 = cycle_offset_id + 15
    fluo.show_time_traces(frame_id = frame_id1)
    fig_name = figure_out + 'cycle' + repr(frame_id1)  + '.svg'
    savefig(fig_name)
    
    inter.polarity = 'bipolar'
    inter.plot_visibility = True
    inter.overlay_transparency = 'tukey'
    # Minimum position
    inter.show_time_traces(frame_id = frame_id0)
    fig_name = figure_out + 'cycle' + repr(frame_id0)  + 'inter.svg'
    savefig(fig_name)
    # Maximum position
    inter.show_time_traces(frame_id = frame_id1)
    fig_name = figure_out + 'cycle' + repr(frame_id1)  + 'inter.svg'
    savefig(fig_name)
    
    # close('all')
   
# (0.0) Path and folders to dave the data/figures
## (1.0) Loading and basic postprocessing/detrending
# (1.1) File names of the main datasets

folder_out = '../output/'
folder_name = '../output/' 

###### ###=============================================================================
sub_folder = ''

file_name = 'single_focus_us_stim_kalman_pretreated'
brain_outline_file = '../data/masks/masks_single/m1_brain_outline_file.png'
file_in = folder_name + sub_folder + file_name + '.mat'
figure_out = folder_out + sub_folder + file_name 
video_out = folder_out + sub_folder + file_name + '_cycle.mp4'

us_start_id = 1 # Manual fine adjustment in case the automated tool is shifted
cycle_offset_id = 19 # Purely for visualization purposes
stim_period = 10 # between stimulations, regardless of the steering
period_shift = 0 # For visualization of the raw data
stim_duration = 150e-3 # (s)
midline_id = 118
hemisphere_sign = -1 # plus or minus

## Check the OA reference figure to know where is the US delivery spot. Select also
## the contralateral points [AUD, Control, STIM]

# hTUS paper figure 2
selected_point = array([[143, 192], [150, 161], [75, 147], [98, 147]])  

basic_processing(file_in, figure_out, video_out, stim_period, us_start_id,
                 selected_point, [1], cycle_offset_id, midline_id,
                 hemisphere_sign, brain_outline_file, unattended_flag = 1)

######=============================================================================
sub_folder = ''
file_name = 'multiple_foci_us_stim_kalman_pretreated'
brain_outline_file = '../data/masks/masks_single/m1_brain_outline_file.png'
file_in = folder_name + sub_folder + file_name + '.mat'
figure_out = folder_out + sub_folder + file_name 
video_out = folder_out + sub_folder + file_name + '_cycle.mp4'

us_start_id = 2 # Manual fine adjustment in case the automated tool is shifted
cycle_offset_id = 19 # Purely for visualization purposes
stim_period = 10 # between stimulations, regardless of the steering
period_shift = 3 # For visualization of the raw data
stim_duration = 150e-3 # (s)
midline_id = 118
hemisphere_sign = -1 # plus or minus

## Check the OA reference figure to know where is the US delivery spot. Select also
## the contralateral points
#selected_point = array([[120, 170], [98, 147], [89, 132], [103, 132]])  

# hTUS paper figure 2
selected_point = array([[143, 192], [98, 147], [90, 132], [106, 132], [98, 137]])  

basic_processing(file_in, figure_out, video_out, stim_period, us_start_id,
                 selected_point, [], cycle_offset_id, midline_id,
                 hemisphere_sign, brain_outline_file, unattended_flag = 0)

