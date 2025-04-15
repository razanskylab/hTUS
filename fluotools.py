# Set of classes and functions to operate and modify postprocessed image stacks
# post processed previously in matlab.

from ast import IsNot
from pylab import *
from chungotools import *
import pickle
import h5py
import csv
import time
import matplotlib.image as img
#import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from scipy.fft import fft, fftfreq
from matplotlib.animation import FuncAnimation
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from scipy.signal import savgol_filter
import skimage.segmentation as seg
from scipy.stats import sem

class MeanCycle:
  pass

class ImStack:

    def __init__(self, file_name):
        self.file_name = file_name
        self.plot_visibility = False
        self.overlay_transparency = 'linear'
        self.c_map_name = 'PRGn'
        self.polarity = 'bipolar' # Monopolar, bipolar, asymm [low, high]
        self.vis_threshold = 100  # Threshold for visualization of images 
        

    def load_image_stack(self):
        # (1.0) Loading the data
        dataFile = h5py.File(self.file_name, 'r')
        print('File opened ... ')
        # (1.1) Extracting the data
        PostPro = dataFile['PostPro']
        refImage = dataFile['PostPro']['BaseLine'][:,:]
        self.ref_image = refImage.transpose()
        dF_F = dataFile[PostPro['dF_F'][0, 0]][:,:,:]
        self.data = dF_F.transpose() # To keep the original data
        self.image_stack = dF_F.transpose()
        self.sampling_freq = dataFile['param']['Data']['SamplingFreq'][0, 0]
        try:
           self.stim_frames = sort(dataFile['param']['StimulationTimings'][:][0])
           self.stim_frames = self.stim_frames\
                         + dataFile['param']['Crop']['BaseLineRange'][:][0] + 1 # Minus croped frames
           self.stim_frames = self.stim_frames.astype('int')
           print(repr(self.stim_frames))
        except KeyError:
           self.stim_frames = array([0])
        

        # Brain outline mask
        try:
           lol = img.imread(self.brain_outline_file)[:,:,3]==0
 
           print('########################### read brain_outline_file #####################')
        except FileNotFoundError:
        #except:
           lol = ones(self.ref_image.shape)
        
        self.brain_outline_mask = lol


        # Brain segmentation mask
        try:
           lol = img.imread(self.brain_segment_file) ####
           print('########################### read brain_segment_file #####################')

        except AttributeError:
        #except:
           lol = ones(self.ref_image.shape)
        
        self.brain_segment_mask = lol
            
        print('data loaded')
    
    def reinit_stack(self):
        self.image_stack = self.data
    
    def spatial_filter_gauss(self, gauss_sigma = [1., 1.]):
        try:
            if (self.gauss_sigma != gauss_sigma) and (gauss_sigma != [1., 1.]):
                self.gauss_sigma = gauss_sigma # (s)
        except AttributeError:
            self.gauss_sigma = gauss_sigma
            
        for n_frame in range(self.image_stack.shape[-1]):
            self.image_stack[:,:, n_frame] = gaussian_filter(
                self.image_stack[:,:, n_frame], 
                sigma = gauss_sigma)
        print('Spatial gaussian filtering done!')
     
    def substrac_spatial_average(self):
        # Substracts the mean of the last array dimension
        # For image stacks x,y,t, it substracts the spatial average
        imStack = self.image_stack
        dimensions = imStack.shape
        imAux = reshape(imStack, 
                        (prod(dimensions[0:-1]), dimensions[-1]), order="F")
        globalMean = nanmean(imAux, axis=0)
        megaMean = tile(globalMean, (prod(dimensions[0:-1]), 1));
        megaMean = reshape(megaMean, dimensions);
        self.image_stack = imStack - megaMean;
        print('Spatial average removed!')
    
    def get_time_traces(self):
 
        try:
            lol = self.stimulation_start_id
        except AttributeError:
            self.stimulation_start_id = 0
            
        try:
            lol = self.spot_size
        except AttributeError:
            self.spot_size = [5, 5] 
            
        try:
            lol = self.selected_point
        except AttributeError:
            self.selected_point = array([[self.spot_size[0], 
                                          self.spot_size[1]]]) 
        

        self.time = (arange(0, self.image_stack.shape[-1]) 
                    - self.stimulation_start_id)/self.sampling_freq

    
        self.time_traces = zeros(
            (self.selected_point.shape[0], self.image_stack.shape[-1]))



        for n_point in range(self.selected_point.shape[0]):
            self.time_traces[n_point, :] = get_roi_time_trace(
                                            self.image_stack, 
                                            self.spot_size, 
                                            self.selected_point[n_point, :])

        try:
            self.Cycle.time_traces = zeros(
                (self.selected_point.shape[0], self.Cycle.mean.shape[-1]))
            self.Cycle.sem_traces = zeros(
                (self.selected_point.shape[0], self.Cycle.mean.shape[-1]))
            for n_point in range(self.selected_point.shape[0]):
                self.Cycle.time_traces[n_point, :] = get_roi_time_trace(
                                            self.Cycle.mean, 
                                            self.spot_size, 
                                            self.selected_point[n_point, :])
                self.Cycle.sem_traces[n_point, :] = get_roi_time_trace(
                                            self.Cycle.sem, 
                                            self.spot_size, 
                                            self.selected_point[n_point, :])
        except AttributeError:
            pass    
        
        print('Time traces extracted!')
        
    def get_freq_traces(self):
        
        try:
            lol = self.spot_size
        except AttributeError:
            self.spot_size = [5, 5] 
            
        try:
            lol = self.selected_point
        except AttributeError:
            self.selected_point = array([[self.spot_size[0], 
                                          self.spot_size[1]]]) 
        
   
        self.freq_traces = zeros(
            (self.selected_point.shape[0], self.fft_stack.shape[-1]))
        for n_point in range(self.selected_point.shape[0]):
            self.freq_traces[n_point, :] = get_roi_time_trace(
                                            abs(self.fft_stack), 
                                            self.spot_size, 
                                            self.selected_point[n_point, :])
        print('Frequency traces extracted!')    

    def mean_cycle(
        self, 
        offset = 0, pulse_duration = 150e-3, prf = 1/30, 
        cycle_offset_id = 0):
        # Calculates the temporal mean, limit values and standard deviation for every 
        # pixel using the information of the stimulation.
        # Useful for statistical analysis and dealing with noisy data.
        # If specific cycles should be excluded, self.exclude_stim list should be provided.
        # Input:
        #   self.image_stack: 3D or 2D array where the last dimension will be averaged.
        #   offset: time in (s) before the first stimulation
        #   pulse_duration: (s) single stimulation duration
        #   prf: pulse repetition frequency (Hz)
        #   cycle_offset_id: Samples to shift the cycle for visualisation. 
        #                    Positive ==> Shift to the right  
        #   self.exclude_stim: list with the cycle number to be excluded from the mean.
        # Output: 
        #   Cycle: structure containing
        #        .time 
        #        .mean
        #        .max
        #        .min
        #        .std
        #
        # 01.04.2020 hector.estrada@posteo.org
        
        try:
            if (self.emission_offset != offset) and (offset != 0):
                self.emission_offset = offset # (s)
        except AttributeError:
            self.emission_offset = offset # (s)
            
        try:
            if ((self.emission_pulse_duration != pulse_duration) and 
                (pulse_duration != 150e-3)):
                self.emission_pulse_duration = pulse_duration # (s)
        except AttributeError:
            self.emission_pulse_duration = pulse_duration # (s)

        try:
            if (self.emission_rep_frequency != prf) and (prf !=  1/30):
                self.emission_rep_frequency = prf # (Hz)
        except AttributeError:
            self.emission_rep_frequency = prf # (Hz)
            
        try:
            if ((self.cycle_offset_id != cycle_offset_id) and 
                (cycle_offset_id != 0)):
                self.cycle_offset_id = cycle_offset_id # samples
        except AttributeError:
            self.cycle_offset_id = cycle_offset_id # samples
            
        
        # (1.0) Calculating the dimensions, cycle length, and reshaping
        self.dimensions = list(self.image_stack.shape)
        im_stack = reshape(self.image_stack, 
                            (prod(self.dimensions[0:-1]), self.dimensions[-1]),
                            order="F")
        self.time = arange(0, self.dimensions[-1])/self.sampling_freq \
                    - self.emission_offset # (s)
        xtra_time = self.cycle_offset_id/self.sampling_freq        
        n_emissions = int((self.time[-1] + xtra_time )
                          *self.emission_rep_frequency)
        print(n_emissions, ' cycles detected')
        print('self.sampling_freq', self.sampling_freq)
        print('self.emission_rep_frequency', self.emission_rep_frequency)
        cycle_length = int(round(self.sampling_freq/self.emission_rep_frequency))

        print('cycle_length', cycle_length)
        cycle_start = int(round(self.sampling_freq
                                              *self.emission_offset) 
                                              - self.cycle_offset_id)
            
        # (1.1) Cropping the input to the useful size
        time = self.time[cycle_start:cycle_start 
                                         + n_emissions*cycle_length]
        im_stack = im_stack[:, cycle_start:cycle_start
                                                    + n_emissions*cycle_length]
        
        self.dimensions[-1] = cycle_length

        # (1.2) Arrange the array size
        im_stack = reshape(im_stack, 
                           (prod(self.dimensions[:-1]), cycle_length, n_emissions),
                           order="F")

        # Are there any cycles to exlude?
        try:
          exclusion_flag = self.exclude_stim
        except AttributeError:
          exclusion_flag = []
        
        print('im_stack.shape before ' + repr(im_stack.shape))
        if len(exclusion_flag) != 0:
          im_stack = delete(im_stack, exclusion_flag, axis=2)
        
        print('im_stack.shape after'+ repr(im_stack.shape))

        # (2.0) Calculating the stats
        self.Cycle = MeanCycle()
        self.Cycle.n_emissions = n_emissions
        self.Cycle.mean = squeeze(nanmean(im_stack, axis=2))
        self.Cycle.max = squeeze(nanmax(im_stack, axis=2))
        self.Cycle.min = squeeze(nanmin(im_stack, axis=2))
        self.Cycle.std = squeeze(nanstd(im_stack, axis=2))
        self.Cycle.sem = squeeze(sem(im_stack, axis=2))
        self.Cycle.time = time[0:cycle_length]


            
        # (3.0) Reshaping the arrays 
        self.Cycle.mean = reshape(self.Cycle.mean, self.dimensions, order="F")
        self.Cycle.max = reshape(self.Cycle.max, self.dimensions, order="F")
        self.Cycle.min = reshape(self.Cycle.min, self.dimensions, order="F")
        self.Cycle.std = reshape(self.Cycle.std, self.dimensions, order="F")
        self.Cycle.sem = reshape(self.Cycle.sem, self.dimensions, order="F")




        # time trace from the averaged sequences
        self.time_traces_avg = zeros((self.selected_point.shape[0],
                self.Cycle.mean.shape[-1]))
        for n_point in range(self.selected_point.shape[0]):
            self.time_traces_avg[n_point,:] = get_roi_time_trace(
                                                self.Cycle.mean,
                                                self.spot_size,
                                                self.selected_point[n_point,:])

    
        
        print('Cycle average performed!')
    
    def temporal_smoothing(self, filt_param):
        ## Performs Savitsky Golay filtering along the time dimension
        # filt_param = array([window_length, polyorder])
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.savgol_filter.html
        # Reshape for speed
        
        A = self.Cycle.mean
        lol = A.shape
        A = reshape(A, (lol[0]*lol[1], lol[2]))
        A = savgol_filter(A, filt_param[0], filt_param[1], axis=1)
      
        self.Cycle.mean = reshape(A, (lol[0], lol[1], lol[2]))
        print('Temporal smoothing ready!')
    
    def define_mask(self, threshold = 0.3):
        try:
            if (self.mask_threshold != threshold) and (threshold != 0.3):
                self.mask_threshold = threshold # relative to maximum
        except AttributeError:
            self.mask_threshold = threshold # relative to maximum
        
        return self.ref_image > self.mask_threshold*nanmax(self.ref_image)
    
    def nan_mask(self):
        for n_frame in range(self.image_stack.shape[-1]):
            self.image_stack[:,:, n_frame] = self.image_stack[:,:, n_frame]* \
            self.brain_outline_mask
        self.image_stack[self.image_stack==0] = NaN
        print('NaN masking performed!')

    def show_time_traces(self, frame_id = 1, tlim = [0, -1]):
        # Plots a given frame together with several timetraces at its right side.
        # Input:
        #   time: 1D array containing the time for the time traces.
        #   image_id: Nx2 indexes of the images where the time traces have been selected.
        #   time_traces: 2D Array with the values at given time points.
        #   frame_id: Index of the frame movie_frame so that time[frame_id] yields the
        #             time when the frame takes place.
        #   movie_frame: 2D image at point frame_id.
        #   mask: 2d binary array to mask movie_frame.
        #   background: image to be shown in gray scale as background for the movie frame.
        #   rgba_colormap: rgba custom colormap.
        #   tlim: limits the x axis display to only the selected range. [0, -1]
        #         is the full range.
        # Output:
        #   fig: figure handles
        #   im2: movie_frame handles
        #   ln: line art handles of the line indicating the image time.
        # 05.05.2020 hector.estrada@posteo.org
        #### TODO: implement range visualization for Cycle.
        
        
        try:
            cycle_flag = 1
            lol = self.Cycle
        except AttributeError:
            cycle_flag = 0
        
        try:
            if (self.see_frame_id != frame_id) and (frame_id != 1):
                self.see_frame_id = frame_id # index
        except AttributeError:
            self.see_frame_id = frame_id # index
        
        asymm_flag = 0
        # (1.0) Prepare the figure
        pageSetup(300, 0.30)
        fig = figure()
        fig.set_tight_layout(True)
        
        # Prepare the mask
        self.mask = self.define_mask()*self.brain_outline_mask


        
        # (2.0) Plotting
        if cycle_flag:
            time_traces = self.Cycle.time_traces
            time = self.Cycle.time
            movie_frame = self.Cycle.mean[:, :, self.see_frame_id]*100
            Mip = amax(self.Cycle.mean, axis = 2)*self.mask
            mip = amin(self.Cycle.mean, axis = 2)*self.mask
            low = time_traces - self.Cycle.sem_traces
            high = time_traces + self.Cycle.sem_traces

        else:
            time_traces = self.time_traces
            time = self.time
            movie_frame = self.image_stack[:, :, self.see_frame_id]*100
            Mip = amax(self.image_stack, axis = 2)*self.mask
            mip = amin(self.image_stack, axis = 2)*self.mask
        
        mip = -robust_max(-mip, 0.9995)
        Mip = robust_max(Mip, 0.9995)
        
        if self.polarity in ('asy', 'asymmetric'):
            M_global = array([mip, Mip])*self.vis_threshold
            if (mip < 0):
                self.polarity = M_global
                asymm_flag = 1
            else:
                # Because sometimes you don't know this in advance
                self.polarity = 'monopolar'
                self.c_map_name  = 'inferno'
        else:
            aux = max([Mip, -mip])
            M_global = array([-aux, aux])*self.vis_threshold
            
        # Local maximum from the time traces
        MM = nanmax(abs(time_traces))*100
        #MM = nanmax(time_traces[0, :]) * 100
        
        # (2.1) Time traces
        ax0 = subplot(1, 2, 2)
        # Additional ticks on the emissions
        if not cycle_flag and self.stim_frames[0]:
            for nn in range(self.stim_frames.shape[0]):
                ax0.plot([time[self.stim_frames[nn]], 
                          time[self.stim_frames[nn]]], 
                          [MM, (2*self.selected_point.shape[0] - 1)*MM], '0.8')
       
        ## Range
        if cycle_flag:
            for n_point in range(self.selected_point.shape[0]):
                ax0.fill_between(time, low[n_point, :]*100 + 2*n_point*MM, 
                                 high[n_point, :]*100 + 2*n_point*MM,
                                 alpha = 0.2,
                                 color = rcParams['axes.prop_cycle'].by_key()['color'][n_point])
            
        # Time traces
        #print('tlim',len(tlim), tlim[0])

        for n_point in range(self.selected_point.shape[0]):
            ax0.plot(time, time_traces[n_point, :]*100 + 2*n_point*MM)
        
        ln, = ax0.plot([time[self.see_frame_id], time[self.see_frame_id]], 
                       [-MM, (2*self.selected_point.shape[0] - 1)*MM], 
                       color=(0.3, 0.3, 0.3))
        ax0.plot([time[tlim[0]], time[tlim[0]]], 
                      [0, round(MM)], 'k', linewidth = 4)
        
        
        
        # show the ranges in the plot
        #ax0.vlines(-10/20,-1,1, colors='r')
        #ax0.vlines(0/20,-1,1, colors='r')
        #ax0.vlines(4/20,-1,1, colors='r')
        #ax0.vlines(14/20,-1,1, colors='r')
        
        ax0.set_xlim(time[tlim[0]], time[tlim[1] - 1])
        yticks(ticks = 2*MM*arange(self.selected_point.shape[0]), labels = [])
        # Adding vertical scale
        time_range = time[tlim[1] - 1] - time[tlim[0]]
        text_x = time[tlim[0]] - time_range*0.035

        text(text_x, round(MM)*0.4, str(int(round(MM))))
        
        grid('on')
        xlabel('Time (s)')
        ylabel('$\Delta$F/F (%)')
        # Removing the box around the axes
        ax0.spines["top"].set_visible(False)
        ax0.spines["right"].set_visible(False)
        ax0.spines["bottom"].set_visible(False)
        ax0.spines["left"].set_visible(False)
        
        # (2.2) Overlay figure
        rgba_colormap = define_colormap(
            name = self.c_map_name, 
            transparency = self.overlay_transparency, 
            polarity = self.polarity)

        ax1 = subplot(1, 2, 1)
        im1 = ax1.imshow(self.ref_image*self.mask, cmap='gray', origin='upper')
        im2 = ax1.imshow(movie_frame*self.mask, cmap=rgba_colormap, 
                         vmin= M_global[0], vmax= M_global[1], 
                         origin='upper')
        if (len(self.brain_segment_mask.shape)>2):
            im3 = ax1.imshow(self.brain_segment_mask[:,:,0], 
                             alpha = self.brain_segment_mask[:,:,3],
                             cmap='Wistia')
        axis('off')
        for n_point in range(self.selected_point.shape[0]):
            scatter(self.selected_point[n_point, 1],
                    self.selected_point[n_point, 0],
                    s = 15,
                    facecolors='none',
                    edgecolors=rcParams['axes.prop_cycle'].by_key()['color'][n_point])

        divider = make_axes_locatable(ax1)
        cax1 = divider.append_axes("right", size="5%", pad=0.05)
        colorbar(im2, cax=cax1)

        if self.plot_visibility:
            show(block=False)
        if asymm_flag:
            self.polarity = 'asymmetric'
        return fig, im2, ln

       
    def get_metrics(self, time_ranges = array([0, 1])):
        """
            Calculates area under the curve, maximum and minimum of cycle-averaged time traces
            and activated area (sum of pixels above a certain defined threshold) within the 
            time range defined by time_ranges
            Input:
            time_ranges: array([[t0, t1], [t2 t3], ...,  [tN-1, tN]]) containing the time points 
                        (s) where the metrics should be extracted.
            self.F0: Fluorescence threshold to be considered as activation. If negative, the
                     value corresponds to the peak-normalized threshold.
            self.Rmax: Maximum radius for hard and soft threshold. Not considered in
                       flood fill calculations.
            Output:
                self.metrics: N points x M time regions x 4 (number of metrics)
                self.activation_area: N points (mm^2)
                self.activation_mip: N x image size. delta F/F in (%)
                
           
        """
        # Default pixel size if not defined
        try:
            self.dx = self.dy
        except AttributeError:
            self.dx = 1.
            self.dy = 1. # pixel size in mm
        # Default activation threshold for activated area
        try:
            mk = self.F0
        except AttributeError:
            self.F0 = 0.5 # in %
        # Default maximum radius for activated area
        try:
            mk = self.Rmax
        except AttributeError:
            self.Rmax = 3 # (mm)
                
        M = shape(time_ranges)[0]
        N = shape(self.selected_point)[0]
        self.metrics = zeros((N, M, 4))
        self.activated_area = zeros((N, M))
        self.activation_mip = []
        
        normal_flag = 0
        if self.F0 < 0:
            normal_flag = 1
            
        for m_time in range(M):
            t1 = argmin(abs(time_ranges[m_time, 0]-self.Cycle.time)) # First time point
            t2 = argmin(abs(time_ranges[m_time, 1]-self.Cycle.time)) # First time point
            peak_fluo = amax(self.Cycle.mean[:,:,t1:t2], axis=2)*100*self.brain_outline_mask # (%)
            # Normalized peak fluo
            if normal_flag:
                peak_fluo = peak_fluo/max(peak_fluo)
                
            # self.activation_mip.append(peak_fluo) # Full 2D MIP 
            # mask_activates = peak_fluo >= self.F0 # Absolute threshold in (%)
            for n_points in range(N):
                self.metrics[n_points, m_time, 0] = sum(self.Cycle.time_traces[n_points, t1:t2])
                self.metrics[n_points, m_time, 1] = max(self.Cycle.time_traces[n_points, t1:t2])
                self.metrics[n_points, m_time, 2] = min(self.Cycle.time_traces[n_points, t1:t2])
                tmax = argmax(self.Cycle.time_traces[n_points, t1:t2])
                self.metrics[n_points, m_time, 3] = self.Cycle.time[t1+tmax]
                # For area calculation
                #xx = arange(0,shape(self.image_stack)[0])*self.dx
                #yy = arange(0,shape(self.image_stack)[1])*self.dy
                x0 = int(self.selected_point[n_points, 0])
                y0 = int(self.selected_point[n_points, 1])
                #x, y = meshgrid(xx, yy, indexing = 'ij' )
                #radius = sqrt((x - xx[x0])**2 + (y - yy[y0])**2)
                #radial_mask = radius <= self.Rmax # Hard
                #radial_mask = exp(-0.5*radius**2/self.Rmax**2) # Soft. Rmax becomes the standard deviation
                #mask_activates = peak_fluo*radial_mask >= self.F0
                
                # Flood fill
                masked = peak_fluo
                if normal_flag:
                    masked[peak_fluo < -self.F0] = 0
                    masked[peak_fluo >= -self.F0] = 1
                else:
                    masked[peak_fluo < self.F0] = 0
                    masked[peak_fluo >= self.F0] = 1
                
                self.activated_area[n_points, m_time] = 0
                mask_activates = zeros(shape(peak_fluo))
                if masked[x0, y0]:
                    mask_activates = seg.flood(masked, (x0, y0), tolerance = 0.0)  # Experiment with tolerance
                    self.activated_area[n_points, m_time] = sum(mask_activates)*self.dx*self.dy
                self.activation_mip.append(mask_activates)
                
                #dx = (max_pixel_range[0] - 1)/2 # odd number required
                #dy = (max_pixel_range[1] - 1)/2
                #ROIx = arange(int(self.selected_point[n_points, 0] - dx), int(self.selected_point[n_points, 0] + dx + 1))
                #ROIy = arange(int(self.selected_point[n_points, 1] - dy), int(self.selected_point[n_points, 1] + dy + 1))
                #self.metrics[n_points, m_time, 3] = sum(mask_activates[ROIx, ROIy])*pixel_size


class LargeParameterList:

    def __init__(self):
        self.folder_in = []
        self.folder_out = []
        self.file_name = []
        self.pressure = []
        self.sub_folder = []
        self.mask_file = []
        self.us_start_id = []
        self.cycle_offset_id = []
        self.stim_period = []
        self.stim_duration = []
        self.time_crop = []
        self.midline_id = []
        self.hemisphere_sign = []
        self.exclude_stim = []
        self.points = []
        self.metrics = MeanCycle()
        self.traces = MeanCycle()
        self.general = MeanCycle()
        
    def get_one(self, N):
        lol = LargeParameterList()
        # Intermediaries
        folder_in = self.folder_in[N]
        folder_out = self.folder_out
        sub_folder = self.sub_folder[N]
        file_name = self.file_name[N]
        mask_file = self.mask_file[N]
        points = self.points[N]
        # Go direct
        lol.general = self.general
        lol.us_start_id = self.us_start_id[N]
        lol.cycle_offset_id = self.cycle_offset_id[N]
        lol.stim_period = self.stim_period[N]
        lol.stim_duration = self.stim_duration[N]
        lol.midline_id = self.midline_id[N]
        lol.hemisphere_sign = self.hemisphere_sign[N]
        lol.exclude_stim = self.exclude_stim[N]
        
        lol.time_crop = self.time_crop
        lol.pressure = self.pressure[N]
        # File paths
        lol.file_in = folder_in + sub_folder + file_name
        lol.figure_out = folder_out + sub_folder + file_name
        lol.brain_outline_file = folder_in + sub_folder + mask_file
        # Points with mirror symmetry
        n_points = shape(points)[0]
        lol.selected_point = zeros((2*n_points, 2))
        for nn in range(n_points):
            lol.selected_point[2*nn, :] = array([points[nn, 0], lol.midline_id - (points[nn, 1] - lol.midline_id)])
            lol.selected_point[2*nn + 1, :] = points[nn,:]
        
        return lol
    
    def load_param(self, this_file):
        try:
            with open(this_file, "rb") as f:
                return pickle.load(f)
        except Exception as ex:
            print("Error during unpickling object (Possibly unsupported):", ex)
            
    def load_csv(self, csv_file_name):
        try:
            gg = []
            with open(csv_file_name, "r") as csvfile:
                reader_variable = csv.reader(csvfile, delimiter="\t")
                for row in reader_variable:
                    gg.append(row)
            # Finding the structure
            self.folder_out = gg[1][0]
            for n_data in range(1, shape(gg)[0]):
                self.folder_in.append(gg[n_data][1])
                self.sub_folder.append(gg[n_data][2])
                self.file_name.append(gg[n_data][3])
                self.mask_file.append(gg[n_data][4])
                self.pressure.append(float(gg[n_data][5]))
                self.us_start_id.append(int(gg[n_data][6]))
                self.cycle_offset_id.append(int(gg[n_data][7]))
                self.stim_period.append(float(gg[n_data][8]))
                self.stim_duration.append(float(gg[n_data][9]))
                self.midline_id.append(int(gg[n_data][10]))
                self.hemisphere_sign.append(int(gg[n_data][11]))
                self.exclude_stim.append(eval(gg[n_data][12]))
                self.points.append(eval(gg[n_data][13]))
                    
        except Exception as ex:
            print("Error during importing object (Possibly unsupported):", ex)
        
    def save_param(self, this_file):
        try:
            with open(this_file, "wb") as f:
                pickle.dump(self, f, protocol=pickle.HIGHEST_PROTOCOL)
                print('File ' + this_file + ' saved succesfully')
        except Exception as ex:
            print("Error during pickling object (Possibly unsupported):", ex)
        
def get_roi_time_trace(im_stack, spot_size = [5, 5], selected_point = [1, 1]):
    # function get_roi_time_trace(dF_F, param)
    # Calculates the average on a given region defined by param for all the 
    # different filtered images from dF_F.
    # Input:
    #      spot_size: [1 x 2] in samples. It must be an odd number.
    #      selected_point: [1 x 2] in samples
    #      im_stack: cell array        
    # Output:
    #      time_trace: array with the time traces
    # hector.estrada@posteo.org
      
    # (1.0) Defining the ROI
    dx = (spot_size[0] - 1)/2 # odd number required
    dy = (spot_size[1] - 1)/2
    ROIx = [int(selected_point[0] - dx), int(selected_point[0] + dx + 1)]
    ROIy = [int(selected_point[1] - dy), int(selected_point[1] + dy + 1)]

    time_trace = mean_roi(im_stack[ROIx[0]:ROIx[1], ROIy[0]:ROIy[1], :])
    return time_trace

def mean_roi(input_array):
    # Function to calculate the mean along the first array dimensions of the 
    # input_array
    # Imput:
    #       input_array: Set of images arranged in an array as 
    # Output: 
    #       vector_mean: vector containing the mean 
    # hector.estrada@posteo.org

    # (1.0) Arranging the array dimensions
    dimensions = input_array.shape

    # (2.0) Reshaping and allocating 
    Y = reshape(input_array, (prod(dimensions[0:-1]), dimensions[-1]), order="F") # Quasi flattened 
    vector_mean = zeros((1, dimensions[-1]))

    # (3.0) Calculation
    vector_mean = nanmean(Y, axis=0);
    return vector_mean

def define_colormap(name = 'PRGn', points = 255, transparency = 'linear', 
                    polarity = 'bipolar'):
    
    # Constant for logistic
    lamb = 15
    #if polarity in (['mono', 'monopolar', 'bi', 'bipolar']):
    if isinstance(polarity, str):
        if transparency in ('log', 'logistic'):
            if polarity in ('mono', 'monopolar'):
                alph = abs(2*(1./(1 + exp(-lamb*linspace(0, 1, points))) - 0.5))
            else:
                alph = abs(2*(1./(1 + exp(-lamb*linspace(-1, 1, points))) - 0.5))
        if transparency in ('cos', 'cosine'):
            if polarity in ('mono', 'monopolar'):
                alph = 0.5*(cos(linspace(0, 1, points)*pi + pi) + 1)
            else:
                alph = 0.5*(cos(linspace(0, 1, points)*2*pi) + 1)
        if transparency in ('tukey', 'tuk'):
            if polarity in ('bi', 'bipolar'):
                a = 0.75 #0.55 #0.65 #0.75
                M = round(0.4*points) #0.65 # 0.5 #0.4
                alph = ones(points)
                alph[round((points - M)/2):round((points - M)/2)+M] = 1 - tukeywin(M, a) 
            else:
                alph = 0 #### Not implemented yet
                print('Not implemented yet!')
        if transparency in ('lin', 'linear'):
            if polarity in ('mono', 'monopolar'):
                alph = abs(linspace(0, 1, points)) # Linear transparency
            else:
                alph = abs(linspace(-1, 1, points)) # Linear transparency
        ## For debugging purposes
        #figure(655)
        #plot(alph)
        #show(block=False)
        #lol = input("Check the plot, if it is fine, press ENTER")        
        if transparency in ('flat'):
            alph = ones((points))    
        else:        
            alph = (alph - nanmin(alph))/(nanmax(alph) - nanmin(alph))
        return rgbaCmap(name, alph, 0)
    else:
    #if polarity in ('asy', 'asymmetric'):
        # Obtain color edges
        alph = array([1, 0, 1])
        jojo = rgbaCmap(name, alph, 0)
        vrgba = zeros((3, 5))
        norm_value = array([polarity[0], 0, polarity[1]])
        norm_value = (norm_value - min(norm_value)) \
                    /(max(norm_value) - min(norm_value))
        vrgba[:, 0] = norm_value
        vrgba[:, 1:] = jojo.colors
        new_cmap = new_rgba_colormap(vrgba, points, 0)
        
        if transparency in ('cos', 'cosine'):
            nl_alpha = 0.5*(cos(new_cmap.colors[:, 3]*pi + pi) + 1)
            nl_alpha = (nl_alpha - nanmin(nl_alpha))/(nanmax(nl_alpha) \
                                 - nanmin(nl_alpha))
            # Debugging
            #figure(655)
            #plot(new_cmap.colors[:, 3], nl_alpha)
            #xlabel('Linear')
            #ylabel('cosine')
            #show(block=False)
            #lol = input("Check the plot, if it is fine, press ENTER")
            #close()
            new_cmap.colors[:, 3] = nl_alpha
            
        if transparency in ('log', 'logistic'):
            nl_alpha = 2*(1./(1 + exp(-lamb*new_cmap.colors[:, 3])) - 0.5)
            nl_alpha = (nl_alpha - nanmin(nl_alpha))/(nanmax(nl_alpha) \
                                 - nanmin(nl_alpha))
            # Debugging
            #figure(655)
            #plot(new_cmap.colors[:, 3], nl_alpha)
            #xlabel('Linear')
            #ylabel('Logistic')
            #show(block=False)
            #lol = input("Check the plot, if it is fine, press ENTER")
            #close()
            new_cmap.colors[:, 3] = nl_alpha
            
        if transparency in ('flat'):
            new_cmap.colors[:, 3] = ones(new_cmap.colors.shape[0])
            
        new_cmap._init()
        return new_cmap

def robust_max(x, threshold):
    # Function to calculate the maximum in an array without considering outliers 
    # using the cumulative sum of the histogram
    # Imput:
    #       x: array of arbitrary dimensions 
    # Output: 
    #       M: maximum 
    # hector.estrada@posteo.org
    
    lol = histogram(x.flatten(), bins = 100)
    nn = cumsum(lol[0])
    nn = (nn-nn.min())/(nn.max() - nn.min())
    
    return lol[1][argmax(nn > threshold)]
