% Global script to run the preprocessing over a whole folder
addpath(genpath(pwd))

% To run this script you need to download two datasets shown in Figs. 2A,2B, and Fig. 3B, 3C from url 

file_name = {'single_focus_us_stim';
            'multiple_foci_us_stim'};
            
% Edit the following folders to match your directory
folder_name_tif = '../data/raw/';
folder_name_mat = '../output/';
save_folder = '../output/';

StartTimings = [499; 526];% The stimulation starts on this frames.

for n = 1:length(file_name)
    param.File.FullName = [folder_name_tif, file_name{n}, '_kalman.tif'];
    param.File.FullNameMat = [folder_name_mat, file_name{n}, '.mat'];
    param.File.SaveFolder = save_folder; %  Name of the folder to be created to save the data
    param.File.SaveFile = [file_name{n}, '_kalman_pretreated.mat']; % Name of the file where to save the data

    param.stimNb        = 20; % total nb of stimulations. % 30(sequential)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    param.Movie.FileName = '';
    param.Movie.jumpFrames = 2;
    param.Data.SamplingFreq = 20; % (Hz)
    param.SpatialFilter.Fourier.freqBands = [];
    param.checkROIFlag = 0;
    
    %% Neda_dips
    param.Crop.LastUsefulFrame = [4600];
    param.Crop.MaskThreshold = [0];
    param.Crop.xy = [7  211    9  217];
    param.Rotate.MidLine = [152.0524      67.98429;
                            41.53403      141.2984];
    param.Rotate.Cropxy = [45  254   30  272];
    param.Visualize.SelectedPoint = [113   93];

    %% Yiming_230404
    param.Crop.BaseLineRange = [26  118];
    % param.Crop.LastUsefulFrame = [4551];
    % param.Crop.MaskThreshold = [0];
    % param.Crop.xy = [6  186   19  201];
    % param.Rotate.MidLine = [139.9519      67.51337;
    %                         29.36898      141.8877];
    % param.Rotate.Cropxy = [40  212   33  215];
    % param.Visualize.SelectedPoint = [81  65];

    param.Filter.Mode = 'smooth';
    param.Filter.Smooth.Order = [2];
    param.Filter.Smooth.Length = [11];

    param.RelativeChange.Method{1, 1} = 'MovingBaseline';
    param.RelativeChange.Method{1, 2} = [1, 0]; % Fourier filtered array id according to freqBa
    param.MovBaseline.WindowSize = 500; % Sample length
    param.MovBaseline.Percentile = 10; % (%)

    param.Filter2D.GaussianFilterSize = []; % Empty turns it off
    param.Filter2D.MedianFilterSize = []; % Empty turns it off

    param.Visualize.SpotSize = [5, 5]; % odd number of pixels
    param.Visualize.SelectedFrame = 2000;
    param.Visualize.CropLastFrames = 1;
    param.Visualize.SeeFigure = 0;
    param.Visualize.EmissionSequence = [];
    param.Visualize.checkAntiDivergence = 0;
    param.Visualize.KeepBaseline = 1;
    param.Visualize.CheckMovingBaseline = prod(param.Visualize.SelectedPoint);

    param.StimulationTimings = StartTimings(n) + [0:200:19*200].';
    %%=========================================================================
%     try
        process_image_stack(param);
%     catch error
%         disp(['process_image_stack error acq ', file_name{n}]);
%     end

    clear param
end


