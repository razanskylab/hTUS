function process_image_stack(param)
% Performs basic filtering, crop, and rotation to image stacks acquired using
% fluorescence. It returns basic png sequences and plots of the important
% image areas. I can also save postprocessed data and the used parameters in a
% mat file.
% Input:
% param.File.FullName: Full file name including path
% param.File.SaveFolder: Name of the folder to be created to save the data
% param.File.SaveFile: Name of the file where to save the data
%
% param.Movie.FileName: Full path of the png sequence to build the movie upon.
%                       If empty, no movie is generated or calculated.
%
% param.Data.SamplingFreq: Frequency at which the image sequence was acquired. (Hz)
%
% param.Crop.BaseLineRange: Range index where no excitation happens and taken as
%                           reference for relative fluorescence changes.
% param.Crop.LastUsefulFrame: Last freane to be used in the processing.
% param.Crop.MaskThreshold: Relative cropping the image's dynamic range (0,1).
% param.Crop.xy: 1x4 First crop on the raw image for faster processing
%                [Idx1, Idx2, Idy1, Idy2]
%
% param.Rotate.MidLine:  2 x 2 Midline coordinates to rotate the image
%                        [Idx1, Idy1; Idx2, Idy2]
% param.Rotate.Cropxy : 1 x 4 Cropping the image after rotation using Id
%                       [Idx1, Idx2, Idy1, Idy2]
%
% param.Filter.Mode : 'fourier' or 'smooth'. It uses either frequency domain
%                     band filtering or Savitzky-Golay smoothing.
%
% param.Filter.Smooth.Order
% param.Filter.Smooth.Length
%
% param.Filter.Fourier.freqBand: Check FD_filtering.m for details
% param.Filter.Fourier.filter: Check FD_filtering.m for details
% param.Filter.Fourier.snr: Check FD_filtering.m for details
%
% param.Filter2D.GaussianFilterSize: For 2D image processing. Empty turns it off
% param.Filter2D.MedianFilterSize: For 2D image processing. Empty turns it off
%
% param.Visualize.SelectedPoint: ID image coordinates
% param.Visualize.SpotSize:  Odd number of pixels
% param.Visualize.SelectedFrame: Frame selected for visualisation
%
% preBaseLine calculated the basic relative change first and performs filtering
% afterwards.
%
% Output:
%   If param.File.SaveFile is not empty, the postprocessed data is saved as follows
%     PostPro.dF_F = dF_F;
%     PostPro.slowdF_F = slowdF_F;
%     PostPro.fastdF_F = fastdF_F;
%     PostPro.ROIdF_F = timetrace0;
%     PostPro.ROIslowdF_F = timetrace1;
%     PostPro.ROIfastdF_F = timetrace2;
%
% hector.estrada@posteo.org

% (1.0) Load tiff file
[folderName, fileName, ext] = fileparts(param.File.FullName);
try
  infoImage = imfinfo(param.File.FullName);
  mImage = infoImage(1).Width;
  nImage = infoImage(1).Height;
  numberImages = param.Crop.LastUsefulFrame - param.Crop.BaseLineRange(1);% before the baseline are removed
  % Constrain for fourier-based analysis
  if mod(numberImages, 2) ~= 0
    numberImages = numberImages - 1;
  end

  im = zeros(nImage, mImage, numberImages, 'uint16');
  for i = 1:numberImages
      im(:,:,i) = imread(param.File.FullName, 'Index', ...
                        i + param.Crop.BaseLineRange(1) - 1, 'Info', infoImage);
  end
  im = double(im);
catch error
  warning("Error during tif file read. " + error.message);
end
disp('Image stack loaded')

% (1.1) Spatial filtering to smooth images
if ~isempty(param.Filter2D.GaussianFilterSize)
  for i = 1:numberImages
    aux = squeeze(im(:,:,i));
    im(:,:,i) = imgaussfilt(aux, param.Filter2D.GaussianFilterSize);
  end
end
if ~isempty(param.Filter2D.MedianFilterSize)
  for i = 1:numberImages
    aux = squeeze(im(:,:,i));
    im(:,:,i) = medfilt2(aux, param.Filter2D.MedianFilterSize);
  end
end

% (1.2) Cropping the image further to make the following steps faster
if ~isempty(param.Crop.xy)
  im = im(param.Crop.xy(1):param.Crop.xy(2), ...
          param.Crop.xy(3):param.Crop.xy(4), :);
end
tstart = tic;

% (1.3) Motion correction
% TFORM = IMREGTFORM(MOVING, FIXED, TRANSFORMTYPE, OPTIMIZER, METRIC)
% TFORM = IMREGCORR(MOVING, FIXED, TRANSFORMTYPE)
% % Example
%     -------
%     % This example solves a registration problem in which the cameraman 
%     % image is synthetically scaled and rotated.
%   
%       fixed  = imread('cameraman.tif');
%       theta = 20;
%       S = 2.3;
%       tform = affine2d([S.*cosd(theta) -S.*sind(theta) 0; 
%                         S.*sind(theta)  S.*cosd(theta) 0; 
%                         0               0              1]);
%       moving = imwarp(fixed,tform);
%       moving = moving + uint8(10*rand(size(moving)));
%   
%       tformEstimate = imregcorr(moving,fixed);
%   
%       figure, imshowpair(fixed,moving,'montage');
%   
%       % Apply estimated geometric transform to moving. Specify 'OutputView' to
%       % get registered moving image that is the same size as the fixed image.
%       Rfixed = imref2d(size(fixed));
%       movingReg = imwarp(moving,tformEstimate,'OutputView',Rfixed);
%   
%       figure, imshowpair(fixed,movingReg,'montage');
%       figure, imshowpair(fixed,movingReg,'falsecolor');


% (1.4) Background removal
if ~isempty(param.SpatialFilter.Fourier.freqBands)
  param.SpatialFilter.Fourier.samplingFreq = 2;
  % (1.1.1) Quasi-flattening the image stack for processing
  dimensions = size(im);
  imAux = reshape(im, [prod(dimensions(1:end - 1)), dimensions(end)]);
  % (1.1.2) Even number of samples for Fourier transform
  last = prod(dimensions(1:end - 1));
  if mod(prod(dimensions(1:end - 1)), 2) ~= 0
    last = prod(dimensions(1:end - 1)) - 1;
  end
  % (1.1.3) Spatial filtering
  [imFilt] = FD_filtering(imAux(1:last, :).', param.SpatialFilter.Fourier);
  % (1.1.4) back to the original array
  imAux(1:last, :) = imFilt.';
  im = reshape(imAux, dimensions);
  clear dimensions imAux imFilt
end

% (1.5) Quasi-flattening the image stack for processing
nFrames = size(im, 3);
dimensions = size(im);
im = reshape(im, [prod(dimensions(1:end - 1)), dimensions(end)]);
tstart = tic;

% (1.6) Previous scaling for preBaseline mode
nMethods = size(param.RelativeChange.Method, 1);
baseLineMean = mean(im(:, 1:diff(param.Crop.BaseLineRange) + 1), 2);
pre_baseline_flag = 0;
for nM = 1:nMethods
  switch param.RelativeChange.Method{nM, 1}
  case {'preBaseLine'}
    pre_baseline_flag = 1; % To avoid normalizing multiple times    
  end
end

if pre_baseline_flag
    norm_factor = repmat(baseLineMean, [1, nFrames]);
    im = im./norm_factor; % Calculating dF/F
end

switch param.Filter.Mode
  case {'fourier'}
    % (2.0) Fourier filtering
    param.Filter.Fourier.samplingFreq = param.Data.SamplingFreq;
    imFiltered = cell(size(param.Filter.Fourier.freqBands, 1), 1);
    [imFiltered{:}] = FD_filtering(im, param.Filter.Fourier);
  case{'smooth'}
    % (2.0) Smoothing filters for slow processes
    nBands = length(param.Filter.Smooth.Length);
    imFiltered = cell(nBands, 1);
    for nn = 1:nBands
      imFiltered{nn} = single(sgolayfilt(im, ...
                              param.Filter.Smooth.Order(nn), ...
                              param.Filter.Smooth.Length(nn), [], 2));
    end
  otherwise
    error('There is no such option for param.Filter.Mode');
end
disp('Image stack filtered')
telapsed = toc(tstart)

% (2.2) Mask calculation
%  permanent mask to avoid signal clipping.
% Tested using test/test_PCmask 07.07.20
[clippingRows, ~] = find(im <= 0);
clippingRows = unique(clippingRows);
mask0 = ones(size(im));
mask0(clippingRows, :) = NaN;


if ~isempty(param.Crop.MaskThreshold)
  mask = baseLineMean > param.Crop.MaskThreshold*max(baseLineMean(:));
end

% (3.0) Calculating dF/F
if param.Visualize.KeepBaseline
  param.Crop.BaseLineRange = -param.Crop.BaseLineRange;
end
dF_F = cell(nMethods, 1);

tstart = tic;
for nM = 1:nMethods
  switch param.RelativeChange.Method{nM, 1}
  case {'preBaseLine'}
    iD = param.RelativeChange.Method{nM, 2};
    dF_F{nM} = imFiltered{iD};
    disp([num2str(nM), ' is processed using ', ...
                                        param.RelativeChange.Method{nM, 1}])
  case {'BaseLine'}
    iD = param.RelativeChange.Method{nM, 2};
    dF_F{nM} = relative_change(imFiltered{iD}, param.Crop.BaseLineRange, ...
                          0, param.RelativeChange.Method{nM, 3});
    disp([num2str(nM), ' is processed using ', ...
                                        param.RelativeChange.Method{nM, 1}])
    case {'BaseLineFit'}
      error('This method is still too slow; it might take a 24 hours');
      megaMask = repmat(mask, [1, 1 , size(imFiltered{3}, 3)]); % Because we don't want to fit non sense
      iD = param.RelativeChange.Method{nM, 2}(1);
      idRef = param.RelativeChange.Method{nM, 2}(2);
      baseLineFit = fit_extrapolate(imFiltered{idRef}.*megaMask, ...
                                    param.Crop.BaseLineRange, ...
                                          param.RelativeChange.Method{nM, 3});
      dF_F{nM} = relative_change(imFiltered{iD}, param.Crop.BaseLineRange, ...
                                                              0, baseLineFit);
      disp([num2str(nM), ' is processed using ', ...
                                          param.RelativeChange.Method{nM, 1}])
    case {'MovingBaseline'}
      iD = param.RelativeChange.Method{nM, 2}(1);
      idRef = param.RelativeChange.Method{nM, 2}(2);
      if (idRef == 0) 
        % Using percentiles to calculate the baseline
        wSize = param.MovBaseline.WindowSize;
        nChuncks = floor(size(im, 2)/wSize);
        baseLineSeed = zeros(size(im, 1), 2*nChuncks - 1);
        middleTon = zeros(1, 2*nChuncks - 1);
        for nn = 1:2*nChuncks - 1
            first = (nn - 1)*wSize/2  + 1;
            last = first + wSize;
            middleTon(nn) = 0.5*(first + last);
            baseLineSeed(:, nn) = prctile(im(:, first:last), ...
                       param.MovBaseline.Percentile, 2, "Method","approximate");
        end
        movBaseline = interp1(middleTon, baseLineSeed.', [1:size(im, 2)],'pchip').';
        size(movBaseline)
      else
        % Using filtered data to calculate the baseline
        movBaseline = imFiltered{idRef};
      end

      % For concept testing
%        if (idRef == 0) 
%          % Using percentiles to calculate the baseline
%          wSize = 500;
%          nChuncks = floor(size(im, 2)/wSize);
%  %          movBaseline = zeros(1, size(im, 2));
%          baseLineSeed = zeros(1, nChuncks);
%          middleTon = zeros(1, nChuncks);
%          for nn = 1:2*nChuncks
%              first = (nn - 1)*wSize/2  + 1;
%              last = first + wSize;
%              middleTon(nn) = 0.5*(first + last);
%  %              if last > size(im, 2)
%  %                  last = size(im, 2);
%  %              end    
%              baseLineSeed(nn) = prctile(im(param.Visualize.CheckMovingBaseline, first:last), 10, 2, "Method","approximate");
%          end
%          movBaseline = interp1(middleTon,baseLineSeed,[1:size(im, 2)],'pchip');
%        else
%          % Using filtered data to calculate the baseline
%          movBaseline = imFiltered{idRef};
%        end
%        figure()
%        plot(middleTon, baseLineSeed, 'o')
%        hold on
%        plot(movBaseline)

      if (iD ~= 0)
        dF_F{nM} = relative_change(imFiltered{iD}, param.Crop.BaseLineRange, ...
                                                param.Crop.MaskThreshold, ...
                                                movBaseline);
      else % no filtering, moving baseline applied directly on raw data
        dF_F{nM} = relative_change(im, param.Crop.BaseLineRange, ...
                                                param.Crop.MaskThreshold, ...
                                                movBaseline);      
      end

      disp([num2str(nM), ' is processed using ', ...
                                          param.RelativeChange.Method{nM, 1}])
      
      % For educational and debugging purposes
      if param.Visualize.CheckMovingBaseline
        raw = im(param.Visualize.CheckMovingBaseline, :);
        if (iD ~= 0)
          signal = imFiltered{iD}(param.Visualize.CheckMovingBaseline, :);
        else
          signal = raw;
        end
        reference = movBaseline(param.Visualize.CheckMovingBaseline, :);
        
        jj = figure()
        plot(raw, 'k', 'LineWidth', 2)
        hold on
        plot(signal, 'b', 'LineWidth', 1.0)
        plot(reference, 'r', 'LineWidth', 1.5)

        xlabel('Time (Samples)')
        ylabel('Voxel intensity')
        legend('Raw', 'Filtered', 'Baseline')
        figOut = [param.File.SaveFolder,'/',param.File.SaveFile(1:end-4), '.fig'];
        savefig(jj, figOut)
      end
    case {'DontDiverge'}
      % (8.0) Fine tune anti-divergence parameter
      if isfield(param.Visualize,'checkAntiDivergence') 
        if param.Visualize.checkAntiDivergence
          param = check_filter_magnitude(im, param); % Updates anti-divergence parameter 
        end
      end
      iD = param.RelativeChange.Method{nM, 2};
      dF_F{nM} = relative_change(imFiltered{iD}, param.Crop.BaseLineRange, ...
                            -param.RelativeChange.Method{nM, 4}, ...
                             param.RelativeChange.Method{nM, 3});
      disp([num2str(nM), ' is processed using ', ...
                                          param.RelativeChange.Method{nM, 1}])
    case {'MoveDontDiverge'}
      % (8.0) Fine tune anti-divergence parameter
      if isfield(param.Visualize,'checkAntiDivergence') 
        if param.Visualize.checkAntiDivergence
          param = check_filter_magnitude(im, param); % Updates anti-divergence parameter 
        end
      end
      iD = param.RelativeChange.Method{nM, 5}(1);
      idRef = param.RelativeChange.Method{nM, 5}(2);
      dF_F{nM} = relative_change(imFiltered{iD}, param.Crop.BaseLineRange, ...
                            -param.RelativeChange.Method{nM, 4}, ...
                             imFiltered{idRef});
      disp([num2str(nM), ' is processed using ', ...
                                          param.RelativeChange.Method{nM, 1}])                                      
    otherwise
      error('This method is not implemented yet');
  end
end
disp('Image stack referenced')
telapsed = toc(tstart)

% (3.1) Masking NaN and zero values
% Masking for clipping voxels
if param.Visualize.KeepBaseline
  for nM = 1:nMethods
    dF_F{nM} = dF_F{nM}.*mask0;
  end
end
  
  
if ~isempty(param.Crop.MaskThreshold)
  mask1 = ones(size(mask));  
  mask1(mask==0) = NaN;  
  megaMask = repmat(mask1, [1, size(dF_F{1}, 2)]);
  for nM = 1:nMethods
    dF_F{nM} = dF_F{nM}.*megaMask;
  end
end

% (3.2) Final reshape and cropping out last useless frames
for nM = 1:nMethods
  new_dimensions = [dimensions(1:end - 1), size(dF_F{nM}, 2)...
                                              - param.Visualize.CropLastFrames];
  baseLineMean = reshape(baseLineMean, dimensions(1:end - 1));                                            
  dF_F{nM} = reshape(dF_F{nM}(:, 1:end - param.Visualize.CropLastFrames),...
                                                                new_dimensions);
end

% (4.0) Rotation and cropping
tstart = tic;
if ~isempty(param.Rotate.MidLine)
  % (4.1) Calculate the angle based on the brain's midline
  dY = (param.Rotate.MidLine(2,2) - param.Rotate.MidLine(1,2));
  dX = (param.Rotate.MidLine(2,1) - param.Rotate.MidLine(1,1));

  rotAngle = atan2(dY, dX)*180/pi + 90;

  if ~isempty(param.Rotate.Cropxy)
    % Rotating and croping the baseline
    auxdF_F = baseLineMean;
    newSize = [param.Rotate.Cropxy(2) - param.Rotate.Cropxy(1) + 1, ...
               param.Rotate.Cropxy(4) - param.Rotate.Cropxy(3) + 1];
    baseLineMean = zeros(newSize);
    baseLineMean = rotate_crop(auxdF_F, rotAngle, ...
                                               'bicubic', param.Rotate.Cropxy);
    clear auxdF_F
    for nM = 1:nMethods
      auxdF_F = dF_F{nM};
      newSize = [param.Rotate.Cropxy(2) - param.Rotate.Cropxy(1) + 1,...
                 param.Rotate.Cropxy(4) - param.Rotate.Cropxy(3) + 1,...
                size(auxdF_F, 3)];
       dF_F{nM} = zeros(newSize);
      for nn = 1:newSize(3)
        % (4.3) Rotate and crop unfiltered fluorescence changes
        dF_F{nM}(:, :, nn) = rotate_crop(squeeze(auxdF_F(:,:,nn)), rotAngle, ...
                                               'bicubic', param.Rotate.Cropxy);
      end
      clear auxdF_F
    end
  end
end
disp('Image stack rotated and cropped again')
telapsed = toc(tstart)

if param.checkROIFlag
  param.Visualize.SelectedPoint = check_roi(dF_F, param);
end
param.Visualize.SelectedPoint

% (5.0) Calculation on selected points performed after rotation
dx = (param.Visualize.SpotSize(1) - 1)/2;
dy = (param.Visualize.SpotSize(2) - 1)/2;
ROIx = param.Visualize.SelectedPoint(1) - dx:param.Visualize.SelectedPoint(1) + dx;
ROIy = param.Visualize.SelectedPoint(2) - dy:param.Visualize.SelectedPoint(2) + dy;

timeTrace = cell(nMethods, 1);
for nM = 1:nMethods
  timeTrace{nM} = mean_roi(dF_F{nM}(ROIx, ROIy, :));
end

% (6.0) Saving the postprocessed data
% (6.1) Creating the folder
folderFlag = ~isempty(param.File.SaveFolder) & ...
 (~isempty(param.File.SaveFile) || ~isempty(param.Movie.FileName));

if folderFlag
  try
    [status, msg, msgID] = mkdir(param.File.SaveFolder);
    disp(msg)
  catch error
    warning("Error creating the folder. " + error.message);
  end
end

% (6.2) Grouping and saving the data
% Accomodating for python readability
PostPro.BaseLine = baseLineMean;
PostPro.dF_F = dF_F;
PostPro.ROIdF_F = timeTrace;
%  if (nMethods == 1)
%    PostPro.dF_F = dF_F{1};
%    PostPro.ROIdF_F = timeTrace{1};
%  else
%    PostPro.dF_F = dF_F;
%    PostPro.ROIdF_F = timeTrace;
%  end
%  % Cropping out last 50 useless frames
%  if (param.Visualize.CropLastFrames ~= 0)
%    for nM = 1:nMethods
%      PostPro.dF_F{nM} = PostPro.dF_F{nM}(:,:,1:end - param.Visualize.CropLastFrames);
%      PostPro.ROIdF_F{nM} = PostPro.ROIdF_F{nM}(1:end - param.Visualize.CropLastFrames);
%    end
%  end
  
if ~isempty(param.File.SaveFile)
  fileOut = [param.File.SaveFolder,'/',param.File.SaveFile];
  save(fileOut, '-v7.3', 'PostPro', 'param');
  disp(['File ',fileOut,' saved!'])
end

% (7.0) Making a movie
if ~isempty(param.Movie.FileName)
  make_fluo_movie(PostPro, param)
end

disp('Finished!')

% (4.0) Checking the calculation for debugging purposes
% figure()
% aux = squeeze(mean(dF_F, 3));
% MM = max(abs(aux(:)));
% imagesc(aux);
% colormap(gray(129));
% colorbar
% title('Selected frame original')
%
% figure()
% aux = squeeze(dF_F(:,:, param.Visualize.SelectedFrame));
% MM = max(abs(aux(:)));
% imagesc(aux, [-MM, MM]);
% colormap(seismic(129));
% colorbar
% title('Selected frame, relative fluorescence change')
%
% figure()
% subplot(2, 2, 1)
% aux = squeeze(slowdF_F(:,:, param.Visualize.SelectedFrame));
% MM = max(abs(aux(:)));
% imagesc(aux, [-MM, MM]);
% colormap(seismic(129));
% colorbar
% title('dF/F. Selected frame')
%
% subplot(2, 2, 2)
% aux = squeeze(fastdF_F(:,:, param.Visualize.SelectedFrame));
% MM = max(abs(aux(:)));
% imagesc(aux, [-200, 200]);
% colormap(seismic(129));
% colorbar
% title('dF/F smooth. Selected frame')
%
% subplot(2, 2, 4)
% yyaxis left
% plot(timeTrace0)
% hold on
% plot(timeTrace1)
% yyaxis right
% plot(timeTrace2)
% title('dF/F in ROI')

% % (5.0) Plotting maximum at time
% NanMask = double(mask);
% NanMask(NanMask == 0) = NaN;
%
% % (5.1) Relative change
% [maxdF_F, maxFrame] = max(fastdF_F, [], 3);
% [mindF_F, minFrame] = min(fastdF_F, [], 3);
%
% maxdF_F = maxdF_F.*NanMask;
% maxFrame = maxFrame.*NanMask;
% [dFMaxFrame, dFMaxFrameCB] = color_coded_feature(maxdF_F, maxFrame, viridis(129), 2.0);
%
% mindF_F = mindF_F.*NanMask;
% minFrame = minFrame.*NanMask;
% [dFMinFrame, dFMinFrameCB] = color_coded_feature(mindF_F, minFrame, viridis(129), 2.0);
%
% figure()
% subplot(1, 2, 1)
% imagesc(dFMaxFrame)
% colormap(dFMaxFrameCB)
% axis image
% colorbar
% title('Time at maximum dF/F fast')
%
% subplot(1, 2, 2)
% imagesc(dFMinFrame)
% colormap(dFMinFrameCB)
% axis image
% colorbar
% title('Time at minima, dF/F fast')
%
% clear maxFrame, minFrame;
%
% % (5.1) Absolute change
% [maxIm, maxFrame] = max(slowdF_F, [], 3);
% [minIm, minFrame] = min(slowdF_F, [], 3);
%
% maxIm = maxIm.*NanMask;
% maxFrame = maxFrame.*NanMask;
% [imMaxFrame, imMaxFrameCB] = color_coded_feature(maxIm, maxFrame, viridis(129), 2.0);
%
% minIm = minIm.*NanMask;
% minFrame = minFrame.*NanMask;
% [imMinFrame, imMinFrameCB] = color_coded_feature(minIm, minFrame, viridis(129), 2.0);
%
% figure()
% subplot(1, 2, 1)
% imagesc(imMaxFrame)
% colormap(imMaxFrameCB)
% axis image
% colorbar
% title('Time at maxima, raw slow')
%
% subplot(1, 2, 2)
% imagesc(imMinFrame)
% colormap(imMinFrameCB)
% axis image
% colorbar
% title('Time at minima, raw slow')
