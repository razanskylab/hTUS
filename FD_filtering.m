function [varargout] = FD_filtering(imStack, fltParam)
% function [slow, fast] = FD_filtering(imStack, lowFreqBand, highFreqBand, samplingFreq);
% Function to calculate relative changes on an image stack based on a predefined
% range. It also returns a binary mask that is applied to visualize the results.
% Input:
%       imStack: Set of images arranged in an array as width x height x frames
%       fltParam: Struct containing the following parameters:
%          samplingFreq: Sampling rate for the frames dimension.
%          freqBands: nBands x 2 Array with the frequency bands (Hz) to be
%                     filtered. If many bands are used (nBands > 1) the output
%                     consists of nBands filtered image stacks. Do not include
%                     the Nyquist frequency.
%          filter: Cell array containing the type of window/filter and its
%                  parameters for each frequency band.
%                  'tukey': rollOff (number between 0 (rectangle) and 1 (hanning))
%          snr: Scalar. If positive represents the noise level of the spectrum,
%               i.e. the point where the snr = 1. If negative, it represents the
%               snr of the peak value of the spectrum and it is used as
%               relative-to-peak for each. HAs to be defined for each frequency band
% Output:
%       TDIm: Array containing the filtered version of the stack. There are as
%             many output arrays as frequency bands.
% hector.estrada@posteo.org

  % (0.0) Extract parameters for filtering
  nBands = size(fltParam.freqBands, 1);
  notFlattened = (length(size(imStack)) == 3);
  if notFlattened
    nFrames = size(imStack, 3);
    % (1.0) Quasi-flattening the image stack for processing
    dimensions = size(imStack);
    imStack = reshape(imStack, [prod(dimensions(1:end - 1)), dimensions(end)]);
  else
    % Already flattened stack
    nFrames = size(imStack, 2);
  end
  % (1.1) Fourier transform
  imStackFD = fftshift(ifft(imStack, [], 2), 2)*nFrames;

  % (2.0) Preparing variables
  freqFull = [-nFrames/2:nFrames/2-1]*fltParam.samplingFreq/nFrames;
  idLimits = cell(nBands, 1);
  FDWindow = cell(nBands, 1);
  for nn = 1:nBands
    [idLimits{nn}, nFreq] = find_limits(freqFull, fltParam.freqBands(nn,:));
    % (3.0) Defining the passband filter
    switch fltParam.filter{nn}{1}
      case {'tukey'}
        FDWindow{nn} = tukeywin(nFreq, fltParam.filter{nn}{2});
      otherwise
        error('This parameter is not supported at FD_filtering yet')
    end
  end
  for nn = 1:nBands
    % (4.0) Filtering the images
    % FDIm = filt_Wiener(FDWindow, imStackFD, idLimits, freqFull, fltParam.snr(nn)); % doesn't make sense in this context
    FDIm = filt_frequency_domain(FDWindow{nn}, imStackFD, idLimits{nn}, freqFull);
    % (5.0) Back to time domain
    TDIm = real(fft(ifftshift(FDIm, 2), [], 2))/length(freqFull);
    if notFlattened
      % (6.0) Final reshape to the output shape
      varargout{nn} = reshape(TDIm, dimensions);
    else
      varargout{nn} = TDIm;
    end
  end

end

function [idLimits, nFreq] = find_limits(freqFull, freqRange)
% function [idLimits, nFreq] = find_limits(freqFull, freqRange)
  % (1.0) Positive frequencies
  valid = (freqFull >= freqRange(1)) & (freqFull <= freqRange(2));
  gagnan = find(valid);
  idLimits = [min(gagnan), max(gagnan)];
  nFreq = max(gagnan) - min(gagnan) + 1;

  % (2.0) Negative frequencies
  valid = (freqFull <= -freqRange(1)) & (freqFull >= -freqRange(2));
  gagnan = find(valid);
  idLimits(3:4) = [min(gagnan), max(gagnan)];

  % (3.0) Check for overlapping regions
  if idLimits(1) <= idLimits(4)
    idLimits = [idLimits(3), idLimits(2)];% One single region
                                          % including 0 frequency
    nFreq = idLimits(2) - idLimits(1) + 1;
  end
end
