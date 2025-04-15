function preview_image_stack(param)
% function preview_image_stack(param)
% Function to determine optimal parameters for the real
%
% param.Crop.BaseLineRange Use imageJ
% param.Crop.LastUsefulFrame Use imageJ

% param.Crop.MaskThreshold
% param.Crop.xy
% param.Rotate.MidLine
% param.Rotate.Cropxy
% param.Visualize.SelectedPoint
% TODO use the function [xi,yi] = getpts
% [xi,yi] = getline
% rect = getrect


close all
disp('This is the image stack preview. At the end of the procedure, you can');
disp('copy/paste the output of the command line into your script.');

% (1.0) Load tiff file
[folderName, fileName, ext] = fileparts(param.File.FullName);
try
  infoImage = imfinfo(param.File.FullName);
  mImage = infoImage(1).Width;
  nImage = infoImage(1).Height;
  numberImages = length(infoImage);% before the baseline are removed
  im = zeros(nImage, mImage, numberImages, 'uint16');
  for i = 1:numberImages
      im(:,:,i) = imread(param.File.FullName, 'Index', i, 'Info', infoImage);
  end
  im = double(im);
catch error
  warning("Error during tif file read. " + error.message);
end

% (2.0) Cropping the image
figure(1);
aux =  squeeze(mean(im, 3));
imagesc(aux);
axis image
colormap(gray(129));
title('Mean from original')
h = images.roi.Rectangle();
disp('Please draw a rectangle to crop the image')
h.draw();
cropFlag = input('Press enter once you are ready','s');
cropWindow = round(h.Position);
cropWindowId = [cropWindow(2), cropWindow(2) + cropWindow(4), cropWindow(1), ...
                cropWindow(1) + cropWindow(3)];

im = im(cropWindowId(1):cropWindowId(2), cropWindowId(3):cropWindowId(4), :);
aux = squeeze(mean(im, 3));
close;
disp('Image cropping ready!')
disp(' ')

% (3.0) Base-line and last useful frame
fig = figure();
ax = subplot(1,2,1);
imagesc(aux)
axis image
colormap(gray)

disp('Please select a point in the image. When sure of last point press Last, and then choose.')
lastButton = uicontrol(fig,'Style','togglebutton','String','Last');
lastButton.Value = false;

h = images.roi.Point(ax);
h.draw();
selectedPoint = fliplr(round(h.Position));
% (3.1) Calculation on selected points
dx = (param.Visualize.SpotSize(1) - 1)/2;
dy = (param.Visualize.SpotSize(2) - 1)/2;
ROIx = selectedPoint(1) - dx:selectedPoint(1) + dx;
ROIy = selectedPoint(2) - dy:selectedPoint(2) + dy;
timeTrace0 = mean_roi(im(ROIx, ROIy, :));
ax2 = subplot(1,2,2);
plot(timeTrace0)
hold on

  while ~lastButton.Value
  h = images.roi.Point(ax);
  h.draw();
  selectedPoint = fliplr(round(h.Position));
  % (3.1) Calculation on selected points
  dx = (param.Visualize.SpotSize(1) - 1)/2;
  dy = (param.Visualize.SpotSize(2) - 1)/2;
  ROIx = selectedPoint(1) - dx:selectedPoint(1) + dx;
  ROIy = selectedPoint(2) - dy:selectedPoint(2) + dy;
  timeTrace0 = mean_roi(im(ROIx, ROIy, :));
  plot(ax2, timeTrace0)
  end

% (3.1) Base-line range
disp('Now select the region for the base-line')
h2 = images.roi.Rectangle(ax2);
h2.draw();
baseLineFlag = input('Press enter once you are ready ','s');
baseLineRange = round([h2.Position(1), h2.Position(1) + h2.Position(3)]);
disp('Base-line range ready!')
disp(' ')

% (3.2) Last useful frame
disp('Now select the last useful frame')
h3 = images.roi.Point(ax2);
h3.draw();
LastFrameFlag = input('Press enter once you are ready ','s');
lastUsefulFrame = round(h3.Position(1));
disp('Last useful frame ready!')
disp(' ')
close

% (4.0) Mask threshold
baseLine = mean(im(:, :, baseLineRange(1):baseLineRange(2)), 3);
minAmp = min(im(:, :, baseLineRange(1):end),[], 3);
testMaskThreshold = input('Enter the relative mask threshold: (default 0.3) ');
if isempty(testMaskThreshold)
  testMaskThreshold = 0.3;
end
mask = baseLine > testMaskThreshold*max(baseLine(:));
while ~isempty(testMaskThreshold)
  maskThreshold = testMaskThreshold; % Saving the last one
  mask = baseLine > maskThreshold*max(baseLine(:));
  figure()
  imagesc(baseLine.*mask)
  axis image
  colormap(gray)
  testMaskThreshold = input('Does it look good? the last will be saved: ');
end
disp('Mask threshold ready!')
disp(' ')
close all
% % (4.1) Don't diverge option
% if maskThreshold == 0
%   figure()
%   subplot(2,1,1)
%   plot(im(:), 1./im(:), '.')
%   title('1/x')
%   subplot(2,1,2)
%   hist(im(:), 500)
%   title('pixel value distribution')
% end


% (5.0) Brain's midline and rotation
figure()
subplot(1, 2 , 1)
imagesc(baseLine.*mask)
axis image
colormap(gray)
h = images.roi.Line();

disp('Please draw a line along the brain''s midline to rotate the image')
disp('First point at the rear of the brain')
h.draw();
lineFlag = input('Press enter once you are ready','s');
midLine = h.Position;
% (5.1) Calculating the angle and preparing the rotation
rotAngle = atan2((midLine(2,2) - midLine(1,2)), ...
                 (midLine(2,1) - midLine(1,1)))*180/pi + 90
rotatedIm = imrotate(minAmp.*mask, rotAngle, 'bilinear', 'loose');
disp('Midline ready!')
disp(' ')

ax = subplot(1, 2, 2);
imagesc(rotatedIm)
axis image
colormap(gray)
h = images.roi.Rectangle(ax);
disp('Please draw a rectangle to crop the image')
h.draw();
cropFlag = input('Press enter once you are ready','s');
cropWindow = round(h.Position)
cropXY = [cropWindow(2), cropWindow(2) + cropWindow(4), cropWindow(1),...
                cropWindow(1) + cropWindow(3)];
cropRotated = rotatedIm(cropXY(1):cropXY(2), cropXY(3):cropXY(4));
close;

% (5.2) Rotating and croping the whole dataset to properly find the selected point
auxIm = im(:,:, baseLineRange(2):lastUsefulFrame);
newSize = [cropXY(2) - cropXY(1) + 1,...
           cropXY(4) - cropXY(3) + 1,...
           size(auxIm, 3)];

im = zeros(newSize);
for nn = 1:size(auxIm, 3)
  % (4.3) Rotate and crop unfiltered fluorescence changes
  im(:, :, nn) = rotate_crop(squeeze(auxIm(:,:,nn)), rotAngle, ...
                               'bilinear', cropXY);
end

fig2 = figure()
ax = subplot(1,2,1);
imagesc(cropRotated)
axis image
colormap(gray)
disp('Image cropping after rotation ready!')
disp(' ')

% (6.0) Point of interest

disp('Please select a point in the image. When sure of last point press Last, and then choose.')
lastButton2 = uicontrol(fig2,'Style','togglebutton','String','Last');
lastButton2.Value = false;

h = images.roi.Point(ax);
h.draw();
selectedPoint = fliplr(round(h.Position));
% (6.1) Calculation on selected points
dx = (param.Visualize.SpotSize(1) - 1)/2;

dy = (param.Visualize.SpotSize(2) - 1)/2;
ROIx = selectedPoint(1) - dx:selectedPoint(1) + dx;
ROIy = selectedPoint(2) - dy:selectedPoint(2) + dy;
timeTrace0 = mean_roi(im(ROIx, ROIy, :));
ax2 = subplot(1,2,2);
plot(timeTrace0)
hold on
while ~lastButton2.Value
  h = images.roi.Point(ax);
  h.draw();
  selectedPoint = fliplr(round(h.Position));
  % (6.1) Calculation on selected points
  dx = (param.Visualize.SpotSize(1) - 1)/2;
  dy = (param.Visualize.SpotSize(2) - 1)/2;
  ROIx = selectedPoint(1) - dx:selectedPoint(1) + dx;
  ROIy = selectedPoint(2) - dy:selectedPoint(2) + dy;
  timeTrace0 = mean_roi(im(ROIx, ROIy, :));
  plot(ax2, timeTrace0)
end
close
disp('Point of interest ready')

% (7.0) Print the results to be copied in the make_fluo_data
disp('')
disp('You can now copy/paste the following lines into your script')
disp('')

lines{1} = ['param.Crop.BaseLineRange = [',num2str(baseLineRange),'];\n'];
lines{2} = ['param.Crop.LastUsefulFrame = [',num2str(lastUsefulFrame),'];\n'];
lines{3} = ['param.Crop.MaskThreshold = [',num2str(maskThreshold),'];\n'];
lines{4} = ['param.Crop.xy = [',num2str(cropWindowId),'];\n'];
lines{5} = ['param.Rotate.MidLine = [',num2str(midLine(1, :)),';\n'];
lines{6} = ['                        ',num2str(midLine(2, :)),'];\n'];
lines{7} = ['param.Rotate.Cropxy = [',num2str(cropXY),'];\n'];
lines{8} = ['param.Visualize.SelectedPoint = [',num2str(selectedPoint),'];\n'];

for nn = 1:8
  fprintf(1, lines{nn});
%}
end
end
