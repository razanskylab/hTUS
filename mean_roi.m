function [vectorMean] = mean_roi(inputArray);
% function [vectorMean] = mean_roi(volume);
% Function to calculate the mean along the first array dimensions of the inputArray
% Imput:
%       inputArray: Set of images arranged in an array as 
% Output: 
%       vectorMean: vector containing the mean 
% hector.estrada@posteo.org

% (1.0) Arranging the array dimensions
dimensions = size(inputArray);

% (2.0) Reshaping and allocating 
Y = reshape(inputArray, [prod(dimensions(1:end - 1)), dimensions(end)]);% Quasi flattened 
vectorMean = zeros(1, dimensions(end));

% (3.0) Calculation
vectorMean = mean(Y, 1);
