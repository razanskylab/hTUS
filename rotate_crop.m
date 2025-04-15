function [outMatrix] = rotate_crop(inMatrix, rotAngle, interp, cropVector)
% function [outMatrix] = rotate_crop(inMatrix, rotAngle, cropVector)
% Rotates by an angle rotAngle (deg) using imrotate with loose croping, interp
% interplation method (string) and uses cropVector: [Idx1, Idx2, Idy1, Idy2]
% as custom cropping range.
% hector.estrada@posteo.org

rotAux = imrotate(inMatrix, rotAngle, interp, 'loose');
outMatrix = rotAux(cropVector(1):cropVector(2), cropVector(3):cropVector(4));
