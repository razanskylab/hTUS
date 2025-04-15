function [varargout] = relative_change(varargin);
% function [dF_F, mask] = relative_change(imStack, baseLineRange, maskThreshold, imReference);
% Function to calculate relative changes on an image stack based on a predefined
% range. It also returns a binary mask that is applied to visualize the results.
%
% Input:
%       imStack: Set of images arranged in an array as width x height x frames
%       baseLineRange: Range index where no excitation happens. If it is negative, 
%                      the baseline is not removed from the output.
%       maskThreshold: Relative to total stack maximum defining the mask. If no
%                      masking is required, use maskThreshold = 0. If
%                      maskThreshold is negative, the difference is weighted
%                      to counteract divergent values from denominators close to
%                      zero.
%       imReference (optional): Image stack containing a reference to apply moving
%                      reference. Should have the same dimensions as imStack.
%       nStandard (optional): Factor accounting for using the standard deviation
%                             to weight the mean and account for oscillations
%                             present in the baseLine.
%                         reference = mean(baselLine) + nStandard*std(baseline);
% Output:
%       dF_F: Array containing the relative changes. It does include the
%             baseLineRange in the output depending on the sign of baseLineRange.
%       mask: width x height binary array masking dF_F. Retrieved only if
%             maskThreshold is not zero.
% hector.estrada@posteo.org

% (0.0) Getting input parameters
imStack = varargin{1};
baseLineRange = varargin{2};
maskThreshold = varargin{3};
if baseLineRange(1) < 0
  keepBaseLineFlag = 1;
  baseLineRange = abs(baseLineRange);
  noBaselineFrames = 0;
else
  keepBaseLineFlag = 0;
  noBaselineFrames = diff(baseLineRange);
end
movingReferenceFlag = 0;
nStandard = 0;
if length(varargin) == 4
  if prod(size(varargin{4})) > 1
    movingReferenceFlag = 1;
    imReference = varargin{4};
  else
    nStandard = varargin{4};
  end
end

% (1.0) Calculating the global_mean baseline

notFlattened = (length(size(imStack)) == 3);
if notFlattened
  nFrames = size(imStack, 3);
  baseLineMean = mean(imStack(:, :, 1:diff(baseLineRange) + 1), 3);
  baseLineSTD = std(imStack(:, :, 1:diff(baseLineRange) + 1), 1, 3);
  % (2.0) Raw calculation
  if ~movingReferenceFlag
    referenceOffset = repmat(baseLineMean, ...
                                         [1, 1, nFrames - noBaselineFrames]);
    referenceScale = repmat(baseLineMean + nStandard*baseLineSTD , ...
                                         [1, 1, nFrames - noBaselineFrames]);
  else
    referenceOffset = imReference(:,:, noBaselineFrames + 1:nFrames);
    referenceScale = referenceOffset;
  end
  % (3,0) Mask
  if (maskThreshold > 0)
    dF_F = (imStack(:,:, noBaselineFrames + 1:nFrames) - referenceOffset)./ ...
           referenceScale;
    mask = baseLineMean > maskThreshold*max(baseLineMean(:));
    % Commented to avoid redundant masking.
    %  newFrames = size(dF_F, 3);
    %  megaMask = repmat(mask, [1, 1 , newFrames]);
    %  dF_F = dF_F.*megaMask;
    varargout{2} = mask;
  elseif (maskThreshold < 0)
    % % Debugging
    % disp('Hey, you better don''t diverge')
    % x = linspace(min(abs(referenceScale(:))), max(referenceScale(:)), 100);
    % ddv = x./(x.^2 - 1/maskThreshold);
    % figure()
    % plot(x, 1./x,'k.')
    % hold on
    % plot(x, ddv, 'r')
    % axis tight
    % ylim([0 max(ddv)]*5)
    % drawnow
    % pause()
    % close()
    dontDiverge = referenceScale./(referenceScale.^2 - 1/maskThreshold);
    dF_F = (imStack(:,:, noBaselineFrames + 1:nFrames) - referenceOffset).* ...
           dontDiverge;
  else
    dF_F = (imStack(:,:, noBaselineFrames + 1:nFrames) - referenceOffset)./ ...
           referenceScale;
  end
  varargout{1} = dF_F;
else
  %  Input imStack flattened
  nFrames = size(imStack, 2);
  baseLineMean = mean(imStack(:, 1:diff(baseLineRange) + 1), 2);
  baseLineSTD = std(imStack(:, 1:diff(baseLineRange) + 1), 1, 2);
  % (2.0) Raw calculation
  if ~movingReferenceFlag
    referenceOffset = repmat(baseLineMean, [1, nFrames - noBaselineFrames]);
    referenceScale = repmat(baseLineMean + nStandard*baseLineSTD, ...
                                           [1, nFrames - noBaselineFrames]);
  else
    referenceOffset = imReference(:, noBaselineFrames + 1:nFrames);
    referenceScale = referenceOffset;
  end
  % (3,0) Mask
  if (maskThreshold > 0)
    dF_F = (imStack(:, noBaselineFrames + 1:nFrames) - referenceOffset)./ ...
           referenceScale;
    mask = baseLineMean > maskThreshold*max(baseLineMean(:));
%      newFrames = size(dF_F, 3);
%      megaMask = repmat(mask, [1, 1 , newFrames]);
%      dF_F = dF_F.*megaMask;
    varargout{2} = mask;
  elseif (maskThreshold < 0)
    % % Debugging
    % disp('Hey, you better don''t diverge')
    % x = linspace(min(abs(referenceScale(:))), max(referenceScale(:)), 100);
    % ddv = x./(x.^2 - 1/maskThreshold);
    % figure()
    % plot(x, 1./x,'k.')
    % hold on
    % plot(x, ddv, 'r')
    % axis tight
    % ylim([0 max(ddv)]*5)
    % drawnow
    % pause()
    % close()
    dontDiverge = referenceScale./(referenceScale.^2 - 1/maskThreshold);
    dF_F = (imStack(:, noBaselineFrames + 1:nFrames) - referenceOffset).* ...
           dontDiverge;
  else
    dF_F = (imStack(:, noBaselineFrames + 1:nFrames) - referenceOffset)./ ...
            referenceScale;
  end
  varargout{1} = dF_F;

  % (4.0) Visualization for validation and debugging
  % (4.1) Figuring our time
  % timeFull = 0:nFrames - 1;
  % timeShort = diff(baseLineRange) + 1:nFrames;
  % % (4.2) Plotting
  % figure()
  % if ~movingReferenceFlag
  %   subplot(1,2,2)
  %   yyaxis('right')
  %   plot(timeShort, dF_F, 'r', 'LineWidth', 1.5);
  %   yyaxis('left')
  %   plot(timeFull, imStack, 'k-', 'LineWidth', 1.5);
  %   hold on
  %   plot([1 diff(baseLineRange) + 1], referenceOffset(1:2), 'b', 'LineWidth', 1.5);
  %   if (maskThreshold >= 0)% Don't diverge
  %     plot([1 diff(baseLineRange) + 1], referenceScale(1:2), 'g', 'LineWidth', 1.5);
  %   else
  %     plot([1 diff(baseLineRange) + 1], 1./dontDiverge(1:2), 'g', 'LineWidth', 1.5);
  %   end
  %   legend('Input', 'Reference offset', 'Reference scale', '\Delta F/F (%)')
  % else
  %   yyaxis('right')
  %   plot(timeShort, dF_F, 'r', 'LineWidth', 1.5);
  %   yyaxis('left')
  %   plot(timeFull, imStack, 'k-', 'LineWidth', 1.5);
  %   hold on
  %   plot(timeShort, referenceOffset, 'b', 'LineWidth', 1.5);
  %   legend('Input', 'Reference', '\Delta F/F')
  % end
  %
  % if ~movingReferenceFlag
  %   subplot(1,2,1)
  %   hist(imStack(:, 1:diff(baseLineRange) + 1), 200);
  %   hold on
  %   plot([baseLineMean, baseLineMean], [0 5])
  %   plot([baseLineMean + nStandard*baseLineSTD, ...
  %         baseLineMean + nStandard*baseLineSTD], [0 5])
  %   title('Base line distribution')
  % end

end
