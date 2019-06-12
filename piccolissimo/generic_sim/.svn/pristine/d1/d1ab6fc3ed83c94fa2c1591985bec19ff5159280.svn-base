function [segmentedX] = curveSegmentation(X,division,low,high)
segmentation = (low:division:high)';
segmentedX = interp1q(X(:,1),X(:,2),segmentation);