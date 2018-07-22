function I_out=regiongrowing(I,maxdist)
%-------------------------------------------------------------------------
% The function performs "region growing" from a seedpoint (x,y).
% Coordinates (x,y) are taken with getpts.
% 
% Input: image and threshold (I, maxdist).
% Output: logical output image of region (I_out).
%
% Author: D. Kroon, University of Twente.
% Code has been modified, used 8-neighbours and aoutomatic coordinates
% detection.
% ------------------------------------------------------------------------


if(exist('maxdist','var')==0)
    maxdist=0.001; 
end
[y,x]=getpts;
y=round(y(1));
x=round(x(1));

I_out = zeros(size(I)); % Prepare output logical image 
inputImageDimension = size(I); % Dimension of input image

meanOfSegmentedRegion = I(x,y); % The mean of the segmented region
pixelsInRegion = 1; % Number of pixels in region

% Free memory to store neighbours of the segmented region
neighbours = 10000; 
neighbourPosition=0;
listOfNeighbours = zeros(neighbours,3); % stores coordinates and intensity

pixdist=0; % Distance of the region newest pixel to the region mean

% 8-Neighbors location 
neighborLocation=[-1 0; -1 1; 0 1; 1 1; 1 0; 1 -1; 0 -1; -1 -1; -1 0];

% Perform regiogrowing while distance between region and posible new pixels 
% becomes higher than a treshold
while(pixdist<maxdist && pixelsInRegion<numel(I))
    % Add new neighbors pixels
    for j=1:8,
        % Get neighbour coordinates
        xn = x + neighborLocation(j,1); 
        yn = y + neighborLocation(j,2);
        
        % Check if neighbour is inside or outside the image
        inside = (xn>=1) && (yn>=1) && (xn<=inputImageDimension(1))&&(yn<=inputImageDimension(2));
        
        % Add neighbor if inside and not already part of the segmented area
        if(inside && (I_out(xn,yn)==0)) 
            neighbourPosition = neighbourPosition+1;
            listOfNeighbours(neighbourPosition,:) = [xn yn I(xn,yn)]; 
            I_out(xn,yn)=1;
        end
    end

    % Add a new block of free memory
    if(neighbourPosition+10>neighbours), neighbours=neighbours+10000; 
        listOfNeighbours((neighbourPosition+1):neighbours,:)=0;
    end
    
    % Add pixel with intensity nearest to the mean of the region to the region
    dist = abs(listOfNeighbours(1:neighbourPosition,3) - meanOfSegmentedRegion);
    [pixdist, index] = min(dist);
    I_out(x,y)=2; 
    pixelsInRegion=pixelsInRegion+1;
    
    % Calculate the new mean of the region
    meanOfSegmentedRegion = (meanOfSegmentedRegion*pixelsInRegion + listOfNeighbours(index,3))/(pixelsInRegion+1);
    
    % Save the x and y coordinates of the pixel (for the neighbour add proccess)
    x = listOfNeighbours(index,1); 
    y = listOfNeighbours(index,2);
    
    % Remove the pixel from the neighbour (check) list
    listOfNeighbours(index,:)=listOfNeighbours(neighbourPosition,:); 
    neighbourPosition=neighbourPosition-1;
end

% Return the segmented area as logical matrix
I_out=I_out>1;




