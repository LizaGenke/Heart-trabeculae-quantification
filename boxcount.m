function [ fd ] = boxcount( image )
%-------------------------------------------------------------------------
%                                                              MAIA Master
%                                                          Medical sensors  
%                                      Quantification of trabeculae inside 
%                                the heart from MRI usingf ractal analysis
%                                                 Professor: Alain LALANDE
% Authors:
% Daria Zotova
% Elizaveta Genke
% 
% This function applies box-count method to determine fractal properties
% Input: image
% Output: fractal dimension (fd)
% -------------------------------------------------------------------------
logicalImage = logical(image);
width = max(size(logicalImage));    
power = ceil(log(width)/log(2)); %highest power of 2 to determine a box-size
width = 2^power; %largest size of the box
maxBox = zeros(width, width); 
maxBox(1:size(logicalImage,1), 1:size(logicalImage,2)) = logicalImage; 
logicalImage = maxBox;
numOfElements=zeros(1,power+1);

numOfElements(power+1) = sum(logicalImage(:));
    for k=(power-1):-1:0
        size1 = 2^(power-k);
        for i=1:size1:(width-size1+1)
            for j=1:size1:(width-size1+1)
                logicalImage(i,j) = logical(sum(sum(logicalImage(i:1:(i+size1-1),j:1:(j+size1-1)))));
            end
        end
        numOfElements(k+1) = sum(sum(logicalImage(1:size1:(width-size1+1),1:size1:(width-size1+1))));
    end
    
numOfElements = numOfElements(end:-1:1);
boxSize = 2.^(0:power); % box size (1, 2, 4, 8...)
figure 
loglog(boxSize,numOfElements,'s-');
xlabel('box size'); ylabel('number of boxes');
title('Fractal Dimension');
fd = mean(-diff(log(numOfElements))./diff(log(boxSize)));
end

