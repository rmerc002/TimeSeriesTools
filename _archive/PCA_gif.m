function PCA_gif(data,range)
%MATRIXPROFILE_SUBLEN_GIF Summary of this function goes here
%   Detailed explanation goes here
filename = '/Users/cdslug/Dropbox/Education/UCR/Research/0_Projects/Rastamat/test_gifs/test.gif';
% filename = '/Users/cdslug/Dropbox/Education/UCR/Research/0_Projects/Rastamat/test_gifs';

lastX = zeros(size(data));
lastY = zeros(size(data));
for subLen = range
    fprintf('Running sublen=%d\n',subLen);
    dataSegments = zeros(length(data)-subLen+1,subLen);
    for i=1:length(data)-subLen+1
        dataSegments(i,:) = data(i:i+subLen-1);
    end
    dataSegments = normalize(dataSegments,2);
%     indexColors = 1:size(dataSegments,1);
    indexColors = randperm(size(dataSegments,1));
    dataSegments = dataSegments(indexColors,:);
    % UPDATED SECTION
    [coeff,score,latent,tsquared,explained] = pca(dataSegments);
%     close all;
    figure;
%     thisX = score(1:3:end,1);
%     thisY = score(1:3:end,2);
    thisX = score(:,1);
    thisY = score(:,2);
    % END UPDATED SECTION
    
    minLen = min(length(lastX),length(thisX));
    diffSame = sqrt((lastX(1:100:minLen)-thisX(1:100:minLen)).^2);
    diffReversed = sqrt((lastX(minLen:-100:1)-thisX(1:100:minLen)).^2);
    if diffReversed < diffSame
        thisX = thisX(end:-1:1);
    end
    
    minLen = min(length(lastX),length(thisX));
    diffSame = sqrt((lastY(1:100:minLen)-thisY(1:100:minLen)).^2);
    diffReversed = sqrt((lastY(minLen:-100:1)-thisY(1:100:minLen)).^2);
    if diffReversed < diffSame
        thisY = thisY(end:-1:1);
    end
    
    lastX = thisX;
    lastY = thisY;
    scatter(thisX, thisY, 1,indexColors,'.');
    colormap(jet);
%     xlim([-0.02,0.02]);
%     ylim([-.02,0.02]);
    
    
    
    drawnow 
    % Capture the plot as an image 
    frame = getframe();%mainWindow.fig); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,16); 
    % Write to the GIF File 

      if subLen == range(1)
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',1/15); 
      else 
%           if mod(subLen,5) == 0
          imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',1/15); 
%           end 
      end
    
end
% close all;
end

