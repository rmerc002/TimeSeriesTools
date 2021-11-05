function  [fullMP] = fullMP_gif(data)
%MATRIXPROFILE_SUBLEN_GIF Summary of this function goes here
%   Detailed explanation goes here
% filename = 'fullMPInverseSum_data-testthis2_5-1-400.gif';
% startLen = 4;
% endLen = 400;
% stepLen = 4;
% range = startLen:stepLen:endLen;
% dataInv = 1./power(data,1/3);
range = [];
startLen = 10;
index = startLen;
while index < 400
    range = [range,index];
    index = index + max(1,ceil(sqrt(index)));
end
endLen = range(end);
fullMP = zeros(endLen,length(data));
for subLen = range
%       [matrixProfile, ~, ~, ~, mainWindow] = interactiveMatrixProfileVer2(data, subLen, 0);%int32(subLen/2));
%       [matrixProfile, ~, ~, ~, mainWindow] = interactiveMatrixProfileVer_RyanInverseSum(data, subLen, 0);%int32(subLen/2));
      [matrixProfile, ~, ~, ~, mainWindow] = interactiveMatrixProfileVer2_sublenNorm(data, subLen, 0);%int32(subLen/2));
%       [matrixProfile, ~, ~, ~, mainWindow] = interactiveMatrixProfileVer2_average(data, subLen, 0);%int32(subLen/2));

      fullMP(subLen,1:length(matrixProfile)) = matrixProfile';
      drawnow 
      % Capture the plot as an image 
      frame = getframe(mainWindow.fig); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,16); 
      % Write to the GIF File 
      
%           if subLen == startLen 
%               imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',1/30); 
%           else 
%               if mod(subLen,5) == 0
%                     imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',1/30); 
%               end 
%           end
      close all;
end

for i=startLen+1:length(data)
    i2 = min(endLen,length(data)-i);
fullMP(1:i2,i) = interp1(range,fullMP(range,i),1:i2);
end

for i=1:endLen
    fullMP(i,end-(i):end) = nan;
end
fullMP(fullMP<min(min(fullMP(:,100:end-500)))) = nan;

% mask = fullMP>0.51;
% fullMP2 = fullMP - 0.5;
% fullMP2(~mask) = nan;

ax1 = subplot(10,1,1);
plot(data);
xlim([1,length(data)]);

ax2 = subplot(10,1,2:10);
h = surf(fullMP);
view(2);
set(h,'LineStyle','none');
colormap(flipud(jet));
% colormap(jet);
xlim([1,length(data)]);
ylim([1,endLen]);
% colorbar;
xlabel('Timeseries indices');
ylabel('SubLength');

linkaxes([ax1,ax2],'x');
end

