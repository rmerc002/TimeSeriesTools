function  [fullMP, dataNoise] = fullMP_addNoiseAdaptive(dataOrig)
%MATRIXPROFILE_SUBLEN_GIF Summary of this function goes here
%   Detailed explanation goes here
% filename = 'fullMPInverseSum_data-testthis2_5-1-400.gif';
% startLen = 4;
% endLen = 400;
% stepLen = 4;
% range = startLen:stepLen:endLen;
% dataInv = 1./power(data,1/3);
dataOrig = normalize(dataOrig);

subLen = 50;
[matrixProfile, ~, ~, ~, mainWindow] = interactiveMatrixProfileVer2_sublenNorm(dataOrig, subLen, 0);
threshold = (max(matrixProfile)+min(matrixProfile))/4;
indices = matrixProfile < threshold;
indices = [false(ceil(subLen/2),1);indices(1:end-ceil(subLen/2))];%matrix profile indices are offset
% dataNoise(indices) = dataOrig(indices);
dataRange = max(dataOrig(indices))-min(dataOrig(indices));
noise = rand(1,length(dataOrig)).*dataRange-dataRange/2;
noise = noise*1.5;
mpNorm01 = matrixProfile-min(matrixProfile)/(max(matrixProfile)-min(matrixProfile));
mpNorm01 = [ones(1,ceil(subLen/2)), mpNorm01(1:end-ceil(subLen/2))'];
noise(1:length(mpNorm01)) = noise(1:length(mpNorm01)).*mpNorm01;
dataNoise = dataOrig+ noise;


range = [];
startLen = 10;
endLen = ceil(length(dataOrig)/20);
index = startLen;
while index < endLen
    range = [range,index];
    index = index + max(1,ceil(sqrt(index)));
end
range = [range,endLen];
fullMP = zeros(length(range),length(dataNoise));
for rangeIndex = 1:length(range)
    subLen = range(rangeIndex);
%       [matrixProfile, ~, ~, ~, mainWindow] = interactiveMatrixProfileVer2(dataNoise, subLen, 0);%int32(subLen/2));
%       [matrixProfile, ~, ~, ~, mainWindow] = interactiveMatrixProfileVer_RyanInverseSum(dataNoise, subLen, 0);%int32(subLen/2));
      [matrixProfile, ~, ~, ~, mainWindow] = interactiveMatrixProfileVer2_sublenNorm(dataNoise, subLen, 0);%int32(subLen/2));
%       [matrixProfile, ~, ~, ~, mainWindow] = interactiveMatrixProfileVer2_average(dataNoise, subLen, 0);%int32(subLen/2));

      fullMP(rangeIndex,1:length(matrixProfile)) = matrixProfile';
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

% for i=startLen+1:length(dataNoise)
%     i2 = min(endLen,length(dataNoise)-i);
% fullMP(1:i2,i) = interp1(range,fullMP(range,i),1:i2);
% end

for i=1:length(range)%endLen
    fullMP(i,end-(range(i)):end) = nan;
end
% fullMP(fullMP<min(min(fullMP(:,100:end-500)))) = nan;


ax1 = subplot(10,1,1);
plot(dataOrig);
xlim([1,length(dataOrig)]);

ax2 = subplot(10,1,2);
% x = 1:length(dataNoise);
plot(noise);
% dataNoiseTemp = nan(size(dataNoise));
% dataNoiseTemp(indices) = dataNoise(indices);
% plot(dataNoiseTemp);
% hold on;
% dataNoiseTemp = nan(size(dataNoise));
% dataNoiseTemp(~indices) = dataNoise(~indices);
% plot(dataNoiseTemp,'red');
% hold off;
xlim([1,length(dataNoise)]);

ax3 = subplot(10,1,3:10);
[X,Y] = meshgrid(1:10:length(dataNoise), range);
h = surf(X,Y,fullMP(:,1:10:end));
view(2);
set(h,'LineStyle','none');
colormap(flipud(jet));
% colormap(jet);
xlim([1,length(dataNoise)]);
ylim([1,endLen]);
% colorbar;
xlabel('Timeseries indices');
ylabel('SubLength');

linkaxes([ax1,ax2,ax3],'x');
end

