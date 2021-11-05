function  [magicMP,profileIndices] = magicMatrixProfileVisualizer_Anytime(data)
% close all;

data = data/max(max(data));
dataLen = length(data);


%add noise in order to increase contrast where the signal is low.
dataNoise = amplifyDissimilarities(data);
plot(dataNoise);
pause;
close all;
f1 = figure;

%%% Determine subLengths to iterate over
startLen = 10;
endLen = ceil(length(data)/20);
if length(data) < 8000
    endLen = 600;
end
subLenSeries = getSubLenSeries(startLen, endLen);


magicMP = ones(length(subLenSeries),length(data))*0.3;
% magicMP(1,1) = 0;
profileIndices = ones(length(subLenSeries),length(data));

screenXRes = 2000;
step = floor(length(data)/screenXRes);
% magicMPDisplay = nan(length(subLenSeries),floor(length(data)/step));
% magicMPDisplay(1,1) = 0;
% magicMPDisplay(11:20,101:200) = ones(10,100);
% magicMPDisplay(11:20,201:300) = zeros(10,100);


ax1 = subplot(10,1,1);

ax2 = subplot(10,1,2);
plot(data,'Color',[.5,.5,.5]);
xlim([1,screenXRes]);

%give an initial guess 
ax3 = subplot(10,1,3:10);
plotSurface(magicMP,subLenSeries);

linkaxes([ax2,ax3],'x');
adjustContrast = true;

tic
for i = 1:200
    tic
    subLenIndex = randi([1,length(subLenSeries)]);
    subLen = subLenSeries(subLenIndex);
    idx = randi([1,dataLen-subLen]);
    query = data(idx:idx+subLen-1);
    proLen = dataLen - subLen + 1;
    [dataFreq, dataMu, dataSig] = massPre(dataNoise, dataLen, subLen);
    
    
    distProfile = MASS_V2(dataNoise, query);
    max(distProfile)
    min(distProfile)
%     distProfile = abs(distProfile);
%     distProfile = sqrt(distProfile);

%     distProfile(isSkip) = inf;
    excZoneLen = round(subLen * 0.5);
    excZoneStart = max(1, idx - excZoneLen);
    excZoneEnd = min(proLen, idx + excZoneLen);
    distProfile(excZoneStart:excZoneEnd) = inf;

%     disp('distProfile');
%     size(distProfile)
    updatePos = distProfile < magicMP(subLenIndex,1:size(distProfile,2));
    profileIndices(subLenIndex,updatePos) = idx;
    magicMP(subLenIndex,updatePos) = distProfile(updatePos);
    [magicMP(subLenIndex,idx), profileIndices(subLenIndex,idx)] = min(distProfile);

%     distProfile = diagonalDist(...
%             data, idx, dataLen, subLen, proLen, dataMu, dataSig);
%     distProfile = abs(distProfile);
%     distProfile = sqrt(distProfile);
% 
%     pos1 = idx:proLen;
%     pos2 = 1:proLen - idx + 1;
% 
%     updatePos = magicMP(subLenIndex,pos1) > distProfile;
%     profileIndices(subLenIndex,pos1(updatePos)) = pos2(updatePos);
%     magicMP(subLenIndex,pos1(updatePos)) = distProfile(updatePos);
%     updatePos = magicMP(subLenIndex,pos2) > distProfile;
%     profileIndices(subLenIndex,pos2(updatePos)) = pos1(updatePos);
%     magicMP(subLenIndex,pos2(updatePos)) = distProfile(updatePos);

%     matrixProfile(isSkip) = inf;
%     profileIndex(isSkip) = 0;
    
    toc
    
%     [timeseriesIndex, subLength] = ginput(1);
%     
%     if checkContrastButtonClicked(timeseriesIndex, subLength, magicMP, subLenSeries)
%         adjustContrast = ~adjustContrast;
%         continue
%     end
%     
%     timeseriesIndex = uint32(timeseriesIndex);
%     [~,subLengthIndex] = max(subLenSeries((subLenSeries-subLength)<= 0));
%     subLength = subLenSeries(subLengthIndex);
%     
%     subplot(10,1,3:10);
%     if adjustContrast == true
%         threshold = magicMP(subLengthIndex, timeseriesIndex);
%         magicMP_contrast = centerOnThreshold(magicMP, threshold);
%         plotSurface(magicMP_contrast,subLenSeries);
%     end
    
    
    plotSurface(magicMP, subLenSeries);
    
%     %%% Plot nearest neighbor subsequence
%     startIndex1 = timeseriesIndex;
%     startIndex2 = profileIndices(subLengthIndex,startIndex1);
%     plotNearestNeighborSubsequences(data, profileIndices, startIndex1, startIndex2, subLength);
%     
%     %%% Plot original data with highlighted subsequences
%     subplot(10,1,2);
%     plot(data,'Color',[.5,.5,.5]);
%     xlim([1,length(data)]);
%     hold on;
%     plot(startIndex1:startIndex1 + subLength,data(startIndex1:startIndex1 + subLength),'g');
%     plot(startIndex2:startIndex2 + subLength,data(startIndex2:startIndex2 + subLength),'r');
%     hold off;
%     
%     %give some descriptive data
%     subplot(10,1,3:10);
%     xlabel(sprintf('index1 = %d, index2 = %d, subLength = %d, dist = %f',startIndex1, startIndex2, subLenSeries(subLengthIndex), magicMP(subLengthIndex,timeseriesIndex)));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   HELPER FUNCTIONS   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotSurface(magicMP, subLenSeries)
    magicMP = drawContrastButton(magicMP, subLenSeries);
    
    
    startIndex = 1;%2000;
    endIndex = size(magicMP,2);%4000;
%     endIndex = size(magicMP,2)
    step = floor((endIndex-startIndex)/1024);
    plotIndices = startIndex:step:endIndex;
    
    magicMPDisplay = zeros(length(subLenSeries), length(plotIndices));
    floor(size(magicMP,2)/10)
    for index=1:floor(size(magicMP,2)/10)
        sIndex = 1+(index-1)*10;
        eIndex = sIndex + 10 - 1;
        magicMPDisplay(:,index) = min(magicMP(:,sIndex:eIndex),[],2);
    end
    
    [X,Y] = meshgrid(plotIndices, subLenSeries);
%     h = surf(X,Y,magicMPDisplay);
    h = surf(X,Y,magicMP(:,1:step:end));
    view(2);
    set(h,'LineStyle','none');
    c1 = hot(150);
    c1 = c1(1:100,:);
    c2 = flipud(winter(100));
    customColor = [c1;ones(5,3);c2];
    colormap(customColor);
    
%     colormap(flipud(jet));
    % colormap(jet);
    xlim([startIndex,endIndex]);
    ylim([1,subLenSeries(end)]);
%     colorbar;
    xlabel('Timeseries indices');
    ylabel('SubLength');
    drawnow;
    
end

function magicMP = drawContrastButton(magicMP, subLenSeries)
    [x,y] = getButtonBounds(size(magicMP,2), subLenSeries(end));
    series = 1:length(subLenSeries);
    index = min(series(subLenSeries >= y));
    magicMP(index:end, x:end) = 0;
end

function buttonClicked = checkContrastButtonClicked(xInput, yInput, magicMP, subLenSeries)
    buttonClicked = false;
    [x,y] = getButtonBounds(size(magicMP,2), subLenSeries(end));
    if xInput >= x && yInput >= y
        buttonClicked = true;
    end
end

function [x, y] = getButtonBounds(xMax,yMax)
    x = ceil(xMax - xMax/25);
    y = ceil(yMax - yMax/10);
    
end

function [magicMP_contrast] = centerOnThreshold(magicMP, threshold)
    magicMP_contrast = magicMP - threshold;
    %matching indices
    indices = magicMP <= threshold;
    magicMP_contrast(indices) = magicMP_contrast(indices) - min(magicMP_contrast(indices));
    magicMP_contrast(indices) = magicMP_contrast(indices)./max(max(magicMP_contrast(indices)));
    magicMP_contrast(indices) = magicMP_contrast(indices)./2.5;

    %differing indices
    indices = ~indices;
    magicMP_contrast(indices) = magicMP_contrast(indices) - min(magicMP_contrast(indices));
    magicMP_contrast(indices) = magicMP_contrast(indices)./max(max(magicMP_contrast(indices)));
    magicMP_contrast(indices) = magicMP_contrast(indices)./2.5 + 0.6;

end

function plotNearestNeighborSubsequences(data, profileIndices, index1, index2, subLength)
    subplot(10,1,1);
    
    startIndex = index1;
    endIndex = startIndex+subLength;

    startIndex2 = index2;
    endIndex2 = startIndex2 + subLength;
    
    plot(zscore(data(startIndex:endIndex)),'g');
    hold on;
    plot(zscore(data(startIndex2:endIndex2)),'r');
    hold off;
    xlim([1,subLength]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   MagicMP HELPER FUNCTIONS   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dataNoise = amplifyDissimilarities(data)
%%% Purpose: Reduce similarity due to low amplitude noise.
%%% If there are two subsequences, A & B, where A has low distance and B
%%% has high distance, when they are concatenated, the resulting
%%% subsequence cannot have less distance than A. B can only make A worse.
%%% I want to penalize small subsequences with low similarity
%%%
%%% Applicaiton: This encourages matching of non-noise motifs. Visually,
%%% there is a convergence in the form of a peak in the 
%%% Magic Matrix Profile plot.
    
    %Use matrix profile to get subsequence similarity
    subLen = 10;
%     [matrixProfile, ~, ~, ~, ~] = interactiveMatrixProfileVer2_xlogx_sublenNorm(data, subLen, 0);
    [matrixProfile,~] = mpx_correlation(data',ceil(subLen/2),subLen);
%     matrixProfile = 1-matrixProfile;
    

    
    %Determine amplitude of noise to add
    threshold = (max(matrixProfile)-min(matrixProfile))/2 + min(matrixProfile);
    indices = matrixProfile < threshold;
    indices = [false(ceil(subLen/2),1);indices(1:end-ceil(subLen/2))];%matrix profile indices are offset
    dataRange = max(data(indices))-min(data(indices));
%     dataRange = max(data)-min(data);
    
    

    noise = rand(1,length(data));%randFunc(length(data));

    noise = noise.*dataRange-dataRange/2;
    noise = noise*1.5;
    mpNorm01 = matrixProfile-min(matrixProfile);
    mpNorm01 = mpNorm01./max(mpNorm01);

    % mpNorm01(mpNorm01>0.5) = 1.5-mpNorm01(mpNorm01>0.5);
    mpNorm01 = [ones(1,ceil(subLen/2)), mpNorm01(1:end-ceil(subLen/2))'];
    noise(1:length(mpNorm01)) = noise(1:length(mpNorm01)).*mpNorm01;
    dataNoise = data+ noise;
    
%     subplot(3,1,1);
%     plot(data);
%     subplot(3,1,2);
%     plot(matrixProfile);
%     subplot(3,1,3);
%     plot(dataNoise);
%     pause;
end

function subLenSeries = getSubLenSeries(startLen, endLen)
%%% Purpose: Reduce space. matrix profile distances are pretty similar from
%%% from one subLen to the next. By trial and error, adding the square root
%%% of the current subLen seems to produce a good distribution.
    subLenSeries = [];
    index = startLen;
    
    while index < endLen
        subLenSeries = [subLenSeries,index];
        index = index + max(1,ceil(sqrt(index)));
    end
    
    subLenSeries = [subLenSeries,endLen];
end

% This code is created by Abdullah Mueen, Yan Zhu, Michael Yeh, Kaveh Kamgar, Krishnamurthy Viswanathan, Chetan Kumar Gupta and Eamonn Keogh.
% The overall time complexity of the code is O(n log n). The code is free to use for research purposes.
% The code may produce imaginary numbers due to numerical errors for long time series where batch processing on short segments can solve the problem.

function [dist] = MASS_V2(x, y)
%x is the data, y is the query
m = length(y);
n = length(x);

%compute y stats -- O(n)
meany = mean(y);
sigmay = std(y,1);

%compute x stats -- O(n)
meanx = movmean(x,[m-1 0]);
sigmax = movstd(x,[m-1 0],1);

y = y(end:-1:1);%Reverse the query
y(m+1:n) = 0; %aappend zeros

%The main trick of getting dot products in O(n log n) time
X = fft(x);
Y = fft(y);
Z = X.*Y;
z = ifft(Z);

dist = 2*(m-(z(m:n)-m*meanx(m:n)*meany)./(sigmax(m:n)*sigmay));
dist = sqrt(dist);
dist = dist./sqrt(length(y));
end

function distProfile = mass(dataFreq, query, ...
    dataLen, subLen, dataMu, dataSig, queryMu, querySig)
query = query(end:-1:1);
query(subLen+1:(subLen+dataLen)) = 0;
queryFreq = fft(query);
productFreq = dataFreq .* queryFreq;
product = ifft(productFreq);
distProfile = 2 * (subLen - ...
    (product(subLen:dataLen) - subLen * dataMu * queryMu) ./ ...
    (dataSig * querySig));
end

function distProfile = diagonalDist(...
    data, idx, dataLen, subLen, proLen, dataMu, dataSig)
xTerm = ones(proLen - idx + 1, 1) * ...
    (data(idx:idx + subLen - 1)' * data(1:subLen));
mTerm = data(idx:proLen - 1) .* ...
    data(1:proLen - idx);
aTerm = data(idx + subLen:end) .* ...
    data(subLen + 1:dataLen - idx + 1);
if proLen ~= idx
    xTerm(2:end) = xTerm(2:end) - cumsum(mTerm) + cumsum(aTerm);
end
end


function [dataFreq, dataMu, dataSig] = massPre(data, dataLen, subLen)
data(dataLen + 1:(subLen + dataLen)) = 0;
dataFreq = fft(data);
dataCumsum = cumsum(data);
data2Cumsum =  cumsum(data .^ 2);
data2Sum = data2Cumsum(subLen:dataLen) - ...
    [0; data2Cumsum(1:dataLen - subLen)'];
dataSum = dataCumsum(subLen:dataLen) - ...
    [0; dataCumsum(1:dataLen - subLen)'];
dataMu = dataSum ./ subLen;
data2Sig = (data2Sum ./ subLen) - (dataMu .^ 2);
dataSig = sqrt(data2Sig);
end