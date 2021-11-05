function dataOut = slidingNorm(data,windowSize)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
dataOut = zeros(size(data));
for i=1:length(data)-windowSize+1
    dataOut(i:i+windowSize-1) = normalize(data(i:i+windowSize-1),'range');
end
dataOut = dataOut/windowSize;
end

