function  [magicMP, profileIndices,range] = magicMatrixProfile(dataOrig)
dataOrig = normalize(dataOrig);

%add noise in order to increase contrast where the signal is low.
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
magicMP = zeros(length(range),length(dataNoise));
profileIndices = zeros(length(range),length(dataNoise));
for rangeIndex = 1:length(range)
    subLen = range(rangeIndex);

      [matrixProfile, profileIndex, ~, ~, mainWindow] = interactiveMatrixProfileVer2_sublenNorm(dataNoise, subLen, 0);

      magicMP(rangeIndex,1:length(matrixProfile)) = matrixProfile';
      profileIndices(rangeIndex,1:length(profileIndex)) = profileIndex';
      
      close all;
end

% for i=startLen+1:length(dataNoise)
%     i2 = min(endLen,length(dataNoise)-i);
% magicMP(1:i2,i) = interp1(range,magicMP(range,i),1:i2);
% end

for i=1:length(range)%endLen
    magicMP(i,end-(range(i)):end) = nan;
end
% magicMP(magicMP<min(min(magicMP(:,100:end-500)))) = nan;

end

