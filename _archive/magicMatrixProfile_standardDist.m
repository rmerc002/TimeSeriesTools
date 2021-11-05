function  [magicMP, profileIndices,range,dataNoise] = magicMatrixProfile_standardDist(data)
data = normalize(data);

findSimilarities = true;
% findSimilarities = false;
%add noise in order to increase contrast where the signal is low.
subLen = 10;
[matrixProfile, ~, ~, ~, mainWindow] = interactiveMatrixProfileVer2_xlogx_sublenNorm(data, subLen, 0);
threshold = (max(matrixProfile)+min(matrixProfile))/4;
indices = matrixProfile < threshold;
indices = [false(ceil(subLen/2),1);indices(1:end-ceil(subLen/2))];%matrix profile indices are offset
% dataNoise(indices) = data(indices);
dataRange = max(data(indices))-min(data(indices));
if findSimilarities == true
    noise = rand(1,length(data));%randFunc(length(data));
else
    noise = ones(1,length(data));
    noise(1:2:length(data)) = zeros(1,length(1:2:length(data)));
end

noise = noise.*dataRange-dataRange/2;
noise = noise*1.5;
mpNorm01 = matrixProfile-min(matrixProfile);
mpNorm01 = mpNorm01./0.33;%max(mpNorm01);

if findSimilarities == false
    mpNorm01 = 1-mpNorm01;
end
    
% mpNorm01(mpNorm01>0.5) = 1.5-mpNorm01(mpNorm01>0.5);
mpNorm01 = [ones(1,ceil(subLen/2)), mpNorm01(1:end-ceil(subLen/2))'];
noise(1:length(mpNorm01)) = noise(1:length(mpNorm01)).*mpNorm01;
dataNoise = data+ noise;


range = [];
startLen = 10;
endLen = ceil(length(data)/20);
if length(data) < 8000
    endLen = 400;
end
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

%       [matrixProfile, profileIndex, ~, ~, mainWindow] = interactiveMatrixProfileVer2_sublenNorm(dataNoise, subLen, 0);
      [matrixProfile, profileIndex, ~, ~, mainWindow] = interactiveMatrixProfileVer2_xlogx_sublenNorm(dataNoise, subLen, 0);

      drawnow;
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

