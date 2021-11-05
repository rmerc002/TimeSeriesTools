dataR = rand(1,16326)<= 1/64;

dataCorrupted = diff([0,dataTermite]);
dataCorrupted(1) = 0;

basePairs = [-2,-1,0,1,2];
indices = randi([1,5],1,sum(dataR));

dataCorrupted(dataR) = basePairs(indices);

dataSumReduced = zeros(1,size(dataTermite,2));

for i=2:size(dataTermite,2)
dataSumReduced(i) = dataSumReduced(i-1) + dataCorrupted(i);
end

[magicMP, profileIndices,subLenSeries] = magicMatrixProfile(dataSumReduced);
visualizeMMP(dataSumReduced, magicMP, profileIndices, subLenSeries);