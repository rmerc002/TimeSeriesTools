function  [rMagicMP, rProfileIndices,lMagicMP, lProfileIndices,subLenSeries] = magicMatrixProfile_LeftRight(data, startLength, endLength, numLengths)
if size(data,1) > 1
    data = data';
end
data = normalize(data);

%add noise in order to increase contrast where the signal is low.
dataNoise = data;%amplifyDissimilarities(data);
% dataNoise = amplifyDissimilarities(data);


%%% Determine subLengths to iterate over
subLenSeries = [];
% startLength = 3;
% endLength = 5000;
% % endLength = ceil(length(data)/20);
% if length(data) < 8000
%     endLength = 600;
% end
% subLenSeries = startLength:endLength;
subLenSeries = getSubLenSeries(startLength, endLength, numLengths);
% subLenSeries = 10:10:800;
% subLenSeries = 1:50;

rMagicMP = zeros(length(subLenSeries),length(dataNoise)); %preallocate
rProfileIndices = zeros(length(subLenSeries),length(dataNoise)); %preallocate
lMagicMP = zeros(length(subLenSeries),length(dataNoise)); %preallocate
lProfileIndices = zeros(length(subLenSeries),length(dataNoise)); %preallocate
for rangeIndex = 1:length(subLenSeries)
    subLen = subLenSeries(rangeIndex)

    [rMatrixProfile, lMatrixProfile, rProfileIndex, lProfileIndex] = mpxLeftRight(dataNoise', ceil(subLen/2),subLen);

    rMagicMP(rangeIndex,1:length(rMatrixProfile)) = rMatrixProfile;
    rProfileIndices(rangeIndex,1:length(rProfileIndex)) = rProfileIndex;
    lMagicMP(rangeIndex,1:length(lMatrixProfile)) = lMatrixProfile;
    lProfileIndices(rangeIndex,1:length(lProfileIndex)) = lProfileIndex;
    
end

%%% remove values that are zero at the end of the matrix profile due to
%%% subLenth
for i=1:length(subLenSeries)%endLength
    rMagicMP(i,end-(subLenSeries(i)):end) = nan;
    lMagicMP(i,end-(subLenSeries(i)):end) = nan;
end
% magicMP(magicMP<min(min(magicMP(:,100:end-500)))) = nan;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   HELPER FUNCTIONS   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    threshold = (max(matrixProfile)-min(matrixProfile))/2 + min(matrixProfile)
    indices = matrixProfile < threshold;
    sum(indices)
    indices = [false(ceil(subLen/2),1);indices(1:end-ceil(subLen/2))];%matrix profile indices are offset
    dataRange = max(data(indices))-min(data(indices))
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

function subLenSeries = getSubLenSeries(startLength, endLength, numLengths)
%%% Purpose: Reduce space. matrix profile distances are pretty similar from
%%% from one subLen to the next. By trial and error, adding the square root
%%% of the current subLen seems to produce a good distribution.
%     subLenSeries = [];
%     index = startLength;
%     
%     while index < endLength
%         subLenSeries = [subLenSeries,index];
%         index = index + max(1,ceil(sqrt(index)));
%     end
%     
%     subLenSeries = [subLenSeries,endLength];

powerMin = log10(startLength);
powerMax = log10(endLength);
powerStep = (powerMax-powerMin)/numLengths;
powers = powerMin:powerStep:powerMax;
subLenSeries = unique(ceil(power(10,powers)));

end
