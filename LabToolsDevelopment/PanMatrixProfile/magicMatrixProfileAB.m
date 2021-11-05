function  [magicMP, profileIndices,subLenSeries] = magicMatrixProfileAB(dataA,dataB, startLength, endLength, numLengths)
if size(dataA,1) > 1
    dataA = dataA';
end
dataA = normalize(dataA);


%%% Determine subLengths to iterate over
subLenSeries = [];
% startLength = 5;
% endLength = ceil(length(dataA)/20);
% endLength = 3500;
% if length(dataA) < 8000
%     endLength = 600;
% end
% subLenSeries = startLength:endLength;
subLenSeries = getSubLenSeries(startLength, endLength, numLengths);
% subLenSeries = 10:10:800;
% subLenSeries = 1:50;

magicMP = zeros(length(subLenSeries),length(dataA)); %preallocate
profileIndices = zeros(length(subLenSeries),length(dataA)); %preallocate
for rangeIndex = 1:length(subLenSeries)
    subLen = subLenSeries(rangeIndex)


%     [mpa, mpb, mpia, mpib] = mpx_ABBA(dataA',dataA',subLen);
%     matrixProfile = mpa;
%     profileIndex = mpia;
%     [matrixProfile, profileIndex] = prescrimp(dataA',dataA', subLen, 1,ceil(length(dataA)/10));


    [matrixProfile,~,profileIndex,~] = abbaJoindetrend(dataA',dataB',subLen, false);
%     [matrixProfile,profileIndex] = mpx(dataA',ceil(subLen/2),subLen);

    magicMP(rangeIndex,1:length(matrixProfile)) = matrixProfile;
    profileIndices(rangeIndex,1:length(profileIndex)) = profileIndex;
    
end

%%% remove values that are zero at the end of the matrix profile due to
%%% subLenth
for i=1:length(subLenSeries)%endLength
    magicMP(i,end-(subLenSeries(i)):end) = nan;
end
% magicMP(magicMP<min(min(magicMP(:,100:end-500)))) = nan;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   HELPER FUNCTIONS   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
