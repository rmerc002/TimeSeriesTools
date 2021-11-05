function  [magicMP, profileIndices,scaleSeries, subLengthSeries] = magicMatrixProfileAB_ScaleDetrend(dataA,dataB, subLen, startScale, endScale, numScales, selfJoin)
if size(dataA,1) > 1
    dataA = dataA';
end
dataA = normalize(dataA);


%%% Determine subLengths to iterate over
scaleSeries = [];
% startScale = 5;
% endScale = ceil(length(dataA)/20);
% endScale = 3500;
% if length(dataA) < 8000
%     endScale = 600;
% end
% scaleSeries = startScale:endScale;
scaleSeries = getSubLenSeries(startScale, endScale, numScales);
subLengthSeries = subLen*ones(size(scaleSeries));
% scaleSeries = 10:10:800;
% scaleSeries = 1:50;

magicMP = zeros(length(scaleSeries),length(dataA)); %preallocate
profileIndices = zeros(length(scaleSeries),length(dataA)); %preallocate
for rangeIndex = 1:length(scaleSeries)
    scale = scaleSeries(rangeIndex)


%     [mpa, mpb, mpia, mpib] = mpx_ABBA(dataA',dataA',subLen);
%     matrixProfile = mpa;
%     profileIndex = mpia;
%     [matrixProfile, profileIndex] = prescrimp(dataA',dataA', subLen, 1,ceil(length(dataA)/10));

    dataBScaled = interp1(dataB,linspace(1,length(dataB),length(dataB)/scale)); 
    [matrixProfile,~,profileIndex,~] = abbaJoindetrend(dataA',dataBScaled',subLen, selfJoin);
%     profileIndices = ceil(profileIndices.*scale);
%     profileIndices = min(profileIndices, length(dataB));

    magicMP(rangeIndex,1:length(matrixProfile)) = matrixProfile;
    profileIndices(rangeIndex,1:length(profileIndex)) = profileIndex;
    
end

%%% remove values that are zero at the end of the matrix profile due to
%%% subLength
magicMP(:,end-subLen+1:end) = nan;

% magicMP(magicMP<min(min(magicMP(:,100:end-500)))) = nan;

profileIndices = ceil(profileIndices.*scaleSeries');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   HELPER FUNCTIONS   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function scaleSeries = getSubLenSeries(startScale, endScale, numScales)
%%% Purpose: Reduce space. matrix profile distances are pretty similar from
%%% from one subLen to the next. By trial and error, adding the square root
%%% of the current subLen seems to produce a good distribution.
%     scaleSeries = [];
%     index = startScale;
%     
%     while index < endScale
%         scaleSeries = [scaleSeries,index];
%         index = index + max(1,ceil(sqrt(index)));
%     end
%     
%     scaleSeries = [scaleSeries,endScale];

powerMin = log10(startScale);
powerMax = log10(endScale);
powerStep = (powerMax-powerMin)/numScales;
powers = powerMin:powerStep:powerMax;
scaleSeries = unique(power(10,powers));

end
