function  [magicMP, profileIndices,scaleSeries, subLengthSeries] = magicMatrixProfileAB_SmoothDetrend(dataA,dataB, subLen, startValue, endValue, numValues, selfJoin)
if size(dataA,1) > 1
    dataA = dataA';
end
dataA = normalize(dataA);


%%% Determine subLengths to iterate over
scaleSeries = [];
% startValue = 5;
% endValue = ceil(length(dataA)/20);
% endValue = 3500;
% if length(dataA) < 8000
%     endValue = 600;
% end
% scaleSeries = startValue:endValue;
scaleSeries = getSubLenSeries(startValue, endValue, numValues)
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

      dataASmooth = smoothdata(dataA, "movmedian",scale)' + rand(size(dataA'))*0.001; 
      dataBSmooth = smoothdata(dataB, "movmedian",scale)'+ rand(size(dataB'))*0.001;
      
      [matrixProfile,~,profileIndex,~] = abbaJoindetrend(dataASmooth,dataBSmooth,subLen, selfJoin);
%       [matrixProfile,~,profileIndex,~] = mpx_ABBA(dataASmooth,dataBSmooth,subLen);

      %     profileIndices = ceil(profileIndices.*scale);
%     profileIndices = min(profileIndices, length(dataB));

    magicMP(rangeIndex,1:length(matrixProfile)) = matrixProfile;
    profileIndices(rangeIndex,1:length(profileIndex)) = profileIndex;
    
end

%%% remove values that are zero at the end of the matrix profile due to
%%% subLength
magicMP(:,end-subLen+1:end) = nan;

% magicMP(magicMP<min(min(magicMP(:,100:end-500)))) = nan;

% profileIndices = ceil(profileIndices.*scaleSeries');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   HELPER FUNCTIONS   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function scaleSeries = getSubLenSeries(startValue, endValue, numValues)
%%% Purpose: Reduce space. matrix profile distances are pretty similar from
%%% from one subLen to the next. By trial and error, adding the square root
%%% of the current subLen seems to produce a good distribution.
%     scaleSeries = [];
%     index = startValue;
%     
%     while index < endValue
%         scaleSeries = [scaleSeries,index];
%         index = index + max(1,ceil(sqrt(index)));
%     end
%     
%     scaleSeries = [scaleSeries,endValue];

powerMin = log10(startValue);
powerMax = log10(endValue);
powerStep = (powerMax-powerMin)/numValues;
powers = powerMin:powerStep:powerMax;
scaleSeries = unique(ceil(power(10,powers)));

end
