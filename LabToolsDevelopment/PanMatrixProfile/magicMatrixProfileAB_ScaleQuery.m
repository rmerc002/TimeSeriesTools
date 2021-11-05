function  [dataB, magicMP, profileIndices,scaleSeries, subLengthSeries] = magicMatrixProfileAB_ScaleQuery(dataA, queries, startScale, endScale, numScales)
if size(dataA,1) > 1
    dataA = dataA';
end

%%% Determine subLengths to iterate over
scaleSeries = getSubLenSeries(startScale, endScale, numScales);
subLengthSeries = zeros(1,length(scaleSeries));

dataB = [];

magicMP = ones(length(scaleSeries),length(dataA)); %preallocate
profileIndices = zeros(length(scaleSeries),length(dataA)); %preallocate
for queryIndex = 1:length(queries)
    dataBQuery = queries{queryIndex};
    
    %dataB used for Profile Indices
    dataBIndex = length(dataB)+1;
    dataB = [dataB, dataBQuery];
    
    for rangeIndex = 1:length(scaleSeries)
        scale = scaleSeries(rangeIndex)

        subLength = ceil(scale*length(dataBQuery));
        subLengthSeries(rangeIndex) = subLength; %TODO: this doesn't need to be rewritten every loop of scales 
        
        dataBQueryInterp = interp1(dataBQuery,linspace(1,length(dataBQuery),subLength)); 
        
        [tempMP,~] = mpx_AB(dataA',dataBQueryInterp',subLength);
        tempMatrixProfile = tempMP';

        bestIndices = tempMatrixProfile < magicMP(rangeIndex,1:length(tempMatrixProfile));
        magicMP(rangeIndex,bestIndices) = tempMatrixProfile(bestIndices);
        profileIndices(rangeIndex,bestIndices) = dataBIndex;
        
        magicMP(rangeIndex,end-subLength+1:end) = nan;
    end
    
end

%%% remove values that are zero at the end of the matrix profile due to
%%% subLenth


% magicMP(magicMP<min(min(magicMP(:,100:end-500)))) = nan;

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
