function  [magicMP, profileIndices,scaleSeries, subLenSeries] = magicMatrixProfileAB_LengthScale(dataA,dataB, startLength, endLength, numLengths, startScale, endScale, numScales)
if size(dataA,1) > 1
    dataA = dataA';
end


%%% Determine subLengths to iterate over

subLengthSeries = getSubLenSeries(startLength, endLength, numLengths);

scaleSeries = getSubLenSeries(startScale, endScale, numScales);


magicMP = zeros(length(subLengthSeries),length(dataA)); %preallocate
profileIndices = zeros(length(subLengthSeries),length(dataA)); %preallocate

for subLengthIndex = 1:length(subLengthSeries)
    subLength = subLengthSeries(subLengthIndex)
    
    for scaleIndex = 1:length(scaleSeries)
        scale = scaleSeries(scaleIndex)


    %     [mpa, mpb, mpia, mpib] = mpx_ABBA(dataA',dataA',subLen);
    %     matrixProfile = mpa;
    %     profileIndex = mpia;
    %     [matrixProfile, profileIndex] = prescrimp(dataA',dataA', subLen, 1,ceil(length(dataA)/10));

        dataBScaled = interp1(dataB,linspace(1,length(dataB),length(dataB)/scale)); 
        
        [matrixProfile,profileIndex] = mpx_AB(dataA',dataBScaled',subLength);
        profileIndices = ceil(profileIndices.*scale);
        profileIndices = min(profileIndices, length(dataB));

        magicMP(rangeIndex,1:length(matrixProfile)) = matrixProfile;
        profileIndices(rangeIndex,1:length(profileIndex)) = profileIndex;

    end
end

%%% remove values that are zero at the end of the matrix profile due to
%%% subLenth
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
