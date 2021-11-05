function  [PMPA, PMPB, PIA, PIB, subLenSeries] = magicMatrixProfileABBA_Strings(dataA,dataB, startLength, endLength, numLengths)
if size(dataA,1) > size(dataA,2)
    dataA = dataA';
end


%%% Determine subLengths to iterate over
subLenSeries = [];

subLenSeries = getSubLenSeries(startLength, endLength, numLengths);

PMPA = nan(length(subLenSeries),length(dataA)); %preallocate
PIA = zeros(length(subLenSeries),length(dataA)); %preallocate

PMPB = nan(length(subLenSeries),length(dataB)); %preallocate
PIB = zeros(length(subLenSeries),length(dataB)); %preallocate

for rangeIndex = 1:length(subLenSeries)
    subLen = subLenSeries(rangeIndex)


%     [mpa, mpb, mpia, mpib] = mpx_ABBA(dataA',dataA',subLen);
%     matrixProfile = mpa;
%     profileIndex = mpia;
%     [matrixProfile, profileIndex] = prescrimp(dataA',dataA', subLen, 1,ceil(length(dataA)/10));


    [mpa, mpb, pia, pib] = mpx_stringsAB(dataA', dataB', subLen);

    PMPA(rangeIndex,1:length(mpa)) = mpa;
    PIA(rangeIndex,1:length(pia)) = pia;
    PMPB(rangeIndex,1:length(mpb)) = mpb;
    PIB(rangeIndex,1:length(pib)) = pib;
    
    %     [matrixProfile,profileIndex] = mpx(dataA',ceil(subLen/2),subLen);

%     PMPA(rangeIndex,1:length(matrixProfile)) = matrixProfile;
%     pia(rangeIndex,1:length(profileIndex)) = profileIndex;
    
end

%%% remove values that are zero at the end of the matrix profile due to
%%% subLenth
for i=1:length(subLenSeries)%endLength
    PMPA(i,end-(subLenSeries(i)):end) = nan;
    PMPB(i,end-(subLenSeries(i)):end) = nan;
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
