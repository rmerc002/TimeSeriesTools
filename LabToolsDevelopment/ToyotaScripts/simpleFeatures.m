posAnomFilled = posAnom;
posAnomFilled(isnan(posAnomFilled)) = 0;

negAnomFilled = negAnom;
negAnomFilled(isnan(negAnomFilled)) = 0;

posBaseFilled = posBase;
posBaseFilled(isnan(posBaseFilled)) = 0;

negBaseFilled = negBase;
negBaseFilled(isnan(negBaseFilled)) = 0;

EDAnom  = sqrt(sum((posAnomFilled - negAnomFilled) .^ 2));
EDBase  = sqrt(sum((posBaseFilled - negBaseFilled) .^ 2));

posAnomVar = nanstd(posAnom)
negAnomVar = nanstd(negAnom)
posBaseVar = nanstd(posBase)
negBaseVar = nanstd(negBase)

posAnomMean = nanmean(posAnom)
negAnomMean = nanmean(negAnom)
posBaseMean = nanmean(posBase)
negBaseMean = nanmean(negBase)

posAnomComplexity = getComplexity(posAnomFilled)
negAnomComplexity = getComplexity(negAnomFilled)
posBaseComplexity = getComplexity(posBaseFilled)
negBaseComplexity = getComplexity(negBaseFilled)

%%% Need to make histograms of individual cars
pos = posBase;

featureVector = [];
startIndex = 1;
for ii = 2:length(pos)
    if isnan(pos(ii))
        featureVector(end+1) = mean(pos(startIndex:ii-1));
        startIndex = ii+1;
    end
end

figure;
histogram(featureVector);
title("posBase")