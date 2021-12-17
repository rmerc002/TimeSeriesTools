function [sortScore, bestK] = sortScoreOnTwoClassSegmented_OR(posSamples, negSamples, platosPos)
    numPlatos = size(platosPos,1);
    sortScore = 0;
    bestK = 1;

    numSamples = size(posSamples, 1) + size(negSamples, 1);

    sampleKPlatoDistances = ones(numPlatos, numSamples);
    
    labels = [ones(1,size(posSamples,1)), zeros(1,size(negSamples,1))];
    samples = [posSamples; negSamples];

    m = size(platosPos,2);
    numSamples = length(labels);

    sampleDistances = zeros(1,numSamples);
    for sampleIndex = 1:numSamples
        sample = samples(sampleIndex,:);
        minDist = inf;
        for platoIndex = 1:numPlatos
            plato = platosPos(platoIndex,:);
            tempDist = min(real(MASS_V2(sample, plato)))/sqrt(2*m);
            tempDist = min(1, tempDist);
            sampleKPlatoDistances(platoIndex, sampleIndex) = tempDist;
        end
    end
    
    maxSortScore = 0;
    bestK = 1; %%%same as platoIndex
    for platoIndex = 1:numPlatos
        tempDistances = min(sampleKPlatoDistances(1:platoIndex,:),[],1);
        [~, sortedIndices] = sort(tempDistances);
        sortedLabels = labels(sortedIndices);
    
        distanceWeights = linspace(1,0,numSamples);
    
        weightedScores = sortedLabels.*distanceWeights;
        sortScore = mean(weightedScores);

        if sortScore > maxSortScore
            maxSortScore = sortScore;
            bestK = platoIndex;
        end
    end

    tempDistances = sampleKPlatoDistances(bestK,:);
    numPosSamples = size(posSamples,1);
    numNegSamples = size(negSamples,1);
    figure;
    hold on;
    scatter(tempDistances(1:numPosSamples), zeros(1,numPosSamples),'filled');
    scatter(tempDistances(numPosSamples+1:end), zeros(1,numNegSamples),'filled');
    hold off;
end
