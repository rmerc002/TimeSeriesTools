figure;
tiledlayout(2,1);

ax1 = nexttile;
plot(S1Drill{:,250}==407521);

ax2 = nexttile;
plot(S1Drill{:,135});

linkaxes([ax1 ax2],'x');


activityMask = S1Drill{:,250}==407521; %Drinking while standing and seated

%%% record start and end indices for drinking
i = 2;
startEndIndices = [];%[0,0];
while i <= length(activityMask)
    if activityMask(i) == 1 && activityMask(i-1) == 0
        startEndIndices(end+1,:) = [0,0];
        startEndIndices(end,1) = i;
       

    elseif activityMask(i) == 0 && activityMask(i-1) == 1
        startEndIndices(end,2) = i;
        
    end
    i = i + 1;
end


wholeLength = size(S1Drill,1);
halfLength = 28300;%ceil(wholeLength/2);
wholeTS = S1Drill{:,69};
negativeTSRaw = wholeTS(1:halfLength);
negativeTS = negativeTSRaw;
positiveTSRaw = wholeTS(halfLength+1:end);
positiveTS = positiveTSRaw;

%%%fill negative time series desired activity with random walk
for i = 1:size(startEndIndices,1)
    startIndex = startEndIndices(i,1);
    endIndex = startEndIndices(i,2);
    if startIndex > halfLength
        break;
    end
    signalLength = endIndex-startIndex+1;
    randWalkSignal = zscore(getRandWalk(signalLength));
    stdMatch = std(negativeTSRaw(startIndex:endIndex));
    randWalkSignal = randWalkSignal*stdMatch;
    negativeTS(startIndex:endIndex) = randWalkSignal;
    
end

%%%fill 2nd drink activity with random walk in positiveTS
%%%TODO: the issue is the signal doesn't start and end at 0
for i = 1:2:size(startEndIndices,1)
    startIndex = startEndIndices(i,1) - halfLength+1;
    endIndex = startEndIndices(i,2) - halfLength+1;
    if startIndex < 0
        continue;
    end
    signalLength = endIndex-startIndex+1;
    randWalkSignal = zscore(getRandWalk(signalLength));
    stdMatch = std(positiveTSRaw(startIndex:endIndex));
    randWalkSignal = randWalkSignal*stdMatch;
    positiveTS(startIndex:endIndex) = randWalkSignal;
end

figure; 
plot(positiveTS);
hold on;
plot(positiveTSRaw);
hold off;

groundTruthIndices = startEndIndices(21:2:end,1) - halfLength + 1;

subLength = signalLength;

[plato] = PLATO(positiveTS, negativeTS, subLength, groundTruthIndices);

%%%%%%%%%%%%%%%
%%% TESTING %%%
%%%%%%%%%%%%%%%
testTS = S1ADL4{:,69};
labelTS = S1ADL4{:,250}==407521;
testTS(isnan(testTS)) = nanmean(testTS);
distanceProfile = MASS_V2(testTS, plato);
distanceProfile = real(distanceProfile);

K = 10;
[bestIndices] = KLowestDistanceIndices_V05(distanceProfile, subLength, K);

figure;
tiledlayout(3,1);

ax1 = nexttile;
plot(testTS);

ax2 = nexttile;
plot(distanceProfile);
hold on;
scatter(bestIndices,ones(size(bestIndices)));
hold off;
ylim([-0.1, 1.1]);

ax3 = nexttile;
plot(labelTS);

linkaxes([ax1, ax2, ax3],'x');
