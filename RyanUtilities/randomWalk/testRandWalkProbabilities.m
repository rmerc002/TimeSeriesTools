tslengths = getExpDistributedSeries(2,100000,30);
numIterations = 1000;
numHits = zeros(1,length(tslengths));


for lengthIndex = 1:length(tslengths)
    tslength = tslengths(lengthIndex);
    fprintf("testing length %d, %d of %d\n", tslength, lengthIndex, length(tslengths));
    for iterIndex = 1:numIterations
       
       steps = randi([-1,1],2,tslength);
       ts1 = cumsum(steps(1,:));
       ts2 = cumsum(steps(2,:));
       ts2 = ts2(end:-1:1);
       if sum(abs(ts1 - ts2) <= 1)
           numHits(lengthIndex) = numHits(lengthIndex) + 1;
       end
    end
end
   
figure;
scatter(tslengths, numHits/numIterations,'filled');
set(gca,'xscale','log');
xlabel("Time Series Length");
ylabel("Probability of Intersection (1000 samples)");
title("Probability of Opposed Random Walk Intersection");
