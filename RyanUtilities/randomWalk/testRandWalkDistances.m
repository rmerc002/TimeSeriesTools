tslengths = getExpDistributedSeries(2,100000,30);
tsLength = 100000;
numIterations = 100;
distances = zeros(1,length(tsLength));



for iterIndex = 1:numIterations
   steps = randi([-1,1],2,tsLength);
   ts1 = cumsum(steps(1,:));
   ts2 = cumsum(steps(2,:));
   ts2 = ts2(end:-1:1);

   distances = distances + abs(ts1 - ts2);
end
   
figure;
scatter(1:tsLength, distances/numIterations,'filled');
% set(gca,'xscale','log');
xlabel("Time Series Index");
ylabel("Mean Distance Between Random Walks");
title("How far are points on average?");
