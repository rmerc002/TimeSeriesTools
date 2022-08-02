function TemporalDistQuery(ts, TDProfile, maxIndex)% [~,maxIndex] = max(NNSProfile);
% maxIndex = 23799;
% targetIndices = TDProfile == maxIndex;

targetIndices = (TDProfile < maxIndex + 20).*(TDProfile > maxIndex - 20);



figure;
hold on;
plot(ts);
X = 1:length(TDProfile);
Y = zeros(length(TDProfile),1);
scatter(X(targetIndices==1), Y(targetIndices == 1),'filled');
hold off;
set(gca, 'TickDir','out');
box off;


% mm = 24*7;
% queryIndex = 55514;
% queryNNIndex = mpi(queryIndex);
% query = ts(queryIndex:queryIndex + mm - 1);
% queryNN = ts(queryNNIndex:queryNNIndex + mm - 1);
% 
% figure;
% hold on;
% plot(zscore(query));
% plot(zscore(queryNN));
% hold off;
% set(gca, 'TickDir','out');
% box off;
