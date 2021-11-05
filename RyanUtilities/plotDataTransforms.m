function plotDataTransforms(ts)

tsNonNan = ts;
tsNonNan(isnan(ts)) = 0;
tsMean = nanmean(ts);
tsSTD = nanstd(ts);
tsNonNan = (tsNonNan-tsMean)/tsSTD;
tsTransforms = {};
tsTransforms{end+1} = {ts,"Raw"};
tsTransforms{end+1} = {log(ts),"Log"};
tsTransforms{end+1} = {diff(ts),"Diff"};
tsTransforms{end+1} = {diff(diff(ts)),"Diff2"};
tsTransforms{end+1} = {cumsum(tsNonNan),"CumSum(ZScore)"};
tsTransforms{end+1} = {cumsum(cumsum(tsNonNan)),"CumSum2(ZScore)"};
tsNorm = (ts-tsMean)/tsSTD;
tsTransforms{end+1} = {tsNorm.^2,"ZScore^2"};
tsTransforms{end+1} = {tsNorm.^3,"ZScore^3"};

figure;
inset = 0.9;
hold on;
for transformIndex = 1:length(tsTransforms)
   plotIndex = transformIndex-1;
   name = tsTransforms{transformIndex}{2};
   tempTS = tsTransforms{transformIndex}{1};
   tempMin = nanmin(tempTS);
   tempMax = nanmax(tempTS);
   tempRange = tempMax - tempMin;
   tempPlot = -plotIndex + inset*(tempTS - tempMin)/tempRange;
   
   plot(tempPlot);
   text(0, -plotIndex+0.8,name,'FontSize',12);
end
plot([0,length(ts)],[-0.1, -0.1], '--','Color',[0.5,0.5,0.5]);
hold off;

titleFormat = sprintf("Time Series Transforms");
title(titleFormat);
xlabel("Time Index");
set(gca,'xtick',[1,length(ts)],'ytick',[], 'TickDir','out');
box off;
       
end