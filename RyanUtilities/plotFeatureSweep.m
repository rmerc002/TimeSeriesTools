function plotFeatureSweep(ts, startLen, endLen,numSteps)
    %%% Generate sets of plots showing moving window of features
    %%% Features:
    %%%   Mean
    %%%   STD
    %%%   Complexity
    %%%   Max
    %%%   Min
    
    featureFuncs = {};
    featureFuncs{1} = {@movmean,"Moving Mean"};
    featureFuncs{2} = {@movstd, "Moving STD"};
    featureFuncs{3} = {@movcomplexity, "Moving Complexity"};
    featureFuncs{4} = {@movmax, "Moving Max"};
    featureFuncs{5} = {@movmin, "Moving Min"};
    
    subLenSeries = getExpDistributedSeries(startLen, endLen, numSteps);
    inset = 0.9;
    for ffi = 1:length(featureFuncs)
       figure; 
       hold on;
       thisFunc = featureFuncs{ffi}{1};
       funcName = featureFuncs{ffi}{2};
       for mIndex = 0:length(subLenSeries)
          plotIndex = mIndex;
          if mIndex == 0
              tempTS = ts;
          else
              m = subLenSeries(mIndex);
              tempTS = thisFunc(ts, m);
              tempTS = omitNaN(ts, tempTS, m); 
          end
          tempMin = nanmin(tempTS);
          tempMax = nanmax(tempTS);
          tempRange = tempMax - tempMin;
          tempPlot = -plotIndex + inset*(tempTS-tempMin)/tempRange;
          
          plot(tempPlot);
       end
       plot([0,length(ts)], [-0.1,-0.1], '--', 'Color', [0.5, 0.5, 0.5])
       hold off;
       lengthString = string(subLenSeries(1));
       for li = 2:length(subLenSeries)
          m = subLenSeries(li);
          lengthString = lengthString + ", " + string(m); 
       end
       titleFormat = sprintf("%s: m = 0, %s",funcName, lengthString);
       title(titleFormat);
       xlabel("Time Index");
       ylabel("Increasing Window Length m");
       set(gca,'xtick',[1,length(ts)],'ytick',[], 'TickDir','out');
       box off;
    end
end

function newTS = omitNaN(ts, newTS, m)
    for ii = m:length(ts)
        startIndex = ii-m+1;
        endIndex = ii;
        if sum(isnan(ts(startIndex:endIndex)))
           newTS(startIndex:endIndex) = nan;
        end
    end
end