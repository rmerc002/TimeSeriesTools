
data = binaryLabels;
queryIndices = [];

index = 2;
while true
    if data(index) == 1 && data(index-1) == 0
       queryIndices = [queryIndices,index]; 
       index = index + 1000;
    else
        index = index+1;
    end
    if index >= length(data)
       break; 
    end
end


for i = 1:length(queryIndices)
   qi = queryIndices(i);
   fprintf("Iteration: %d, index: %d\n", i, qi);
   motif1 = zscore(motionDataTrain(qi:qi+subLen-1));
   DTWExperiment_Test_DTW_ED;
   w = waitforbuttonpress;
end
