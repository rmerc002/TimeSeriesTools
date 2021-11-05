function [mp,mpi] = mpchunk(data1, w,mpref)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if size(data1,1) > size(data1,2)
data1 = data1';    
end

data2 = data1;

hw = ceil(w/2); %half subsequence size, used for exclusion zone

%initialize output
mp = nan(1,length(data1));
mpi = nan(1,length(data1));

k = 5; %forshadowing knn mpchunk
chunkLength = w*k;
numSubsequencesPerChunk = chunkLength - (w-1);  % -(w-1) because w-1 of the last subsequence needs to be included in the next chunk
numChunks = ceil(size(data1,2)/(numSubsequencesPerChunk));
figure;
drawnow;
for chunkIndex1 = 1:numChunks-1
    startIndex = 1 + (chunkIndex1-1)*numSubsequencesPerChunk;
    endIndex = min(length(data1),startIndex + chunkLength -1);
%     disp(startIndex);
    for chunkIndex2 = chunkIndex1:numChunks-1
        startIndex2 = 1 + (chunkIndex2-1)*numSubsequencesPerChunk;
        endIndex2 = min(size(data2,2),startIndex2 + chunkLength -1);
%         disp(startIndex2);
%         if endIndex < w || endIndex2 < w
%             continue;
%         end
        if chunkIndex1 == chunkIndex2
            
            [temp_mp,temp_mpi] = mpx(data1(1,startIndex:endIndex)',ceil(w/2), w);
            mp(1,startIndex:startIndex + numSubsequencesPerChunk-1) = min(mp(1,startIndex:startIndex + numSubsequencesPerChunk-1), temp_mp(1:numSubsequencesPerChunk)','omitnan');
        else
            [mpa, mpb, mpia, mpib] = mpx_ABBA(data1(1,startIndex:endIndex)', data2(1,startIndex2:endIndex2)', w);
            
            mp(1,startIndex:startIndex + numSubsequencesPerChunk-1) = min(mp(1,startIndex:startIndex + numSubsequencesPerChunk-1) , mpa(1:numSubsequencesPerChunk)','omitnan');
%             mpi(1,startIndex:startIndex + numSubsequencesPerChunk-1) = min(mpi(1,startIndex:startIndex + numSubsequencesPerChunk-1) , mpia(1:numSubsequencesPerChunk),2);
            
           
           mp(1,startIndex2:startIndex2 + numSubsequencesPerChunk-1) = min(mp(1,startIndex2:startIndex2 + numSubsequencesPerChunk-1) , mpb(1:numSubsequencesPerChunk)','omitnan');
%            mpi(1,startIndex2:startIndex2 + numSubsequencesPerChunk-1) = min(mpi(1,startIndex2:startIndex2 + numSubsequencesPerChunk-1) , mpib(1:numSubsequencesPerChunk),2);
        end
        prompt = sprintf("startIndex1: %d, startIndex2: %d", startIndex, startIndex2);
        plot(mp(1:length(mpref))-mpref);
        drawnow;
        x = input(prompt);
        
    end
end

end

