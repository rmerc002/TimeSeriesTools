function [mp, mpi] = mp_dtw_rand(data1, subLen, maxwarp)
mp = nan(1,length(data1));
mpi = nan(1,length(data1));

halfSubLen = ceil(subLen/2);
for iter = 1:1215000%0
    index1 = randi([1,length(data1)-subLen]);
    index2 = randi([1,length(data1)-subLen]);
%     fprintf("index1: %d, index2: %d\n", index1, index2);
    if mod(iter,100000) == 0
        fprintf("index1: %d\n",iter);
    end
    ss1 = data1(index1:index1+subLen-1);
    ss1 = zscore(ss1);
%     for index2 = index1+halfSubLen:length(data1) - subLen
%         if index2 >= length(data1)-subLen - 5
%            continue; 
%         end
        ss2 = data1(index2:min(index2+subLen-1, length(data1)-subLen));
        if length(ss2) < subLen || length(ss1) < subLen || abs(index1-index2) <= halfSubLen
            continue
        end
        ss2 = zscore(ss2);
%         fprintf("length ss1: %d, length ss2: %d\n",length(ss1), length(ss2));
        tempDist = dtw_upd(ss1,ss2,maxwarp);
        if isnan(mp(index1)) || tempDist < mp(index1)
            mp(1,index1) = tempDist;
            mpi(1,index1) = index2;
        end
        
        if isnan(mp(index2)) || tempDist < mp(index2)
            mp(1,index2) = tempDist;
            mpi(1,index2) = index1;
        end
%     end
end

end
