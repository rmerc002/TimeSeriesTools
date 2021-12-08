function subsequences = getRandomSubsequences(ts, m, KNN)
    if isempty(KNN)
        KNN = 10;
    end
    subsequences = nan(KNN,m);
    for ii = 1:KNN
        subsequence = nan(1,200);
        while(sum(isnan(subsequence))> 0)
            startIndex = max(1, floor(rand(1)*length(ts)-m+1))
            endIndex = startIndex + m - 1
            subsequence = ts(startIndex:endIndex);
        end
        subsequences(ii,:) = subsequence;
    end
end

