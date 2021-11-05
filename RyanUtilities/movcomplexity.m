function complexity = movcomplexity(ts, m)
    ts = reshape(ts, length(ts),1);
    complexity = nan(length(ts),1);
    tsDiff = [0;diff(ts)];
    for ii = 1:length(ts)
        startIndex = max(1, ii-m+1);
        endIndex = ii;
        tempSubsequence = tsDiff(startIndex:endIndex);
        if sum(isnan(tempSubsequence)) == 0
            tempComplexity = sum(abs(tempSubsequence));
            complexity(ii) = tempComplexity;
        end
    end
end