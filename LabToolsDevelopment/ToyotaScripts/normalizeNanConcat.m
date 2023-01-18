function pos = normalizeNanConcat(pos)
    startIndex = 1;
    for ii = 2:length(pos)
        if isnan(pos(ii))
            pos(startIndex:ii-1) = zscore(pos(startIndex:ii-1));
            startIndex = ii+1;
        end
    end
end