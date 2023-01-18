function [tsSeg] = segmentConcatenation(ts)
    tsSeg = {};
    startIndex = 1;
    for ii = 2:length(ts)
        if isnan(ts(ii))
            if startIndex == ii
                startIndex = ii+1;
                continue;
            end
            tsSeg{end+1,1} = ts(startIndex:ii-1);
            tsSeg{end,2} = ii-1;
            startIndex = ii+1;
        end
    end
end