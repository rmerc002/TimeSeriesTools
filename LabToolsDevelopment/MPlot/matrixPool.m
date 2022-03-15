function Y = matrixPool(X, poolsize, extrema)
    numBlocks = floor(size(X,1)/poolsize);
    Y = zeros(numBlocks);
    for rr = 1:numBlocks
        for cc = 1:numBlocks
            rrStart = poolsize*(rr-1) + 1;
            rrEnd = rrStart + poolsize - 1;
    
            ccStart = poolsize*(cc-1) + 1;
            ccEnd = ccStart + poolsize - 1;
            if extrema == "max"
                Y(rr,cc) = max(X(rrStart:rrEnd,ccStart:ccEnd),[],'all');
            else
                Y(rr,cc) = min(X(rrStart:rrEnd,ccStart:ccEnd),[],'all');
            end
        end
    end
end