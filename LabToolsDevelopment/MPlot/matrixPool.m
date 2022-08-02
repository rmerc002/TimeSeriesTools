function Y = matrixPool(X, poolsize, extrema)
    numRowBlocks = floor(size(X,1)/poolsize);
    numColBlocks = floor(size(X,2)/poolsize);
    Y = zeros(numRowBlocks, numColBlocks);
    for rr = 1:numRowBlocks
        for cc = 1:numColBlocks
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