function subLenSeries = getSubLenSeries(startLength, endLength, numLengths)
%%% Purpose: Reduce space. matrix profile distances are pretty similar from
%%% from one subLen to the next. By trial and error, adding the square root
%%% of the current subLen seems to produce a good distribution.
    powerMin = log10(startLength);
    powerMax = log10(endLength);
    powerStep = (powerMax-powerMin)/numLengths;
    powers = powerMin:powerStep:powerMax;
    subLenSeries = unique(ceil(power(10,powers)));
end

