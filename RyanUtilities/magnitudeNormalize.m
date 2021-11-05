function ts = magnitudeNormalize(ts)
    [~, sortedIndices] = sort(ts);
    ts(sortedIndices) = linspace(0,1,length(ts));
end