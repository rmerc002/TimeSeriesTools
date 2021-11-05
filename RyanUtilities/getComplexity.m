function complexity = getComplexity(ts)
    n = length(ts);
    complexity = sum(abs(diff(ts)))/n;
end