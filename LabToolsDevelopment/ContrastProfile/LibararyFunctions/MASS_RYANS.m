function [dp] = MASS_RYANS(ts, query)
subLength = length(query);
dp = inf(length(ts) - length(query) + 1);

for i1 = 1:length(dp)
   dp(i1) = norm(zscore(ts(i1:i1+subLength-1)) - zscore(query)); 
end
end