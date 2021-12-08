function tsOut = mov01Norm(ts,m)
%%%It seems better to have a sliding window length smaller than pattern
%%%length
%%% This is probably because it reduces number of rising edges with many
%%% low amplitude regions
ts = reshape(ts,1,length(ts));
tsOut = nan(size(ts));
for i=1:length(ts)-m+1
    startIndex = i;
    endIndex = i+m-1;
%     if sum(isnan(ts(startIndex:endIndex)))
    tsOut(startIndex:endIndex) = sum([tsOut(startIndex:endIndex);normalize(ts(startIndex:endIndex),'range')],'omitnan');
end
tsOut = tsOut/m;
tsOut(isnan(ts)) = nan;
end
