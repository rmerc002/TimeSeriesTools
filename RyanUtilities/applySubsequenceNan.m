function profileNan = applySubsequenceNan(ts, profile, mm)
    nanIndices = isnan(ts);
    profileNan = profile;
    profileNan(nanIndices(1:length(profile))) = nan;
    for ii = 2:length(profile)
        if isnan(ts(ii)) && ~isnan(ts(ii-1))
            startIndex = max(1,ii-mm+1);
            endIndex = ii-1;
            profileNan(startIndex:endIndex) = nan;
        end
    end
end