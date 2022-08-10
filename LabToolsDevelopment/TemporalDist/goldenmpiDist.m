function mp = goldenmpiDist(ts, mpi, mm)
    mp = nan(length(mpi),1);
    for ii = 1:length(mpi)
        if isnan(mpi(ii)) || mpi(ii) > length(ts)-mm+1
            continue;
        end
        startIndex = ii;
        endIndex = startIndex + mm - 1;
        ss1 = ts(startIndex:endIndex);
        tempMean = mean(ss1,'omitnan');
        tempStd = std(ss1,'omitnan');
        ss1 = (ss1-tempMean)/tempStd;

        startIndex = mpi(ii);
        endIndex = startIndex + mm - 1;
        ss2 = ts(startIndex:endIndex);
        tempMean = mean(ss2,'omitnan');
        tempStd = std(ss2,'omitnan');
        ss2 = (ss2-tempMean)/tempStd;

        mp(ii) = norm(ss1-ss2);
    end
end