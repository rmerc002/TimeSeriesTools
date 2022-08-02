function mup = MutationProfile(ts, mm, rr)

%%% requires r > mm
    mup = zeros(length(ts)-rr+1,1);
    updateIndices = getExpDistributedSeries(1,length(ts), 30);
    for ii = 1:length(ts)-rr+1
        if min(abs(updateIndices-ii)) == 0
            fprintf("Index %d\n", ii);
        end
        startIndex = ii;
        endIndex = ii + mm - 1;
        query = ts(startIndex:endIndex);

        startIndex = endIndex+1;
        endIndex = startIndex + rr - mm - 1;

        tempTS = ts(startIndex:endIndex);

        mp = MASS_V2_nan(tempTS, query);
        mup(ii) = min(mp, [], 'omitnan')/sqrt(2*mm);
    end
    mup = movmin(mup, mm);

    figure; 
    tiledlayout(2,1);
    ax1 = nexttile();
    plot(ts);
    
    ax2 = nexttile();
    plot(mup);
    ylim([0,1]);

    linkaxes([ax1, ax2], 'x');
end