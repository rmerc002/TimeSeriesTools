function [TD, TDProfile] = TemporalDist(ts, mm, forcePlot, mp, mpi)% ts = data(:,5);
    %%%TODO: make a class
    %%%     search for a single TemporalDist in the ts
    %%%     search for a range around the TemporalDist
    %%%     visualize a query for a given Temporal Dist
    %%%     interactive temporalDist highlighting
    %%%     launch app for interactivity
    %%%     selection between weighted vs unweighted frequency

    exclusion_mm = ceil(mm/2);

    if nargin < 3
        forcePlot = false;
    end
    if nargin < 4
        [mp, mpi] = mpx_v3(ts, exclusion_mm, mm, false);
    end

    ts = reshape(ts,length(ts),1);
 
    mpOriginal = reshape(mp, length(mp), 1);
    mpi = reshape(mpi, length(mpi), 1);
    
    indices = 1:length(mp);
    indices = indices';
    TDProfile = abs(mpi-indices); %%% Nearest Neighbor Spatial Distances
    
    %%% normalize by dividing by noise equivalent ED, ignore anti-correlated
    mp = mpOriginal/sqrt(2*mm);
    mp = min(1, mp);
    mp = 1-mp;
    
    %%% Nearest Neighbor Spatial Profile
    TD = zeros(length(mp),1);
    for ii = 1:length(mp)
        nnsd = TDProfile(ii);
        if nnsd <= 0 || isnan(nnsd)
            continue;
        end
        nnsd = ceil(nnsd);
        TD(nnsd) = TD(nnsd) + 1;%mp(ii);
    end
%     TD = TD/length(mp); %%% time series length normalize
    distributionNorm = linspace(2,0,length(TD))';
    
    TD = TD./distributionNorm;
    %%%delete last 20%
    startIndex = ceil(0.8*length(TD));
    TD(startIndex:end) = [];
    TD = TD/mean(TD,'omitnan');

    TDhat = mean(TD,'omitnan');
    RMSE = sqrt(mean((TD - TDhat).^2,'omitnan'));

    weights = linspace(1,0,length(TD));
    earlyScore = (weights*TD)/length(TD);
    earlyScore = 2*(earlyScore- 0.5); %%% 0.5 if if there is an even distribution
    
    %%%%%%%%%%%%%%%%%
    %%%   Plots   %%%
    %%%%%%%%%%%%%%%%%
    if forcePlot == true
        figure;
        tiledlayout(2,1);
        
        ax1 = nexttile();
        plot(TD);
        set(gca, 'TickDir','out');
        box off;
        title("TemporalDist");
        

        TDLength = length(TD);
        numBins = 30;
        binIndices = round(linspace(1,TDLength,numBins+1));
        TDBinned = zeros(numBins,1);
        for ii = 1:numBins
            startIndex = binIndices(ii);
            endIndex = binIndices(ii+1)-1;
            TDBinned(ii) = mean(TD(startIndex:endIndex));
        end
        barX = round(movmean(binIndices,2));
        barX(1) = [];
       
        ax2 = nexttile();
        
        diagonalMean = mpx_diagonalMean(ts, ceil(mm/2), mm, false);
        hold on;
%         plot([0,length(TD)],[TDhat,TDhat],'--','Color',[0.4, 0.4, 0.4])
%         bar(barX,TDBinned); %%TODO: adjust for weighting
%     %     h = histogram(TDProfile,100); %%TODO: adjust for weighting

        plot(diagonalMean);
        meanMP = mean(mpOriginal,'omitnan');
        dm = mean(diagonalMean,'omitnan');
        plot([0,length(diagonalMean)], [meanMP, meanMP],"--");
        plot([0,length(diagonalMean)], [dm, dm],"--");
        hold off;
        ylim([0,2*sqrt(mm)]);
        set(gca, 'TickDir','out');
        box off;

        mindm = min(diagonalMean);
        score = (mindm-meanMP)/(dm-meanMP);
        title(sprintf("diagonalMeans, min ED score: %.2f, [0,1]", score));
        
        linkaxes([ax1, ax2], 'x');
    end
end


