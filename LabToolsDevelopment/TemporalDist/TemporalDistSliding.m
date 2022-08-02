function [slidingNearScore] = TemporalDistSliding(ts, mm, rr, forcePlot, mp, mpi)% ts = data(:,5);
    %%%TODO: make a class
    %%%     search for a single TemporalDist in the ts
    %%%     search for a range around the TemporalDist
    %%%     visualize a query for a given Temporal Dist
    %%%     interactive temporalDist highlighting
    %%%     launch app for interactivity
    %%%     selection between weighted vs unweighted frequency

    exclusion_mm = ceil(mm/2);

    if nargin < 4
        forcePlot = false;
    end
    if nargin < 3
        rr = length(ts);
    end
    if nargin < 5
        [mp, mpi] = mpx_radius(ts, exclusion_mm, mm, rr);
    end
    ts = reshape(ts,length(ts),1);
 
    
    
    indices = 1:length(mp);
    indices = indices';
    TDProfile = abs(mpi-indices); %%% Nearest Neighbor Spatial Distances
    
    %%% normalize by dividing by noise equivalent ED, ignore anti-correlated
    mp = mp/sqrt(2*mm);
    mp = min(1, mp);
    mp = 1-mp;
    
    slidingNearScore = zeros(length(mp),1);
    progressIndices = getExpDistributedSeries(1,length(mp),30);
    TD = zeros(rr,1);
    for jj = 1:length(mp)
        if min(abs(jj - progressIndices)) == 0
            fprintf("starting index %d\n",jj);
        end
        %%% Nearest Neighbor Spatial Profile
        if jj + rr - 1 > length(mp)
                break;
        end

        if jj == 1
            startIndex = max(1, jj-rr);
            endIndex = min(length(mp), jj+rr-1);
            for ii = startIndex:endIndex
                nnsd = TDProfile(ii);
                if nnsd <= 0 || isnan(nnsd)
                    continue;
                end
                TD(nnsd) = TD(nnsd) + mp(ii);
            end
        else
            nnsd = TDProfile(jj+rr-1);
            TD(nnsd) = TD(nnsd) - mpi(jj-1) + mp(jj+rr-1);
        end
    %     TD = TD/length(mp); %%% time series length normalize
        %%% The distribution norm doesn't apply except at the ts
        %%% boundaries
%         distributionNorm = linspace(2,0,length(TD))';
        
%         TD = TD./distributionNorm;
        %%%delete last 20%
%         startIndex = ceil(0.8*length(TD));
%         TD(startIndex:end) = [];
        TDnorm = TD/mean(TD,'omitnan');
    
%         TDhat = mean(TD,'omitnan');
%         RMSE = sqrt(mean((TD - TDhat).^2,'omitnan'));
    
%         weights = linspace(1,0,length(TD));
%         nearScore = (weights*TD)/length(TD);
%         nearScore = 2*(nearScore- 0.5); %%% 0.5 if if there is an even distribution
        [maxValue, maxIndex] = max(TDnorm);
        slidingNearScore(jj) = maxIndex;
    end
    
    %%%%%%%%%%%%%%%%%
    %%%   Plots   %%%
    %%%%%%%%%%%%%%%%%
    if forcePlot == true
        figure;
        tiledlayout(2,1);
        
        ax1 = nexttile();
        plot(ts);
        set(gca, 'TickDir','out');
        box off;
        title("Time Series");
        
        ax2 = nexttile();
        plot(slidingNearScore);
        set(gca, 'TickDir','out');
        box off;
        title("Sliding NearScore");
        
        linkaxes([ax1, ax2], 'x');
    end
end


