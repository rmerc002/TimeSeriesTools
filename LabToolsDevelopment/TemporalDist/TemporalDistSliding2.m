function [slidingNearScore] = TemporalDistSliding2(ts, mm, rr, plot_output, mp, mpi)% ts = data(:,5);
    %%%TODO: make a class
    %%%     search for a single TemporalDist in the ts
    %%%     search for a range around the TemporalDist
    %%%     visualize a query for a given Temporal Dist
    %%%     interactive temporalDist highlighting
    %%%     launch app for interactivity
    %%%     selection between weighted vs unweighted frequency

    exclusion_mm = ceil(mm/2);

    if nargin < 4
        plot_output = false;
    end
    if nargin < 3
        rr = length(ts);
    end
    if nargin < 5
        [mp, mpi] = mpx_v3(ts, exclusion_mm, mm);
    end
    ts = reshape(ts,length(ts),1);
 
    subcount = length(ts)-mm+1;
    
    TDMax = 100000;
    
    indices = 1:length(mp);
    indices = indices';
    TDProfile = abs(mpi-indices); %%% Nearest Neighbor Spatial Distances

    weights = linspace(subcount,0.5,subcount); %multiply by weights, or divide by linspace(1/subcount,2,subcounts)
    TDweights = weights(TDProfile);
    numBins = 100;
    binSize = ceil(TDMax/numBins);
    TD2DProfile = zeros(length(ts)-rr+1, numBins);
    for ii = 1:length(ts)-rr+1
        if ii == 1
            tempTDProfile = TDProfile(1:rr);
            tempTDweights = TDweights(1:rr);
            for kk = 1:numBins
                
                binStartIndex = 1+(kk-1)*binSize;
                binEndIndex = binStartIndex + binSize - 1;
                binaryIndices = tempTDProfile > binStartIndex & tempTDProfile < binEndIndex;
                TD2DProfile(ii,kk) = sum(tempTDweights(binaryIndices));
            end
        
        else
            if ii+rr-1 > length(ts)-rr + 1
                break;
            end

            TD2DProfile(ii,:) = TD2DProfile(ii-1,:);

            
            if TDProfile(ii-1) < TDMax
                oldBinIndex = ceil(TDProfile(ii-1)/binSize);
                oldBinValue = TDweights(ii-1);
                TD2DProfile(ii,oldBinIndex) = TD2DProfile(ii,oldBinIndex) - oldBinValue;
            end
            if TDProfile(ii+rr-1) < TDMax
                newBinIndex = ceil(TDProfile(ii+rr-1)/binSize);
                newBinValue = TDweights(ii+rr-1);
                TD2DProfile(ii,newBinIndex) = TD2DProfile(ii,newBinIndex) + newBinValue;
            end
        end
    end

    TD2DProfile = TD2DProfile./rr;

    if plot_output
        figure; 
        imagesc(TD2DProfile');
        ax = gca;
        ax.YDir = 'normal';
        colormap(hot);
    end
    
end


