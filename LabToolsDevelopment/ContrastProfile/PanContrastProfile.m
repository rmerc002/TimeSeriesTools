function [PanCP, PanPlato] = PanContrastProfile(positiveTS, negativeTS, subLengthMin, subLengthMax, forcePlot)
    %%% Input:
    %%%   positiveTS: A time series containing at least two instances of a desired behavior
    %%%   negativeTS: A time series containing zero instances of a desired behavior
    %%%   lengthMin: Smallest subsequence length
    %%%   lengthMax: Largest subsequence length
    %%%   forcePlot: (optional) Display the following two plots
    %%%
    %%% Output:
    %%%   PanCP: Contrast Profile, which indicates subsequence within
    %%%       positiveTS that are less conserved in negativeTS
    %%%   PanPlato: The subsequence that most distinguishes 
    %%%          positiveTS from negativeTS. 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                            Look HERE
    %%% When true, iterates over a fraction of the lengths, 
    %%%  exponentially spaced. Recommended for initial exploration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    onetwoskipafewninetynineonehundred = true;%false;
    
    if nargin == 4
        forcePlot = false;
    elseif nargin ~= 5
        error('incorrect number of input arguments');
    elseif ~isvector(positiveTS) || ~isvector(negativeTS)
        error('first and second arguments must be a 1D vector');
    elseif ~(isfinite(subLengthMin) && floor(subLengthMin) == subLengthMin) || ...
           ~(isfinite(subLengthMax) && floor(subLengthMax) == subLengthMax) || ...
           (subLengthMin < 2) ||...
           (subLengthMax <= subLengthMin)
        error('subsequence length must be an integer value between 2 and the length of the timeseries');
    end

    if onetwoskipafewninetynineonehundred == true
        subLengths = getExpDistributedSeries(subLengthMin, subLengthMax,30);
    else
       subLengths = subLengthMin:subLengthMax;
    end
    
    PanCP = NaN(length(subLengths), length(positiveTS));
    PanPlato = NaN(length(subLengths), subLengths(end));
    tic;
    iterationTimes = NaN(length(subLengths),1);
    for i = 1:length(subLengths)
        tic;
        subLength = subLengths(i);
        [plato, ~, CP] = ContrastProfile(positiveTS, negativeTS, subLength); 
        PanPlato(i,1:length(plato)) = plato;
        PanCP(i,1:length(CP)) = CP;
        iterationTimes(i) = toc;
        if forcePlot == true fprintf("\tCompleted CP_%d. Max contrast value = %.2f, %.2f seconds\n", subLength, max(CP), iterationTimes(i));, end
    end
    if forcePlot == true fprintf("Total Time Elapsed: %.0f seconds\n\n",sum(iterationTimes));, end

    %%%%%%%%%%%%%%%%%%%%%%%
    %%% Plot 3D Surface %%%
    %%%%%%%%%%%%%%%%%%%%%%%
    if forcePlot == true
        %%% Requires PanCP, subLengths

        [maxPanCP,maxPanCPIndex] = max(max(PanCP,[],1));
        [~,maxPanCPSubLengthIndex] = max(max(PanCP,[],2));
        maxPanCPSubLength = subLengths(maxPanCPSubLengthIndex);

        %redundant as far as calculation, but this function has the plots already
        % [plato, plato_twin, CP] = ContrastProfile(pos, neg, maxPanCPSubLength); 

        %surface plot
        figure;
        [X,Y] = meshgrid(1:size(PanCP,2),subLengths);

        hold on;
        surf(X,Y,PanCP)
        view(-45,45);
        shading interp
        grid on;
        scatter3(maxPanCPIndex,maxPanCPSubLength,maxPanCP,'filled','r')
        hold off;

        xlabel("Time");
        zlabel("Contrast Profile")
        ylabel("Subsequence Length")
        title("Pan Contrast Profile Surface, \color{red}Plato \color{black}Index");

        fprintf("%d: subLength with the maximum CP value\n", maxPanCPSubLength);
    end
end
function subLenSeries = getExpDistributedSeries(startLen, endLen,numSteps)
    %%% Purpose: Reduce space. matrix profile distances are pretty similar from
    %%% from one subLen to the next. By trial and error, adding the square root
    %%% of the current subLen seems to produce a good distribution.
    powerMin = log10(startLen);
    powerMax = log10(endLen);
    % numSteps = 80;
    powerStep = (powerMax-powerMin)/numSteps;
    powers = powerMin:powerStep:powerMax;
    subLenSeries = unique(ceil(power(10,powers)));
end