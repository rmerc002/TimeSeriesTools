function [platos_Anytime,platoIndices_Anytime, CPs_Anytime, percentages_Anytime] = ContrastProfileAnytime(positiveTS, negativeTS, subLength, forceplot)
    %%%Output
    %%% platos_Anytime: Matrix of Platos for each percentage step completed.
    %%%   shape: |percentages_Anytime| x subLength 
    %%% platoIndices_Anytime: Vector of indices of plato position within
    %%%   positiveTS
    %%%   shape: |percentages_Anytime| x 1
    %%% CPs_Anytime: Matrix of Contrast Profiles for each percentage step
    %%%   completed.
    %%%   shape: |percentages_Anytime| x |positiveTS|-subLength+1
    %%% percentages_Anytime: Vector of percentage steps completed
    %%%   
    tic;
    initial__pos_subcount = length(positiveTS) - subLength + 1;
    initial__neg_subcount = length(negativeTS) - subLength + 1;
    
    if nargin == 3
        forceplot = false;
    elseif nargin ~= 4
        error('incorrect number of input arguments');
    elseif ~isvector(positiveTS) || ~isvector(negativeTS)
        error('first argument must be a 1D vector');
    elseif ~(isfinite(subLength) && floor(subLength) == subLength) || (subLength < 2) || (initial__pos_subcount < 2) || (initial__neg_subcount < 2)
        error('subsequence length must be an integer value between 2 and the length of the timeseries');
    end
    
    
    dataOrientation = 0; %0 for row, 1 for column
    %Change to row vector for internal use
    if size(positiveTS,1) == 1 
        %save the data orientation for matching output format to input
        %give orientation priority to positiveTS
        dataOrientation = 1; 
        positiveTS = positiveTS';
    end
    if size(negativeTS,1) == 1 
        %give orientation priority to positiveTS, do not save negative
        %orientation if different than positive
        negativeTS = negativeTS';
    end
   
    %%%%%%%%%%%%%%%%%%%%%%
    %%% Anytime Mode 2 %%%
    %%%%%%%%%%%%%%%%%%%%%%
    fprintf("Starting Anytime Self-Join Matrix Profile\n");
    anytimeMode = 2;
    %%% Matrix profile self-join using positiveTS
    [MPs_AA, percentages_Anytime] = interactiveMatrixProfileVer3_AA(positiveTS, subLength, anytimeMode);
    %%% Euclidean distance values above sqrt(2*subLength) are equivalent to
    %%%   anti-correlated values
    MPs_AA = clipMatrixProfileAmplitude(MPs_AA, subLength);
    %%% pad with NaN to make future comparisons between matrix profiles
    padLength = length(positiveTS) - length(MPs_AA) + 1;
    MPs_AA = [MPs_AA,NaN(length(percentages_Anytime),padLength)];

    fprintf("Starting Anytime AB-Join Matrix Profile\n");
    %%% Matrix profile AB-join between positiveTS and negativeTS
    [MPs_AB, percentages_Anytime] = interactiveMatrixProfileVer3_AB(positiveTS, negativeTS, subLength, anytimeMode);
    %%% Euclidean distance values above sqrt(2*subLength) are equivalent to
    %%%   anti-correlated values
    MPs_AB = clipMatrixProfileAmplitude(MPs_AB, subLength);
    %%% pad with NaN to make future comparisons between matrix profiles
    padLength = length(positiveTS) - length(MPs_AB) + 1;
    MPs_AB = [MPs_AB,NaN(length(percentages_Anytime),padLength)];
    
    fprintf("Starting Contrast Profile Calculation\n");
    %%% Contrast Profile
    dim1 = min(size(MPs_AB,1), size(MPs_AA,1));
    dim2 = min(size(MPs_AB,2), size(MPs_AA,2));
    CPs_Anytime = MPs_AB(1:dim1,1:dim2) - MPs_AA(1:dim1,1:dim2);
    %%% Normalize values to the range [0,1]
    CPs_Anytime = normalizeContrastProfileAmplitude(CPs_Anytime, subLength);
    CPs_Anytime(1,:) = zeros(1,size(CPs_Anytime,2));%null hypothesis is that pos and neg are the same, so CP is zero
    
    
    %%% plato is the subsequence in positiveTS corresponding to index with
    %%%   largest contrast profile value
    [~, platoIndices_Anytime] = max(CPs_Anytime,[],2);
    platos_Anytime = [NaN(1,subLength)];
    for i=1:length(platoIndices_Anytime)
        platos_Anytime(i,:) = positiveTS(platoIndices_Anytime(i):platoIndices_Anytime(i) + subLength - 1);
    end
    
    CP = CPs_Anytime(end,:);
    plato = platos_Anytime(end,:);
    MP_AA = MPs_AA(end,:);
    MP_AB = MPs_AB(end,:);
    %%%%%%%%%%%%%
    %%% PLOTS %%%
    %%%%%%%%%%%%%
    if forceplot == true
        redColor = [0.73,0.05,0];
        greenColor = [0,0.73,0.41]; 
        blueColor = [0,0.29,0.73];
        grayColor = [0.75,0.75,0.75];
        lightBlueColor = [0.01, 0.83,0.99];
        platoColor = [129/255, 51/255, 144/255];
        platoTwinColor = [115/255, 170/255, 43/255];

        tsLength = length(positiveTS);
        maxTSLength = max(length(positiveTS),length(negativeTS));
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        %%% Plot Convergence %%%
        %%%%%%%%%%%%%%%%%%%%%%%%
        errorsAnytime = [];
        for i = 1:size(CPs_Anytime,1)
            error = CPs_Anytime(i,:) - CPs_Anytime(end,:);
            errorsAnytime(end+1) = sqrt(mean(error.*error));%RMSE
        end
        
        figure;
        hold on;
        plot(percentages_Anytime * 100, errorsAnytime);
        hold off;
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        %%% Plot Max vs 100% %%%
        %%%%%%%%%%%%%%%%%%%%%%%%
        figure; 
        hold on; 
        area(max(CPs_Anytime,[],1));colororder([1,0,0]); 
        plot(CPs_Anytime(end,:),'g','LineWidth',2); 
        hold off;
        set(gca,'xtick',[1,length(positiveTS)],'ytick',[0,1], 'TickDir','out');
        box off;
        
        xlim([1,length(positiveTS)]);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Plot Plato Convergence %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(gcf, 'Position', [0,100,600,400]);
        
        subplot(4,6,[3:6]);
        plot(positiveTS,'Color',blueColor);
        formattedTitle = sprintf("\\color[rgb]{%f,%f,%f}T^{(+)}\\color{black}: Positive Time Series",blueColor(1),blueColor(2),blueColor(3));
        title(formattedTitle);
        set(gca,'xtick',[1,length(positiveTS)],'ytick',[], 'TickDir','out');
        box off;
        
        subplot(4,6,[9:12]);
        hold on;
        plot(CPs_Anytime(end,:),'Color',grayColor);
        scatter(platoIndices_Anytime(end), 1,10,'MarkerFaceColor',platoColor,'MarkerEdgeColor',platoColor,'MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0); 
        hold off;
        xlim([1,maxTSLength]);
        ylim([0,1]);
        formattedTitle = sprintf("\\color[rgb]{%f,%f,%f}Contrast Profile, \\color[rgb]{%f,%f,%f}Plato",grayColor(1), grayColor(2), grayColor(3), platoColor(1), platoColor(2), platoColor(3));
        title(formattedTitle);
        set(gca,'xtick',[1,length(positiveTS)],'ytick',[0,1], 'TickDir','out');
        box off;
        
        subplot(4,6,[13,14,19,20]);
        hold on;
        plot(errorsAnytime,percentages_Anytime * 100);
        set(gca,'ytick',[0,50,100], 'TickDir','out','XDir','reverse');
        ylabel("Percent Complete");
        title("RMSE");
        xlabel("Error (reversed)");
        
        subplot(4,6,[15:18,21:24]);
        scatter(platoIndices_Anytime,percentages_Anytime * 100,2,'filled');
        xlim([1,length(positiveTS)]);
        title("Convergence of Plato Index");
        
        set(gca,'xtick',[1,length(positiveTS)],'ytick',[], 'TickDir','out');
        box off;
    end
    fprintf("Total time: %.2f seconds\n",toc);
end

function [mp] = clipMatrixProfileAmplitude(mp,subLength)
    mp = min(sqrt(2*subLength),mp,'includenan'); %clip the anti-correlated values
    mp = real(mp);
end

function [cp] = normalizeContrastProfileAmplitude(cp, subLength)
    cp = cp/(sqrt(2*subLength)); %normalize between 0,1
    %discard negative values. These occur when a subsequence in T+ has a
    %closer nearest neighbor in T-. It's not a behavior of interest for
    %this research.
    cp = nanmax(0, cp); 
end
