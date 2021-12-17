function [platos, platoIndices, platoContrast] = RelativeFrequencyContrastProfile_TopK(tsPos, tsNeg, m, maxFreq, K, forcePlot)
    %%% Input:
    %%%   positiveTS: A time series containing at least two instances of a desired behavior
    %%%   negativeTS: A time series containing zero instances of a desired behavior
    %%%   m: The approximate window size of the desired behavior
    %%%   forcePlot: (optional) Display the following two plots
    %%%
    %%% Output:
    %%%   plato: The subsequence that most distinguishes 
    %%%          positiveTS from negativeTS
    %%%   plato_twin: The nearest neighbor of plato within positiveTS
    %%%   CP: Contrast Profile, which indicates subsequence within
    %%%       positiveTS that are less conserved in negativeTS

    initial__pos_subcount = length(tsPos) - m + 1;
    initial__neg_subcount = length(tsNeg) - m + 1;
    
    if nargin == 5
        forcePlot = false;
    elseif nargin ~= 6
        error('incorrect number of input arguments');
    elseif ~isvector(tsPos) || ~isvector(tsNeg)
        error('first and second arguments must be a 1D vector');
    elseif ~(isfinite(m) && floor(m) == m) || (m < 2) || (initial__pos_subcount < 2) || (initial__neg_subcount < 2)
        error('subsequence length must be an integer value between 2 and the length of the timeseries');
    end
    
    tsPos = reshape(tsPos, 1, length(tsPos));
    tsNeg = reshape(tsNeg, 1, length(tsNeg));
%     tsPos(isnan(tsPos)) = nanmean(tsPos);
% 
%     tsNeg(isnan(tsNeg)) = nanmean(tsNeg);
    
    %Change to row vector for internal use
%     if size(tsPos,1) == 1 
%         %save the data orientation for matching output format to input
%         %give orientation priority to positiveTS
%         tsPos = tsPos';
%     end
%     if size(tsNeg,1) == 1 
%         %give orientation priority to positiveTS, do not save negative
%         %orientation if different than positive
%         tsNeg = tsNeg';
%     end
    
    platos = nan(K, m);
    platoIndices = nan(1,K);
    platoContrast = nan(1,K);

    %%% Matrix profile self-join using positiveTS
    [RFMP_AA, RFMP_AA_Indices] = RelativeFrequencyMatrixProfile(tsPos, tsPos, m, maxFreq);
    RFMP_AA = real(RFMP_AA);
    %TODO: need to figure out the bound of MASS_V2 in mpxAAKNN
    
    %%% Euclidean distance values above sqrt(2*m) are equivalent to
    %%%   anti-correlated values
    RFMP_AA = clipMatrixProfileAmplitude(RFMP_AA, m);
    %%% pad with NaN to make future comparisons between matrix profiles
%     padLength = length(positiveTS) - length(MP_AA) + 1;
%     MP_AA = [MP_AA;NaN(padLength,1)];

    RFMP_AB_history = sqrt(2*m)*ones(maxFreq, length(tsPos));
    pastPlato = tsNeg;
    
    [RFMP_AB_history, RFMP_AB_Indices] = RelativeFrequencyMatrixProfile(tsPos, pastPlato, m, maxFreq);
    RFMP_AB_history = real(RFMP_AB_history);
    RFMP_AB_history = clipMatrixProfileAmplitude(RFMP_AB_history, m);
    for ki = 1:K
        %%% Matrix profile AB-join between positiveTS and negativeTS
        if ki >= 2
            MP_AB = mpx_ABBA_v2(tsPos, pastPlato, m);
            MP_AB = real(MP_AB);
            MP_AB = reshape(MP_AB,1,length(MP_AB));
            MP_AB = [MP_AB, nan(1,length(tsPos)-length(MP_AB))];
            for mf = 1:maxFreq
                RFMP_AB_history(mf,:) = nanmin(RFMP_AB_history(mf,:), MP_AB);
            end
        end
        
        
        
        
        
%         %%%If a non-overlapping candidates is less than maxFreq
%         effectiveMaxFreq = min(size(RFMP_AA,1), size(RFMP_AB,1));
%         if maxFreq > effectiveMaxFreq
%             fprintf("Warning!: Choice of maxFreq=%d was too large for dataset, setting to %d",maxFreq, effectiveMaxFreq);
%             RFMP_AA = RFMP_AA(1:effectiveMaxFreq,:);
%             RFMP_AB = RFMP_AB(1:effectiveMaxFreq,:);
%         end



        %%% Contrast Profile
        RFCP = RFMP_AB_history - RFMP_AA;
        %%% Normalize values to the range [0,1]
        RFCP = normalizeContrastProfileAmplitude(RFCP, m);


        %%% plato is the subsequence in positiveTS corresponding to index with
        %%%   largest contrast profile value
%         [maxContrastValues, platoIndices] = max(RFCP,[],2);
%         [maxContrastValue, KBest] = max(maxContrastValues);
%         platoIndex = platoIndices(KBest);
%         plato = tsPos(platoIndex:platoIndex+m-1);

        %%%Alternative Plato calculation
        CPmean = zeros(1,size(RFCP,2));
        for ti = 1:size(RFCP,2)-m+1
            subsequence = tsPos(ti:ti+m-1);
            if sum(isnan(subsequence)) > 0
                continue
            end
           CPmean(ti) = norm(RFCP(:,ti))/sqrt(maxFreq); %% root mean squared
        end
        [maxContrast, platoIndex] = max(CPmean);
        startIndex = platoIndex;
        endIndex = startIndex + m - 1;
        platos(ki,:) = tsPos(startIndex:endIndex);
        platoIndices(ki) = platoIndex;
        platoContrast(ki) = maxContrast;


        exclusionLength = m;
        startIndex = max(1, platoIndex-exclusionLength);
        endIndex = min(length(tsPos), ceil(platoIndex + m - 1 + exclusionLength));
        pastPlato = tsPos(startIndex:endIndex);
        pastPlato = reshape(pastPlato,1,length(pastPlato));
%         RFMP_AA(:,startIndex:endIndex) = nan;
    end
    platoContrast = platoContrast(1,:);
end

function [mp] = clipMatrixProfileAmplitude(mp,m)
    mp = min(sqrt(2*m),mp,'includenan'); %clip the anti-correlated values
    mp = real(mp);
end

function [cp] = normalizeContrastProfileAmplitude(cp, m)
    cp = cp/(sqrt(2*m)); %normalize between 0,1
    %discard negative values. These occur when a subsequence in T+ has a
    %closer nearest neighbor in T-. It's not a behavior of interest for
    %this research.
    cp = nanmax(0, cp); 
end