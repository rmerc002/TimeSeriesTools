function DP = getPlatoDistanceProfile(ts, platos)
        ts = reshape(ts,1,length(ts));
        m = size(platos,2);
        tsNoNan = ts;
        tsNoNan(isnan(ts)) = 0;
        DP = sqrt(4*m)*ones(1,length(ts));
        for platoIndex = 1:size(platos,1)
            plato = platos(platoIndex,:);
            tempDP = MASS_V2(tsNoNan, plato);
            DP = nanmin(DP(1:length(tempDP)), real(tempDP));
        end

        maxVal = sqrt(4*m);
        for nanIndex = 1:length(tsNoNan)-m+1
            if sum(isnan(ts(nanIndex)))
                startIndex = max(1, nanIndex-m+1);
                DP(startIndex:nanIndex) = maxVal;
            end
        end
end