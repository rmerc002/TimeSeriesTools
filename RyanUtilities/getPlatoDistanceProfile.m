function DP = getPlatoDistanceProfile(ts, platos)
        ts = reshape(ts,length(ts),1);
        if size(platos,1) == 1 || size(platos,2) == 1
            platos = reshape(platos,1,length(platos));
        end

        mm = size(platos,2);
        tsNoNan = ts;
        tsNoNan(isnan(ts)) = 0;
        DP = sqrt(4*mm)*ones(length(ts),1);
        for platoIndex = 1:size(platos,1)
            plato = platos(platoIndex,:)';
            tempDP = real(MASS_V2(tsNoNan, plato));
            DP = min(DP(1:length(tempDP)), tempDP,"omitnan");
        end

        maxVal = sqrt(4*mm);
        for nanIndex = 1:length(tsNoNan)-mm+1
            if sum(isnan(ts(nanIndex)))
                startIndex = max(1, nanIndex-mm+1);
                DP(startIndex:nanIndex) = maxVal;
            end
        end
end