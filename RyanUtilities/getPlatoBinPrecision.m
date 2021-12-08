function binPrecision = getPlatoBinPrecision(platos, ts, labels, binSize)
    DP = getPlatoDistanceProfile(ts, platos);
    m = size(platos,2);
    binPrecision = getBinPrecision(DP, labels, m, binSize);
end



