function [posnan, negnan] = singleVehiclePosDiff(data)
    X1 = data(1:end-1,3);
    X2 = data(2:end,3);
    Y1 = data(1:end-1,4);
    Y2 = data(2:end,4);
    posnan = sqrt((X1-X2).^2 + (Y1-Y2).^2);
    posnan = posnan + randn(length(posnan),1)*1e-5;
    
    X1 = data(1:end-1,8);
    X2 = data(2:end,8);
    Y1 = data(1:end-1,9);
    Y2 = data(2:end,9);
    negnan = sqrt((X1-X2).^2 + (Y1-Y2).^2);
    negnan = negnan + randn(length(negnan),1)*1e-5;
end