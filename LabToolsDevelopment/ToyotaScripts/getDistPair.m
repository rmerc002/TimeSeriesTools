function [posnan, negnan] = getDistPair(data)
    X1 = data(:,3);
    X2 = data(:,8);
    Y1 = data(:,4);
    Y2 = data(:,9);
    posnan = sqrt((X1-X2).^2 + (Y1-Y2).^2);
    posnan = posnan + randn(length(posnan),1)*1e-5;
    
    X1 = data(:,8);
    X2 = data(:,13);
    Y1 = data(:,9);
    Y2 = data(:,14);
    negnan = sqrt((X1-X2).^2 + (Y1-Y2).^2);
    negnan = negnan + randn(length(negnan),1)*1e-5;
end