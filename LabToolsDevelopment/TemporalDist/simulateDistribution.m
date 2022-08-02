function TDDistribution = simulateDistribution(n, iterations)
    TDDistribution = zeros(n,1);

    for ii = 1:n
        for jj = 1:iterations
            tempTD = 0;
            while tempTD == 0
                tempTD = abs(ii-randi([1, n]));
            end
            TDDistribution(tempTD) = TDDistribution(tempTD) + 1;
        end
    end
    TDDistribution = TDDistribution/iterations;

    figure; plot(TDDistribution);
end