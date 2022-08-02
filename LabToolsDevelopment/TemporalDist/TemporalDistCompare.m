function TemporalDistCompare(TD1, TD2)
    figure;
    hold on;
    plot(TD2);
    plot(-TD1);
%     scatter(1:length(TD1),TD1-TD2,'_', 'MarkerEdgeColor',[0.5,0.1,0.5]);
    hold off;
    set(gca, 'TickDir','out');
    box off;


    
    TDPercentChangeNeg = 100*max(0, TD1-TD2)./(TD1+ones(length(TD1),1)*1e-2) + ones(length(TD1),1)*1e-5;
    TDPercentChangePos = 100*max(0,TD2-TD1)./(TD1+ones(length(TD1),1)*1e-2) + ones(length(TD1),1)*1e-5;

%     TDPercentChangeNeg = max(0, TD1-TD2) + ones(length(TD1),1)*1e-5;
%     TDPercentChangePos = max(0,TD2-TD1) + ones(length(TD1),1)*1e-5;
    
    %%% Plots increased amplitude in TD1 as a loss (red)
    %%%       increased amplitude in TD2 as a gain (blue)
    
    figure;
    hold on;
    plot(TDPercentChangePos);
    plot(TDPercentChangeNeg);
    hold off;
%     set(gca, 'YScale', 'log');
    set(gca, 'TickDir','out');
    box off;
end

