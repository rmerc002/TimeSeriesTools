inputPath = "C:\Users\Ryan\Dropbox\Education\UCR\Research\0_Projects\ActiveResearch\Toyota\Data\Car-following V2";

for ii = 1:12

    anomFilePath = fullfile(inputPath, "Data_anom",sprintf("anom_veh%d.txt",ii));
    baseFilePath = fullfile(inputPath, "Data_base",sprintf("base_veh%d.txt",ii));

    dataAnom = importdata(anomFilePath);
    dataBase = importdata(baseFilePath);
    

    %%% Distance between following and anomalous
    [posAnom, negAnom] = getDistPair(dataAnom);
    [posBase, negBase] = getDistPair(dataBase);
    
    sampleRate = 10;
    mm = 5*sampleRate;
    
    MP_AB_Anom = mpx_ABBA_v2(posAnom, negAnom, mm);
    MP_AB_Base = mpx_ABBA_v2(posBase, negBase, mm);

    figure;
    tiledlayout(3,1);
    ax1 = nexttile();
    hold on;
    plot(zscore(posBase));
    plot(zscore(negBase));
    hold off;
    title('Distance Between Vehicles: Base');

    ax1 = nexttile();
    hold on;
    plot(zscore(posAnom));
    plot(zscore(negAnom));
    hold off;
    title('Distance Between Vehicles: Anom');

    ax3 = nexttile();
    hold on;
    plot(MP_AB_Base);
    plot(MP_AB_Anom);
    hold off;
    ylim([0,sqrt(2*mm)]);
    title('Dissimilarity to Base Behavior');

    linkaxes([ax1, ax2, ax3],'x');
   
end

