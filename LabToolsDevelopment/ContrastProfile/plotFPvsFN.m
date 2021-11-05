function plotFPvsFN(fpfnResults, figureSavePath)
    xmax = max(fpfnResults(:,1));
    ymax = max(fpfnResults(:,2));
    
    grayColor = [0.6,0.6,0.6];
    fig = figure;
    scatter(fpfnResults(:,1),fpfnResults(:,2),10,'MarkerFaceColor','r','MarkerEdgeAlpha',0,'MarkerFaceAlpha',.5);
    hold on; 
    plot([0,xmax],[0,0],'--','Color', grayColor);
    plot([0,0],[0,ymax],'--','Color', grayColor)
    
    xlim([-1,xmax]);
    ylim([-1,ymax]);
    set(gca,'xtick',[],'ytick',[]);
    
    figureFileName = "FPvsFN";
    saveas(fig, fullfile(figureSavePath, figureFileName + ".fig"));
    saveas(fig, fullfile(figureSavePath, figureFileName +".png"));
    close(fig);
end