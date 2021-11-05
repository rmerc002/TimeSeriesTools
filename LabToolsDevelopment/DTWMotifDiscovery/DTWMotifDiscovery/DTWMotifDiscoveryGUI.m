function [best_so_far,motiffirst,motifsec] = DTWMotifDiscoveryGUI(TS,subseqlen,maxwarp)
   
    gui_name = mpgui_dtw(TS, subseqlen);
    gui_name.drawgui();    
    tic
    [mp, ~] = mpx_v2(TS, subseqlen, subseqlen);%mpx(TS, subseqlen, subseqlen);
    MPlatency = toc;
    sTime = tic;    
    [lb_t,sp_rates,pr_rate,prrate,best_so_far,motiffirst,motifsec,DTW_time]= dtw_motifGUI(TS,subseqlen, maxwarp,mp,gui_name);
    totalTime =toc(sTime)
    
    totalpairs = length(TS)^2;%(length(TS)-2*subseqlen+1)*(length(TS)-2*subseqlen+2)/2;
    fprintf(['Metadata: \n\nTime series length: %d\nSubsequence Length: %d\nMaxwarp: %d\n',...
    'Pair-wise DTW Computation Percentage (After all pruning): %f\nTime to finish ED MP computation: %d seconds\nTime to finish DTW MP computation: %d seconds\n'],length(TS),subseqlen,maxwarp,pr_rate/totalpairs,MPlatency,DTW_time);
    
    
    for i= 1:length(lb_t)
        prunedparis(i) = (length(TS)-prrate(i))^2;%(length(TS)-prrate(i)-2*subseqlen+1)*(length(TS)-prrate(i)-2*subseqlen+2)/2;
        fprintf(['Time to finish lb_keogh (sampling rate:%d): %d seconds --- prunning rate: %f\n'],sp_rates(i),lb_t(i),1-(prunedparis(i)/totalpairs));           
    end    
end