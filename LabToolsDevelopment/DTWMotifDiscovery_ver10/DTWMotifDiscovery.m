function [ED_motif_distance,first_ed_motif,second_ed_motif,DTW_motif_distance,first_dtw_motif,second_dtw_motif] = DTWMotifDiscovery(TS,subseqlen,maxwarp,suppress)
% Main function for computing Scalable Warping Aware Matrix Profile (SWAMP)

% Inputs are time series: TS,
% subsequence length (integer): subseqlen,
% minimum index offset between subsequence pairs (integer): minlag,
% maximum warping offset (integer): maxwarp,
% flag to enable/disable gui: suppress

% Outputs are ED motif pair and their distance and DTW motif pair and their
% distance

if nargin ==3 
   suppress = false;    
elseif nargin < 3   
   error('incorrect number of input arguments');  
elseif ~isvector(TS)
   error('first argument must be a 1D vector');
elseif ~(isfinite(subseqlen) && floor(subseqlen) == subseqlen) || (subseqlen < 2)  
   error('subsequence length must be an integer value between 2 and the length of the timeSeries');
end      
minlag = floor(subseqlen/2);
[mp, ~] = mpx_v2(TS, minlag, subseqlen);
mp_ed = mp;
[ED_motif_distance,first_ed_motif]=min(mp);
mp(max(first_ed_motif-minlag+1,1):first_ed_motif+minlag-1) = inf;
[~,second_ed_motif] = min(mp);

tic 
fprintf('\nThe input time series was of length %d, with a subsequence length of %d and a warping window of %d.\n',length(TS), subseqlen, maxwarp);

if ~suppress
    gui_name = mpgui_dtw(TS, subseqlen);
    gui_name.drawgui(); 
    gui_name.plotProfile(mp_ed)
    gui_name.plotData(first_ed_motif, second_ed_motif)
    gui_name.plotMotifed(first_ed_motif, second_ed_motif)
    gui_name.plotMotifdtw(first_ed_motif, second_ed_motif)
    gui_name.drawgui()
    if gui_name.shouldHalt
        return
    end
    [DTW_motif_distance,first_dtw_motif,second_dtw_motif] = dtw_motifGUI(TS,subseqlen, maxwarp,ED_motif_distance,first_ed_motif, second_ed_motif,gui_name);
else
    [DTW_motif_distance,first_dtw_motif,second_dtw_motif] = dtw_motif(TS,subseqlen, maxwarp,ED_motif_distance,first_ed_motif, second_ed_motif);
end

totalTime =toc;
fprintf('\nFinished! (Elapsed time is %.6f seconds)\n', totalTime);
fprintf('The best Euclidean distance motif was the pair [%d,%d] with a distance of %.3f\n', first_ed_motif, second_ed_motif,ED_motif_distance);   
fprintf('The best DTW distance motif was the pair [%d,%d] with a distance of %.3f\n', first_dtw_motif, second_dtw_motif,DTW_motif_distance);

end