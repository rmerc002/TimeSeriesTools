function [best_so_far,first_motif,second_motif] = dtw_motifGUI(ts,subseqlen,maxwarp,initial_best_so_far,first_motif,second_motif,gui_name)
% SWAMP Algorithm - Phase I

% Inputs are time series: ts,
% subsequence length (integer): subseqlen,
% maximum warping offset (integer): maxwarp,
% the seed value for the best motif distance: initial_best_so_far,
% pair of ED motif: first_min, sec_min,
% the gui handle for plotting the data and the motifs: gui_name

% Outputs are motif pair(first_min,sec_min) and their distance(best_so_far)

minspacing = floor(subseqlen/2);
totalpairs = length(ts)^2;
fprintf('\nBeginning Phase |...\n');

%%%%%%%%%%%%%%%%%%% first best-so-far=ED motif distance %%%%%%%%%%%%%%%%%%%%%%%%%%%
best_so_far = initial_best_so_far;
fprintf('best-so-far is %.3f\n',best_so_far);    

%%%%%%%%%%%%%%%%%%% second best-so-far=DTW distance between ED motifs %%%%%%%%%%%%%%%%%%%%%%%%%%%
aa = ts(first_motif:first_motif+subseqlen-1);
bb = ts(second_motif:second_motif+subseqlen-1);
aa = (aa -  mean(aa)) ./ std(aa,1);
bb = (bb - mean(bb)) ./ std(bb,1);
dist = dtw_upd(aa', bb', maxwarp);
if dist < best_so_far 
    best_so_far = dist;    
    gui_name.plotData(first_motif, second_motif)
    gui_name.plotMotifdtw(first_motif, second_motif)
    gui_name.drawgui()
    if gui_name.shouldHalt
        return
    end
end
    
%%%%%%%%%%%%%%%%%%%% Hierarchically downsample(i:1 to 1:1), compute Lower bound MPs and prune TS %%%%%%%%%%%%%%%%%%%%%%
i=64;
dnc = zeros(1,length(ts));
ignores = zeros(1,length(ts));
lb_t = [];
k = 1;

while (length(ts)/i) >= 0.75
    if (round(length(ts)/i) == 1)
        ds = ts;
        dncs = dnc;
    else
        ds = PAA_updated(ts,i);%newpaa(a,i); 
        range = PAA_updated([1:length(ts)],i);
        dncs = PAA_updated(dnc,i); %newpaa(dnc,i);%
    end
   % tic;    
    if round(subseqlen*length(ds)/length(ts)) < 4 || round(subseqlen*length(ds)/length(ts)) < maxwarp
        i = i*2;       
        continue;
    end
%%%%%%%%%%%%%% 1. Compute Downsampled Lower bound %%%%%%%%%%%%%%%%%%
    if (round(length(ts)/i) == 1)
        [mpa_lbk, mpi_lbk, ignores] = LB_Keogh_mp_earlyabandon_updated(ds, max(round(subseqlen*length(ds)/length(ts)),4), max(round(subseqlen*length(ds)/length(ts)),4), maxwarp, dncs, initial_best_so_far);
    else
        [mpa_lbk, mpi_lbk] = LB_Keogh_mp_updated(ds, max(round(subseqlen*length(ds)/length(ts)),4), max(round(subseqlen*length(ds)/length(ts)),4), maxwarp, dncs);    
    end
    
%%%%%%%%%%%%%% 2. Upsample the vector %%%%%%%%%%%%%%%%%%
%     mpa_lbk_stretched = repmat(mpa_lbk', floor(length(ts)/length(ds)), 1);
    scale = sqrt(length(ts)/i); %sp_rates(k));%
%     mpa_lbk_stretched = scale*mpa_lbk_stretched(:);    
%     len_stretch = min(length(ts)-subseqlen+1,length(mpa_lbk_stretched));
%     mpa_lbk_stretched = mpa_lbk_stretched(1:len_stretch);
    
    mpa_lbk_stretched = stepw_upsamp(mpa_lbk, length(ts)-subseqlen+1);
    mpa_lbk_stretched = scale*mpa_lbk_stretched;
        
   %mpi_lbk = (mpi_lbk - 1)*ceil(length(ts)/length(ds)) + 1;
    %mpi_lbk_stretched = repmat(mpi_lbk', round(length(ts)/length(ds)), 1);    
    %mpi_lbk_stretched = mpi_lbk_stretched(:);       
%     mpi_lbk_stretched = mpi_lbk_stretched(1:len_stretch);
    mpi_lbk_stretched = stepw_upsamp(mpi_lbk, length(ts)-subseqlen+1);
    mpi_lbk_stretched = (mpi_lbk_stretched - 1)*floor(length(ts)/length(ds)) + 1;
%%%%%%%%%%%%%% 3. Update Best so far if needed %%%%%%%%%%%%%%%%%%
    if all(mpa_lbk_stretched == Inf)
        return;
    end
    [~,samind] = min(mpa_lbk_stretched);
    temp = mpa_lbk_stretched;
    
    %%%%%%%% Compute DTW distance for k neighbors before and after the min location instead of just the min %%%%%%%%%
    for neigh = (max(1,samind-(floor(length(ts)/length(ds)))):min(length(ts)-subseqlen+1,samind+(floor(length(ts)/length(ds)))))         
        samind_sec = mpi_lbk_stretched(neigh);
        if samind_sec == Inf || abs(samind_sec - samind) <=subseqlen
            continue;
        end
        aa = ts(neigh:neigh+subseqlen-1);
        %samind_sec+subseqlen-1
        bb = ts(samind_sec:samind_sec+subseqlen-1);
        aa = (aa - mean(aa)) ./ std(aa,1);
        bb = (bb - mean(bb)) ./ std(bb,1);
        dist = dtw_upd(aa', bb', maxwarp, best_so_far);    
        if dist < best_so_far        
            best_so_far = dist;
            first_motif = neigh;
            second_motif = samind_sec;  
            gui_name.plotData(first_motif, second_motif)
            gui_name.plotMotifdtw(first_motif, second_motif)
            gui_name.drawgui()
            %pause(1);
            if gui_name.shouldHalt
                return
            end
        end
        mpa_lbk_stretched = temp;
    end
    
               
%%%%%%%%%%%%%% 4. Prune time series if needed %%%%%%%%%%%%%%%%%%  
    mpa_lbk_stretched = temp;    
    dnc(mpa_lbk_stretched > best_so_far) = 1; %%% pruning vector  
    dnc(ignores==1) = 1;   
    prunedpairs = (length(ts)-(sum(dnc)))^2;
    mpa_lbk_stretched(dnc==1)=nan;           
    lb_t(k)=toc;
    fprintf('Completed resolution %d to 1... %.3f%% of pairs have been eliminated. (Elapsed time is %.6f seconds)\n', round(length(ts)/i),100-(prunedpairs/totalpairs)*100,lb_t(k));    
    i = i*2;    
    k = k+1;  
    
    if gui_name.shouldHalt
        return
    end
end
mpa_lbk_stretched(dnc==1)=nan;

%%%%%%%%%%%%%% If all were pruned, return the ED motifs %%%%%%%%%%%%%%%%%%  
if sum(dnc)>=length(ts)-subseqlen          
    gui_name.plotData(first_motif, second_motif)
    gui_name.plotMotifdtw(first_motif, second_motif)
    gui_name.drawgui()
    if gui_name.shouldHalt
            return
    end
    return
end

%%%%%%%%%%%%%%%% Compute DTW Matrix Profile for the pruned TS %%%%%%%%%%%
[mp_sorted,~] = sort(mpa_lbk_stretched);
[best_so_far,first_motif,second_motif,dnc] = dtw_mpGUI(ts', mpi_lbk, mp_sorted, subseqlen, minspacing, maxwarp,dnc,best_so_far,gui_name,first_motif,second_motif);

end

