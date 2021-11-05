function [best_so_far,first_dtw_motif,second_dtw_motif,dnc] = dtw_mpGUI(ts,mp_ind,mp_sorted,subseqlen,minlag,warpmax,dnc,best_so_far,gui_name,first_dtw_motif,second_dtw_motif)
% SWAMP Algorithm - Phase II

% Inputs are time series: ts,
% Matrix Profile Index using LB Keogh: mp_ind,
% sorted residual Matrix Profile: mp_sorted,
% subsequence length (integer): subseqlen,
% minimum index offset between subsequence pairs (integer): minlag,
% maximum warping offset (integer): warpmax,
% vector of true and false booleans showing which subsequences to prune: dnc
% the lowest DTW motif distance found so far: best_so_far,
% the gui handle for plotting the data and the motifs: gui_name,
% the best motif pair found so far: first_dtw_motif,second_dtw_motif

% Outputs are motif pair(motiffirst,second_dtw_motif) and their distance(best_so_far)
whentoabandon2 = [];

mu = movmean(ts, [0 subseqlen-1], 'Endpoints', 'discard');
sig = movstd(ts, [0 subseqlen-1], 1, 'Endpoints', 'discard');
subcount = length(ts) - subseqlen + 1;

tr= (dnc==1);
mp_ind(tr) = -1;

fprintf('\nBeginning Phase ||...\n');
fprintf('best-so-far is %.3f\n', best_so_far);
prev = sum(dnc);

for i = 1 : subcount      
    id = (mp_ind(i));    
    if id==-1 || dnc(id)==1
        continue;
    end  
    neighbors = [id + minlag : subcount];  
    mp_ind(dnc==1) = -1;
    neigh = mp_ind(id);
    neighbors(neighbors==neigh)= [];
    neighbors = [neigh,neighbors];
    
    for j = neighbors       
        idp = j;        
        if idp==-1 || dnc(idp)==1
            continue;
        end       
        
        a = (ts(id : id + subseqlen - 1) - mu(id))./sig(id);
        b = (ts(idp : idp + subseqlen - 1) - mu(idp))./sig(idp);
        
        %%%%%%%%%%%%%% Compute LB_Kim lower bound distance and skip DTW computation if higher than best-so-far %%%%%%%%%%%%%%%%%%          
        lb_Kim = max(abs(a(1)-b(1)),abs(a(end)-b(end)));            
        if lb_Kim >= best_so_far                 
             continue;
             
        %%%%%%%%%%%%%% Compute LB_Keogh lower bound distance and skip DTW computation if higher than best-so-far %%%%%%%%%%%%%%%%%%  
        else
            Ua = movmax(a,[warpmax warpmax]);
            La = movmin(a,[warpmax warpmax]);                   
            [LB_Keogh, wb] = lb_upd(b,Ua,La,best_so_far);   
            whentoabandon2 = [whentoabandon2; wb];
            if sqrt(LB_Keogh) >=best_so_far               
                continue;
            end
        end        
        
        %%%%%%%%%%%%%% Compute the DTW distance if the previous two checks failed %%%%%%%%%%%%%%%%%%         
        [dist,ea] = dtw_upd(a,b,warpmax,best_so_far);    
        whentoabandon2 = [whentoabandon2; ea];
        %dist = dtw_itakura(a, b, warpmax);
        
        %%%%%%%%%%%%%% Update best-so-far if needed %%%%%%%%%%%%%%%%%%
        if dist < best_so_far
            best_so_far = dist;
            first_dtw_motif = id;
            second_dtw_motif = idp;     
            gui_name.plotData(id, idp)
            gui_name.plotMotifdtw(id, idp)
            gui_name.drawgui()  
            %pause(1);
            if gui_name.shouldHalt
                return 
            end
            dnc(mp_ind(mp_sorted>best_so_far)>0) = 1;            
            mp_ind(mp_sorted>best_so_far) = -1;
            mp_sorted(mp_sorted>best_so_far) = -1;  
            prrate = sum(dnc);
            prunedpairs = ((length(ts)-prrate)^2)/(length(ts)^2);
            
            update_time = toc;
            fprintf('best-so-far updated to %.3f, %.3f%% of pairs have been eliminated. (Elapsed time is %.6f seconds) \n', best_so_far, 100-(prunedpairs)*100, update_time);
            prev = prrate + prev;
        end
        if gui_name.shouldHalt
            return 
        end
    end
end

end

