function [mp_lbk, mp_ed, mpi_lbk, mpi_ed, t_lb_keogh] = mp_compare(timeseries, subseqlen, minspacing, maxwarp, dnc)
% Provides the matrix profile output based on LB_Keogh and based on 
% Inlining LB_Keogh within this file resulted in absurdly slow performance
% in MATLAB 2019a. I have reverted that inlining step. 
tic
[mp_lbk, mpi_lbk] = LB_Keogh_mp_updated(timeseries, subseqlen, minspacing, maxwarp, dnc);
t_lb_keogh = toc;
%disp(t_lb_keogh)
%disp('lb finished');
[mp_ed, mpi_ed] = mpx(timeseries, minspacing, subseqlen); % parameters have always been in the order minlag followed by subseqlen for all versions I can remember. 
                                                          % I think Eamonn dislikes that, so the new LB_keogh code is the opposite.
                                                         
% figure();
% plot(mp_ed);
% hold on;
% plot(mp_lbk);
% hold off;
% legend('Euclidean distance', 'LB\_Keogh');

end