function [warped_T] = add_warping_one_time_series(T,p)
    i = 1:length(T);
    ni = randperm(length(T),floor(length(T) *p));
    i(ni) = [];
    warped_T = smooth(resample(T(i), length(T), length(i)),1);
end