%  DTW_upd.m Help file for DTW_interf.m MEX-file.
%  
%  This returns the dynamic time warping distance for x,y using the L2
%  norm.
% 
%  Inputs are 
%  x:  first time series
%  y:  second time series
%  w:  max warping window
%  threshold: maximum allowable distance
%
%  This will return either dtw distance or the first partial calculation
%  which exceeds the threshold. It will not return a warping path.
%  
%