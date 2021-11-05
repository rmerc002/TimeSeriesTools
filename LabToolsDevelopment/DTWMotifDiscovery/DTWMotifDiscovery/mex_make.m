% build file

% This won't cover all of them, but Clang supports most of the basic gcc
% flags

cc = mex.getCompilerConfigurations('C++');
devel = cc.Manufacturer;
files = 'dtw_upd.cpp';
if strcmp(devel, 'GNU')
   flags = 'CXXFLAGS="-O3 -march=native -funroll-loops -fwrapv -ffp-contract=fast"';
    
elseif strcmp(devel, 'Microsoft')
   % Matlab makes it inconvenient to check enable architectural features
   % Can't easily access setbv or mayIuse if supported. I may revisit this
   % if I add explicit vectorization Auto vectorization and pragma simd
   % probably won't do much for dtw, particularly without some manual
   % unrolling in key areas.
   
   % side side note, MSVC does not have a specific flag for floating point
   % contraction
   
   flags = '/O2 /GL';
    
end


mex('-v', '-R2018a', files, flags);

