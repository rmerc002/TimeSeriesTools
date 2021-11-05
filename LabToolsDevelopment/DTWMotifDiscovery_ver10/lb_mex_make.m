

% do not run this from another directory or you won't be able to find the
% necessary files.

filelist = 'lb_upd.c';
cc = mex.getCompilerConfigurations('C');
if length(cc) > 1
    % use default compiler
    cc = cc(1);
elseif isempty(cc)
    error('This package requires a C++ compiler to be registered with Matlab. Currently there is no default compiler set.');
end
if strcmp(cc.Manufacturer, 'GNU') || strcmp(cc.Manufacturer, 'Apple')  
    % Flags appropriate for Clang or GCC on x86_64 platforms (most matlab
    % users). Feel free to make adjustments. These are based on timings and
    % inspecting the output assembly. Unroll loops tends to result in
    % better scheduling of certain arithmetic heavy functions.
    
    % These will enable AVX, AVX2, and FMA3 if available.
    flags = 'CXXFLAGS="-O3 -march=native -funroll-loops -fwrapv -fopenmp -ffp-contract=fast"';

    %  Todo: retest -fprefetch-loops-arrays, this should not really matter a lot less
    %  with a good choice of tile size

elseif strcmp(cc.Manufacturer, 'Microsoft')
    % Haven't tested this one.. also need to 
    
    % If you're using Windows with MSVC and these settings are problematic,
    % you're probably using an old cpu model or an older OS. 
    % remove the /Arch:AVX2 and it should work
    flags = 'CXXFLAGS="-I../include/ CXXFLAGS=/O2 /Arch:AVX2 /fp:fast /GL"';
else
    % If you are hitting this section and want to hack it yourself, code performance relies on 
    % - partial loop unrolling (enabled with unknown trip count)
    % - simd auto-vectorization options, and disabled strided load patterns on targets supporting AVX. 
    % - contraction a*b + c into fma expressions
    
    %If this makes no sense to you, send an email with the full
    % output of mex.getCompilerConfigurations and any available info on your target platform. 
   % warning('Unfortunately we have not tested this exact configuration. Mathworks does not provide a reference for compiler vendor name as they use it here.');
end


mex('-v', '-R2018a', filelist, flags);
