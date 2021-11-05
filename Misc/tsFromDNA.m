function [timeSeries] = tsFromDNA(seq)
% This is a front end to dnatots written as a one off for the .mat files
% found at
% https://www.biorxiv.org/content/biorxiv/early/2018/08/20/394932.1.full.pdf
% https://github.com/grandhawa/MLDSP/tree/master/DataBase

if ~isvector(seq) || ~iscell(seq)
    error('expected a 1D cell array');
end

rows = 0;

for i = 1 : length(seq)
    rows = max(rows, length(seq{i}));
end

timeSeries = NaN(rows, length(seq));

for i = 1 : length(seq)
    timeSeries(1 : length(seq{i}), i) = dnatots(seq{i});
end

plotSpecial(timeSeries, rows);

end


function plotSpecial(timeSeries, longest)
fg = figure();
ax = axes(fg);
plot(ax, timeSeries);
ax.Box = 'off';
ax.YTick = [];
ax.XTick = [0 longest];
ax.XLim = [0 longest];
end

function ts = dnatots(dna)

  % weights = [2,-1,1,-2]; %prev
weights = [1,-2,2,-1]; %all sum(abs)

ts = zeros(length(dna), 1);
if dna(1) == 'A'
    ts(1) = weights(1);
elseif dna(1) == 'C'
    ts(1) = weights(2);
elseif dna(1) == 'G'
    ts(1) = weights(3);
elseif dna(1) == 'T'
    ts(1) = weights(4);
else
    error('unrecognized character ''%c'' at index %g', ts(1), 1);
end

for i = 2 : length(dna)
    if dna(i) == 'A'
        ts(i) = ts(i - 1) + weights(1);
    elseif dna(i) == 'C'
        ts(i) = ts(i - 1) + weights(2);
    elseif dna(i) == 'G'
        ts(i) = ts(i - 1) + weights(3);
    elseif dna(i) == 'T'
        ts(i) = ts(i - 1) + weights(4);
    else
        error('unrecognized character ''%c'' at index %g', ts(i), i);
    end
end

end
