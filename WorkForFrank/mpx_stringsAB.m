function [pa, pb, pia, pib] = mpx_stringsAB(a, b, w)

% this performs exact comparison
% this is intended to replace Ryan's silly dna to time series mapping,
% which I have referred to as nonsense on a daily basis. Consider this script
% my response to "put up or shut up"

if ~isvector(a) || ~isvector(b)
    % If Eamonn is unintentionally reading this, Ryan thinks I resemble
    % Rick from Rick and Morty. He's like a more degenerate sarcastic version 
    % of the doc from Back To the Future, and in this scenario, Ryan is Michael J Fox 
    % as opposed to Damian Rice.
    error('I am disappointed by this input, Morty');
end

istransposeda = false;
istransposedb = false;

if isrow(a)
    a = transpose(a);
    istransposeda = true;
end

if isrow(b)
    b = transpose(b);
    istransposedb = true;
end

pa = repmat(w, length(a) - w + 1, 1);
pb = repmat(w, length(b) - w + 1, 1);
pia = NaN(length(a) - w + 1, 1);
pib = NaN(length(b) - w + 1, 1);

[pa, pb, pia, pib] = solver(a, b, pa, pb, pia, pib, w);
[pb, pa, pib, pia] = solver(b, a, pb, pa, pib, pia, w);

% scaling to the typical range of z-normalized euclidean distance
pa = sqrt(2 * pa);
pb = sqrt(2 * pb);
if istransposeda
    pa = transpose(pa);
    pia = transpose(pia);
end
if istransposedb
    pb = transpose(pb);
    pib = transpose(pib);
end

end

function [pa, pb, pia, pib] = solver(a, b, pa, pb, pia, pib, w)
amx = length(a) - w + 1;
bmx = length(b) - w + 1;

for ia = 1 : amx
    mx = min(amx - ia + 1, bmx);
    h = sum(a(ia : ia + w - 1) ~= b(1 : w));
    for ib = 1 : mx
        if h < pa(ib + ia - 1)
            pa(ib + ia - 1) = h;
            pia(ib + ia - 1) = ib;
        end
        if h < pb(ib)
            pb(ib) = h;
            pib(ib) = ib + ia - 1;
        end
        if ia + ib <= mx
            h = h - (a(ia + ib - 1) ~= b(ib)) + (a(ia + ib - 1 + w) ~= b(ib + w)); 
        end
    end
end
end
