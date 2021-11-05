function [pa, pia] = mpx_strings(a,minlag, w)

% this performs exact comparison
% this is intended to replace Ryan's silly dna to time series mapping,
% which I have referred to as nonsense on a daily basis. Consider this script
% my response to "put up or shut up"

if ~isvector(a) 
    % If Eamonn is unintentionally reading this, Ryan thinks I resemble
    % Rick from Rick and Morty. He's like a more degenerate sarcastic version 
    % of the doc from Back To the Future, and in this scenario, Ryan is Michael J Fox 
    % as opposed to Damian Rice.
    error('I am disappointed by this input, Morty');
end

istransposeda = false;

if isrow(a)
    a = transpose(a);
    istransposeda = true;
end


pa = repmat(w, length(a) - w + 1, 1);
pia = NaN(length(a) - w + 1, 1);

[pa, pia] = solver(a, pa, pia, minlag, w);

% scaling to the typical range of z-normalized euclidean distance
pa = sqrt(2 * pa);
if istransposeda
    pa = transpose(pa);
    pia = transpose(pia);
end

end

function [pa, pia] = solver(a, pa, pia, minlag, w)
amx = length(a) - w + 1;

for ia = 1 : amx
    mx = amx - ia + 1;
    h = sum(a(ia : ia + w - 1) ~= a(1 : w));
    for ib = 1 : mx
        if h < pa(ib + ia - 1) && minlag < abs(ia - 1)
            pa(ib + ia - 1) = h;
            pia(ib + ia - 1) = ib;
        end
        if h < pa(ib) && minlag < abs(ia - 1)
            pa(ib) = h;
            pia(ib) = ib + ia - 1;
        end
        if ia + ib <= mx
            h = h - (a(ia + ib - 1) ~= a(ib)) + (a(ia + ib - 1 + w) ~= a(ib + w)); 
        end
    end
end
end
