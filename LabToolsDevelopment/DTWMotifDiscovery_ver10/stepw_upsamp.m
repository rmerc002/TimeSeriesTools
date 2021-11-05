function [stepw] = stepw_upsamp(ts, upsamp_len)

% Takes a time series and an upsampling length
% It will make a stepwise function approximation of this

if ~isvector(ts) || upsamp_len < length(ts)
    error('invalid parameters');
elseif length(ts) == upsamp_len
    stepw = ts;
    return;
end

upsamp_factor = upsamp_len / length(ts);
fractional = upsamp_factor - floor(upsamp_factor);

if fractional == 0
    stepw = zeros(upsamp_len, 1);
    for i = 1 : length(ts)
        stepw((i - 1) * upsamp_factor + 1 : i * upsamp_factor) = ts(i);
    end
else
    if fractional < 0.05
        per = 20;
    else
        per = round(1 / fractional);
    end
    % don't compensate more than a factor of 20. There might be a better
    % way to do this, but this avoids blowing up memory on very long arrays
    % this could be made into a parameter
    per = per * upsamp_len;
    stepw = zeros(per * length(ts), 1);
    for i = 1 : length(ts)
        stepw((i - 1) * per + 1 : per * i) = ts(i);
    end
end

stepw = interp1(linspace(1, upsamp_len, length(stepw)), stepw, (1 : upsamp_len)');

end

