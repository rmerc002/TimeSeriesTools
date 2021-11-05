function [PAA] = PAA_updated(timeseries, numcoeffs)
% This is not a very good approach. I do not recommend it. I recommend
% truncating sections that do not fill PAA sections, because it maintains
% correct correspondence between PAA output and time series input. This
% does __not__. 

assert(isvector(timeseries));
N = length(timeseries);
if N < 1 || N < numcoeffs
    error('invalid input');
end

per_section = (N / numcoeffs);
if floor(per_section) == per_section
    PAA = movmean(timeseries, [0 per_section-1], 'Endpoints', 'discard');
else
    PAA = NaN(numcoeffs, 1);
    paa_index = 1;
    ts_index = 1;
    carry = 1;
    while ts_index <= N && paa_index <= numcoeffs
        if carry > 0
            if carry > 1
                error('invalid carry, something is wrong');
            end
            p = carry * timeseries(ts_index);
            ts_index = ts_index + 1;
            post_carry = per_section - carry;
            full = floor(post_carry);
            fractional = post_carry - full;
        else
            p = 0;
            full = floor(per_section);
            fractional = per_section - full;
        end
        if ts_index + full - 1 > N
            warning('possible loss of precision in determining intervals');
            % just run it to the end. It may not be correct. Floating point
            % arithmetic for this kind of thing is a bit tenuous and not
            % always well conditioned
            
            % This isn't the only way it can fail. This just happens to be
            % easily detectable. Notice, skips fractional component, which
            % is definitely corrupted by roundoff if you're hitting this
            % point.
            full = N - ts_index + 1;
            p = p + sum(timeseries(ts_index : ts_index + full - 1));
            PAA(paa_index) = p;
            break;
        end
        ts_index = ts_index + full;
        if ts_index <= N && fractional > 0
            p = p + fractional * timeseries(ts_index);
            carry = 1 - fractional; 
        else
            if fractional < 0
                % again, not the only way it can fail, just simple and
                % obvious
                warning('loss of precision in carry value due to rounding');
            end
            carry = 0;
        end
        PAA(paa_index) = p;
        paa_index = paa_index + 1;
    end
    PAA = PAA ./ per_section;
end
end
