% The goal of this function is to make a small space-delimited string
% representation of a vector of integers. We represent each integer as the raw
% bytes (in hexadecimal) of its magnitude prepended with '-' if the value was
% negative.
%
% This function properly handles all representable integers of integral data
% types (uintX and intX with X = 8,16,32,64) and all double precision integers
% <= 2^52. If the input contains a double precision integer > 2^52, this
% function issues a warning about potential inaccuracies. If you need to
% represent integers this large, use the int64 or uint64 data type. (Compare
% what happens when the input to this function is 2^53+1 and int64(2)^53+1).

function out = notmydec2hex(in)
nonnegative = in >= 0;
if isinteger(in)
    intclass = class(in);
    if intclass(1) ~= 'u'
        intmins = in == intmin(intclass);
        if any(intmins(:))
            % the most negative signed integer of n bits is -2^(n-1), but we
            % can't fit its magnitude into the largest signed integer 2^(n-1)-1.
            % after `abs`, the result is truncated to 2^(n-1)-1, so cast the
            % result to unsigned and add 1.
            in = feval(['u' intclass], abs(in));
            in(intmins) = in(intmins)+1;
        else
            in = abs(in);
        end
    end
else
    if any(in(:) > 1/eps)
        warning(message('MATLAB:dec2hex:TooLargeArg')); % hijack dec2hex's warning message
    end
    in = abs(in);
end

out = sprintf(' -%x', in);
minuses = find(out == '-');
% remove leading space too
out([1 minuses(nonnegative)]) = [];
