%DAT2NUM
%
%   EXAMPLE: nums = dat2num(vals, codes, offset, dataType [, bookend]);
%
% Takes the values from a PLEXON data file and extracts the appropriate number
% of values following each instance of a 7/8000 level code and converts them to
% a numeric representation for MATLAB. 'vals' should be the vector of event
% values and ecodes. 'codes' should be the vector of ecodes to find in 'vals'.
% 'offset' should be the offset value added during transmission from REX to the
% PLEXON data stream. 'dataType' should be a string defining how REX represented
% the data represented by 'codes' before they were droped into the PLEXON data
% stream. 'dataType' can be u/short, u/int, u/long, u/llong, float, or double
% (case insensitive). 'bookend' should be set to 1 if the values to be converted
% are bookended with ecodes. If 'bookend' is set to 1 and there are not exactly
% 2 occurrences of the ecode then this function skips to the next ecode or
% returns an empty cell array if called with only one ecode. 'nums' is a cell
% array with the same number of elements as 'code' (i.e., one vector for each
% ecode in 'code').

% cah  01/15/08
% zalb 10/10/13 cleaned up a bit for a small performance boost
% zalb 06/20/14 adding more data types

function nums = dat2num(vals, codes, OFFSET, dataType, bookend)
lowCutoff = 2000; % vals below 2000 shouldn't occur via the low priority code dropping
% b/c they'll get confused with real 1000 level ecodes.

if any(strcmpi(dataType, {'int' 'char'}))
    bytes_wide = 1;
elseif any(strcmpi(dataType, {'short' 'ushort'}));
    bytes_wide = 2;
elseif any(strcmpi(dataType, {'float' 'long' 'ulong' 'uint'}))
    bytes_wide = 4;
elseif any(strcmpi(dataType, {'double' 'llong' 'ullong'}))
    bytes_wide = 8;
else
    error('Unspecified or inadequate data type specifier.');
end

do_bookending = nargin >= 5 && bookend;

% find out of range codes (in this case codes that could conflict with 1000
% level codes) and remove them.
vals = vals(vals >= lowCutoff);
nums = cell(size(codes));
for ii = 1:length(codes)
    codeIdxs = find(vals == codes(ii))';
    numCodes = length(codeIdxs);
    
    % ignore one bad code and continue to the rest
    if numCodes == 0
        nums{ii} = nan;
        continue
    end
    
    % Now extract the appropriate number of values following (or between)
    % occurrences of the ecode and subtract the offset value.
    if do_bookending
        if numCodes ~= 2
            error('\n *** Bookending expected 2 <%d> ecodes but found %d ***\n', codes(ii), numCodes);
        end
        bytes = vals(codeIdxs(1)+1 : codeIdxs(2)-1) - OFFSET;
    else
        % generate the indices of the bytes (they're in the positions after the `codeIdxs`)
        ind = zeros(bytes_wide, numCodes);
        for jj = 1:bytes_wide
            ind(jj,:) = ind(jj,:) + codeIdxs + jj;
        end
        bytes = vals(ind) - OFFSET;
    end
    
    % now convert to a numeric representation that makes sense for MATLAB (i.e.,
    % characters, doubles, or u/int64s)
    nums{ii} = uint2num(bytes, dataType);
end
