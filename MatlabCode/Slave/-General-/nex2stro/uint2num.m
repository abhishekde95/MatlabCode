%UINT2NUM   Convert array of bytes to numbers
% NUMS = UINT2NUM(A, TYPE) takes a vector of bytes A and the desired numeric
% TYPE and returns the results as NUMS. TYPE must be one of the following
% strings: 'char', 'short', 'ushort', 'int', 'uint', 'long', 'ulong', 'llong',
% 'ullong', 'float', or 'double'. There will be an error if the number of
% elements of A is not a multiple of the number of bytes of the output data
% type. By convention, 'int' is assumed to represent small integral values and
% not bytes; in this case UINT2NUM returns the data unchanged. Use a
% different data type if you have codes (defined as value+VALOFFSET) less than
% 2000 or greater than 6999.
%
% This function will cast the output as double precision unless the requested
% data type is a 64-bit integer.

% cah  01/15/08
% cah  03/15/08 error checking for doubles less than 15
% gdlh 08/20/08 adding support for floats
% zalb 10/14/13 using typecast to reinterpret the bytes directly
% zalb 06/20/14 adding more data types

function nums = uint2num(bytes, datatype)

if strcmpi(datatype, 'int')
    nums = bytes;
    return
elseif strcmpi(datatype, 'char')
    nums = char(bytes(:)');
    return
end

from_labels = {'float' 'double' 'long' 'short' 'ushort' 'uint' 'ulong' 'llong' 'ullong'};
to_labels = {'single' 'double' 'int32' 'int16' 'uint16' 'uint32' 'uint32' 'int64' 'uint64'};

desired_datatype = to_labels(strcmpi(from_labels, datatype));

if isempty(desired_datatype)
    error('Unrecognized data type given "%s"', datatype);
end

try
    nums = typecast(uint8(bytes), desired_datatype{1});
    % don't cast 64bit ints since a double can't represent all integers > 2^53
    if ~strcmp(desired_datatype{1}(end-1:end), '64')
        nums = double(nums);
    end
catch e
    error('Unexpected number of bytes for data type %s (%s)', datatype, strtrim(e.identifier));
end
