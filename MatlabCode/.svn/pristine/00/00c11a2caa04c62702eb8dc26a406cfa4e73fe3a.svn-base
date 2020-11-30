%%%%%%%%%%%%%%%%%%%%%%%%%
% Integrate this with GetVals?
function [out, p] = GetValsIfPossible(codes, p, type)

VALOFFSET = 4000;
out = [];

try
    if (length(codes) == 2 && codes(1) == codes(2))  % Bookending with double codes
        % GDLH 6/12/09 For backward compatibility assuming anything
        % bookended is double until otherwise specified.
        if (nargin < 3)
            type = 'double';
        end
        L = (p.events == codes(1));
        if (sum(L) < 2)
            return;
        else
            idxs = find(L);
            out = p.events(idxs(end-1)+1:idxs(end)-1)-VALOFFSET;
            
            % GDLH 1/21/16
            % Removing negative numbers from bookended code drops *only* if
            % they are not ints. Patrick drops bookended negative ints
            % and the code below was removing them. This is a bit of hack
            % anyway - why are we ever getting negative bytes at all?
            if (~strcmp(type,'int'))
                Lbad = out < 0 | out > 255;
                if any(Lbad)
                    warning(['Trimmed out the following erroneous "bytes" for code %d:\n\t' ...
                        sprintf('%d ', out(Lbad)') '\n'], codes(1));
                    out(Lbad) = [];
                end
            end
            out = uint2num(out,type);
            p.lastprocessed_t = p.times(max(idxs));
        end
    else   % Single codes or lists of codes
        if (nargin < 3)
            type = 'int';
        end
        [L, idxs] = ismember(codes, p.events);
        if (any(idxs == 0))
            return;
        end
        if (strcmpi(type,'int'))
            if (max(idxs)+1 <= length(p.events))
                out = p.events(idxs+1)-VALOFFSET;
                p.lastprocessed_t = p.times(max(idxs)+1);
            end
        elseif (strcmpi(type,'long'))
            if (max(idxs)+4 <= length(p.events))
                out = p.events(idxs+[1:4])-VALOFFSET;
                out(out == -VALOFFSET) = []; % GDLH added 1/3/15
                out = [2^0 2^8 2^16 2^24]*out;
                p.lastprocessed_t = p.times(max(idxs)+4);
            end
        elseif (strcmpi(type,'float'))
            if (max(idxs)+4 <= length(p.events))
                out = p.events(idxs+[1:4])-VALOFFSET;
                out(out == -VALOFFSET) = []; % GDLH added 1/3/15
                out = uint2num(out,'float');
                p.lastprocessed_t = p.times(max(idxs)+4);
            end
        elseif (strcmpi(type,'double'))
            if (max(idxs)+8 <= length(p.events))
                out = p.events(idxs+[1:8])-VALOFFSET;
                out(out == -VALOFFSET) = []; % GDLH added 1/3/15
                out = codes2num(out);
                p.lastprocessed_t = p.times(max(idxs)+8);
            end
        end
    end
catch ME
    disp(getReport(ME));
    keyboard
end
