% This function returns true if the machine is hooked up to at least one
% *Pixx device. The result is stored for fast retrieval on subsequent
% calls. Run `clear all` to remove the stored result.

function tf = IsVPixx()
persistent rslt
if isempty(rslt)
    if exist('Datapixx') ~= 3 %#ok<EXIST> % no Datapixx mex file
        warning('No Datapixx mex file found to drive VPixx devices!');
        rslt = 0;
    elseif Datapixx('IsReady') % the device is open; don't close it at the end
        rslt = FoundPixxDevice();
    else
        try % mex file throws an error when there's no VPixx device found
            % evalc consumes the error message printed by the mex file
            % this way users aren't confusedby the harmless message
            [nil,is_open] = evalc('Datapixx(''Open'')');
            if is_open
                rslt = FoundPixxDevice();
                Datapixx('Close');
            end
        catch
            rslt = 0;
        end
    end
end
tf = rslt;

function rslt = FoundPixxDevice()
cmds = {'IsPropixxCtrl' 'IsPropixx' 'IsDatapixx' 'IsDatapixx2' 'IsViewpixx' 'IsViewpixx3D'};
rslt = 0;
for k = 1:length(cmds)
    if Datapixx(cmds{k})
        fprintf('** %s: Connected to %s **\n\n', mfilename, cmds{k}(3:end));
        rslt = 1;
        break
    end
end