%called by nex2db, at least, if not others. 
function [rfX, rfY, bf] = dbRFSorter(PID, stro)
%RF info is named randomly, depending on the paradigm. converts to standard
%format.
if(PID==212)
    RFX = stro.sum.exptParams.rfPosX;
    RFY = stro.sum.exptParams.rfPosY;
elseif(PID==211)
    RFX = stro.sum.exptParams.RFPosX;
    RFY = stro.sum.exptParams.RFPosY;
elseif (PID==157 || PID==213)
    RFX = stro.sum.exptParams.stim_x;
    RFY = stro.sum.exptParams.stim_y;
elseif (PID==158)
    RFX = stro.sum.exptParams.targ_x;
    RFY = stro.sum.exptParams.targ_y;
else
    RFX = stro.sum.exptParams.rf_x;
    RFY = stro.sum.exptParams.rf_y;
end
RFs = {RFX RFY};
for a = 1:length(RFs)
    if isnan(RFs{a})
        RFs{a} = 'NULL';
    elseif RFs{a} >= -1000 && RFs{a} <= 1000
        RFs{a} = int2str(RFs{a});
    else
        RFs{a} = 'Junk';
    end
end
if sum(strcmp(RFs, 'Junk')) > 0
    rfX = 'Junk';
    rfY = 'Junk';
    bf = 1;
else
    bf = 0;
    rfX = RFs{1};
    rfY = RFs{2};
end
end