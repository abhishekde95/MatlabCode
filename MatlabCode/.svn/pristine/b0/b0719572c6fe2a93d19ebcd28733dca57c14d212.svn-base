%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Go through the eventList returned by PL_GetTS and put the information 
% in the 'p' structure.
%
% Fields of the p structure are:
%   p.times             Vector containng the times of the events
%   p.events         	Vector containing the event codes
%   p.spikes            Cell array containing the spike times
%   p.processnowflag    1 = the p structure contains a full trial of events
%                       0 = does not yet contain a full trial of events
function p = ProcessEventList(p, eventList)

    CODETYPECOL = 1;
    CHANNELCOL = 2;
    EVENTCOL = 3;
    TIMECOL = 4;
    
    SPIKE = 1;
    EVENT = 4;
    
    eventList(:,TIMECOL) = eventList(:,TIMECOL)*1000;   % converting to msec
    L = eventList(:,TIMECOL) > p.lastprocessed_t;
    Levent = eventList(:,CODETYPECOL) == EVENT;
    Lspike = eventList(:,CODETYPECOL) == SPIKE & eventList(:,EVENTCOL) ~= 0; % eventList(:,EVENTCOL) = 0 for unselected waveforms
    uniquespikeidxs = unique(eventList(Lspike,[CHANNELCOL EVENTCOL]),'rows');
    p.times = [p.times; eventList(L&Levent,TIMECOL)];
    p.events = [p.events; eventList(L&Levent,EVENTCOL)];
            
    for i = 1:size(uniquespikeidxs,1)
        Lthisspike = Lspike&eventList(:,CHANNELCOL) == uniquespikeidxs(i,1) & eventList(:,EVENTCOL) == uniquespikeidxs(i,2);
        spiketimevect = eventList(L&Lthisspike,TIMECOL);
        p.spikes{uniquespikeidxs(i,2),uniquespikeidxs(i,1)} = [p.spikes{uniquespikeidxs(i,2),uniquespikeidxs(i,1)}; spiketimevect];
    end
    if (any(p.events == p.processnowcode))
       p.processnowflag = 1;
    end
end