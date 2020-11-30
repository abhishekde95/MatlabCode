%Probably necessary but not sure - at least called by nex2db (if not by more
%functions) to convert calendar string dates into numbers
function recDate = recDateFromCal(recordingDates)
if strcmp(recordingDates, '*')
    recDate = '*';
    return;
elseif size(recordingDates,2) > 11
    rd1 = recordingDates(2:end-1);
    rd2 = strsplit(rd1, {' ', ',', ''''});
    rdate = rd2(2:end-1);
else
    rdate = {recordingDates};
end
    function m = dateconv(month)
        switch month
            case 'Jan'
                m = '01';
            case 'Feb'
                m = '02';
            case 'Mar'
                m = '03';
            case 'Apr'
                m = '04';
            case 'May'
                m = '05';
            case 'Jun'
                m = '06';
            case 'Jul'
                m = '07';
            case 'Aug'
                m = '08';
            case 'Sep'
                m = '09';
            case 'Oct'
                m = '10';
            case 'Nov'
                m = '11';
            case 'Dec'
                m = '12';
            otherwise
                m = 'err';
        end
    end
recCell = {};
for i = 1:length(rdate)
    nsd = rdate{i};
    month = nsd(4:end-5);
    mm = dateconv(month);
    if strcmp(mm, 'err')
        disp('month error');
        keyboard;
        break;
    else
        yyyy = (nsd(end-3:end));
        dd = nsd(1:2);
        recCell{i} = sprintf('%s-%s-%s',yyyy,mm,dd);
    end
end
sloppyRecDate = sprintf('''%s'',', recCell{:});
recDate = sloppyRecDate(1:end-1);
end