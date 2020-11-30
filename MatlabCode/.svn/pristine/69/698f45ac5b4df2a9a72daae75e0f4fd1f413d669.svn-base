% SENDTOREX    Send the contents of variables to REX via UDP
%    SENDTOREX(udpCom, toRex, dataType, message) encodes `toRex` as a string and
%    bundles it into a larger message to send over UDP. 'udpCom' is a structure
%    that specifies the UDP socket details, like the port number (.port) and
%    REX's local IP address (.rexip). 'message' should be an exact copy of the
%    entire function call string before execution.
%
% For example:
%    >> toRex = rand(20);
%    >> message = 'SENDTOREX(udpCom, toRex, ''double'', message)';
%    >> SENDTOREX(udpCom, toRex, 'double', message);
%
% CAH  11/20/07
% ZALB 06/20/14 more data types

function notsendToRex(udpCom, toRex, dataType, msg)
REXBUFF = 8000; % The size of the UDP read buffer on the REX machine in bytes

% Turn the numbers that REX requests into a long string of space delimited
% hexadecimal numbers. Do this differently for integers and doubles.
if strcmpi(dataType, 'double')
    calElement = sprintf(' %bx', toRex); % front pad with a space!
    calElement(1) = [];
elseif strcmpi(dataType, 'float')
    calElement = sprintf(' %tx', toRex);
    calElement(1) = [];
elseif any(strcmpi(dataType, ...
        {'char' 'ushort' 'short' 'uint' 'int' 'integer' ...
        'long' 'ulong' 'ullong' 'uint64' 'llong' 'int64'}))
    calElement = notmydec2hex(toRex);
else
    error('Unspecified or inadequate data type specifier "%s"', dataType);
end

% The first UDP message will be from REX's main state set requesting a number or
% list of numbers. Echo this request (i.e, the EXACT same string that REX sent)
% so that REX can call the appropriate subfunction from with in the UDP
% "listening" state set. If the requested data will fit with this "echo" go
% ahead and send them together. This will force the REX to use the
% small-transfer protocol instead of the large-batch transfer protocol.
maxMsgLength = REXBUFF-1; % buffer size minus \0
numElements = notmydec2hex(numel(toRex));
totalLengthInChar = length(calElement) + length(msg) + length(numElements) + 5; % four for both >>'s one for \0

if totalLengthInChar < maxMsgLength
    pnet(udpCom.sock, 'write', [msg '>>' numElements '>>' calElement]);
    pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
    return % exit here for small-batches
else
    pnet(udpCom.sock, 'write', [msg '>>' numElements '>>']);
    pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
    pause(0.005); % pause to let REX catch up.
end

% Send the pieces in maximally sized chunks
starts = 1:maxMsgLength:length(calElement);
stops = [starts(2:end)-1 length(calElement)];
for a = [starts; stops]
    chunk = calElement(a(1):a(2));
    pnet(udpCom.sock, 'write', chunk)
    pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
    pause(0.005); % avoid REX buffer overflow issues (necessary?)
end

% If the length of calElement is greater than the size of the REX read buffer,
% then REX will have to concatenate the transmitted set of hex numbers into a
% character array in a while loop. In order to break out of this while loop send
% a message that indicates that the transmission has been sent in its entirety.
pnet(udpCom.sock, 'write', 'transferComplete'); % must be identical to that used by the transmission functions on REX
pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
