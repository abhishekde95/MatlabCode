
function [output]=xml2mat(thisxml)


[p xmlname] = fileparts(thisxml);
matname = [xmlname '.mat'];
warning off;
z = myReadExpoXML(thisxml);
matfilename =  [thisxml(1:end-3) 'mat'];
warning on;

tick2msec = GetConversionFactor(z, 'ticks', 'msec');
pix2deg = GetConversionFactor(z, 'pix', 'deg');

% get main variables

% get the slot number that corresponds to the 'eval' slot in the expo program
validslot = GetSlots(z,'eval');

% get the pass numbers that the program went over the target slot
validpasses = GetPasses(z,validslot,0);

% go through those passes and for specific routines in the slot, read off variables of interest
var1 = GetEvents(z, validpasses,'Op Variable','coherence', 0, 'Destination');
var2 = GetEvents(z, validpasses,'Op Variable','eval', 0, 'Source 1');
var3 = GetEvents(z, validpasses,'Op Variable','eval', 0, 'Source 2');
var4 = GetEvents(z, validpasses,'Op Variable','responsetime', 0, 'Destination');
%var3 = GetEvents(z, validpasses,'Op Variable','routine_label3', 0, 'Destination');
output =[var1(:) var2(:) var3(:) var4(:)];
%save var_vec var1 var2 var3
