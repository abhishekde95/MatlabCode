% Allow the user to edit the stimuli manually by typing into the command window
% just before beginning an experiment. I probably spent a lot of time writing
% this and it'll never be used.

function [stims,changed] = user_stim_edit(stims)
changed = false(size(stims,1),1);
reply = ' ';
while ~isempty(reply) && reply ~= 'y' && reply ~= 'n'
    reply = lower(input('Do you want to edit the stimuli above? (default = n): ', 's'));
end

if isempty(reply) || reply == 'n', return; end

reply = '';
fprintf('\nPlease enter the new stimulus components separated by spaces (L M TF)\n');
for k = 1:size(stims,1)
    valid = false;
    while ~valid
        reply = input(sprintf('New stimulus #%d (enter to skip): ', k), 's');
        [valid,new_stim] = valid_stim_str(reply);
    end
    if ~isempty(reply)
        stims(k,:) = [log10(new_stim(3)) mod(cart2pol(new_stim(1), new_stim(2)), pi)];
        changed(k) = true;
    end
end
fprintf('\nHere is the new set of stimuli:\n');
print_stimuli(stims);

function [valid,new_stim] = valid_stim_str(s)
if isempty(s)
    valid = true;
    new_stim = [];
    return
end
new_stim = str2num(s); %#ok<ST2NM>
valid = ~isempty(s) && length(new_stim) == 3;
