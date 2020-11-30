function print_stimuli(stimuli)
if size(stimuli, 2) == 2
    stimuli = [cos(stimuli(:,2)) sin(stimuli(:,2)) 10.^stimuli(:,1)];
end

[L,M] = arrayfun(@(x) rat(x,5e-2), abs(stimuli(:,1)./stimuli(:,2)));
L = L.*sign(stimuli(:,1));
M = M.*sign(stimuli(:,2));
fprintf('\n');
fprintf('#%d -> (% .5f, % .5f) ~ (% 4dL : % 4dM) @ %5.2f Hz\n', [(1:size(stimuli,1))' stimuli(:,1:2) L M stimuli(:,3)]');
fprintf('\n');
