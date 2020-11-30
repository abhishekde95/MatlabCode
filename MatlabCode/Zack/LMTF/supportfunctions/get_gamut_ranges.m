% This function makes the supports for each QUEST function based on the current
% display's gamut. We extend the upper end to twice the maximum cone contrast to
% allow the QUEST procedure to suggest contrasts outside the physically
% allowable range. From the help text of Psychtoolbox's QuestCreate.m:
%
%   Don't restrict [standard deviation] or "range" by the limitations of what
%   you can display. Keep open the possibility that threshold may lie outside
%   the range of contrasts that you can produce, and let Quest consider all
%   possibilities.
%
% Remember that the contrast that Quest suggests is just that, a suggestion.
% Other parts of this paradigm will truncate Quest's suggested contrast to the
% edge of the gamut and update the Quest function at the presented contrast.

function [ranges,maxccs] = get_gamut_ranges(LM)
[~,maxccs] = gamut_extent(LM);
ranges = cell(size(LM,1),1);
for k = 1:length(ranges)
    ranges{k} = 1e-4:5e-4:2*maxccs(k); % gamut edge in the middle of the support
end
