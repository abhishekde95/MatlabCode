function new_mask_changes = user_defined_mask_changes(mask_changes,ind)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if (nargin < 1)
    len = size(mask_changes,2);
else 
    len = ind;
end

new_mask_changes = zeros(size(mask_changes,1),numel(len));
for i = 1:numel(len)
    new_mask_changes(:,i) = mask_changes(:,len(i));
end


end

