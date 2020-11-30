% Exclude points on a surface whose theta or Z components lie outside the
% specified bounds and scale r. `restrict` = [lower r scale (<= 1), upper r
% scale (>= 1), lower theta bound (>= 0), upper theta bound (<= pi), lower Z
% bound, upper Z bound]. This script assumes that theta wraps around from pi to
% 0.

function out = restrict_samples(subsampled, restrict)
if ~any(restrict)
    out = subsampled;
    return
end

[sub_th, sub_r, sub_z] = cart2pol(subsampled(:,1), subsampled(:,2), subsampled(:,3));

r_scale = restrict(1:2); th_bounds = restrict(3:4); z_bounds = restrict(5:6);

purge = false(size(sub_th));

if any(th_bounds)
    purge(mod(sub_th, pi) < th_bounds(1) | mod(sub_th, pi) > th_bounds(2)) = true;
end

if any(z_bounds)
    purge(sub_z < z_bounds(1) | sub_z > z_bounds(2)) = true;
end

sub_th(purge) = []; sub_r(purge) = []; sub_z(purge) = [];

if any(r_scale)
    a = log10(r_scale(1));
    b = log10(r_scale(2));
    sub_r = sub_r .* 10.^(a + (b-a) .* rand(size(sub_r)));
end

[newX, newY, newZ] = pol2cart(sub_th, sub_r, sub_z);
out = [newX newY newZ];
