function new_idx = find_appro_index(gl)
new_idx = [];

N = numel(gl.orig_color_ID);
for jj = 1: N
    flag = 0;
    ii = 0;
    while (flag == 0)
        ii = ii + 1;
        flag = strcmp(gl.orig_color_ID(ii),gl.color_ID(jj));
    end
    new_idx = [new_idx ii];
end

end

