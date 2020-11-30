function smatrix = calc_serialcorrelation(colordirchoiceidxs,colorstimpresent)
trialsids = zeros(size(colordirchoiceidxs));
trialsids(colordirchoiceidxs & colorstimpresent) = 1; % Hits
trialsids(~colordirchoiceidxs & colorstimpresent) = 2; % Misses
trialsids(colordirchoiceidxs & ~colorstimpresent) = 3; % Correct Rejects
trialsids(~colordirchoiceidxs & ~colorstimpresent) = 4; % False Alarms
smatrix = zeros(4,4);
N = numel(colordirchoiceidxs)-1;
for jj = 1:4
    for ii = 1:4
        smatrix(ii,jj) = numel(find(trialsids(find(trialsids(1:N)==ii)+1)==jj));
    end
end
smatrix = smatrix(:);
end

