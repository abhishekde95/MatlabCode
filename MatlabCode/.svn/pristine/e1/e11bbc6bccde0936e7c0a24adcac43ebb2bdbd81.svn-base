function nextX = computeNextStim_WpersistentVariables(x, r, stimDim, support)

persistent dbArray;

if isempty(dbArray)
    dbArray.x = x;
    dbArray.r = r;
else
    dbArray.x = [dbArray.x; x];
    dbArray.r = [dbArray.r; r];
end

nextX = rand;