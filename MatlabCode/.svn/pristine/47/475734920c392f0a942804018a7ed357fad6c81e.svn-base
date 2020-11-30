function d = calcdprime(Hit,FA,stimpresent,stimabsent)
% A function for calculating d-prime

% Assuming FA never reaches stimabsent
if Hit==0
    d = norminv(1-(0.5./stimpresent))-norminv(max([0.5./stimabsent FA./stimabsent]));
else
    d = norminv(min([1-(0.5./stimpresent) Hit./stimpresent]))-norminv(max([0.5./stimabsent FA./stimabsent]));
end
end

