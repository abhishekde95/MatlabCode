function AUC = calcAUC(prediction,response)
Hits = sum(prediction == 1 & response == 1);
Miss = sum(prediction == 0 & response == 1);
FA = sum(prediction == 0 & response == 1);
CR = sum(prediction == 0 & response == 0);

pHits = Hits/(Hits+Miss);
pFA = FA/(FA+CR);
AUC = 0.5*(1+pHits-pFA);
end

