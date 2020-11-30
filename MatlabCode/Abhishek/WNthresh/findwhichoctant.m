function whichoctants = findwhichoctant(lmsvals)
% Intakes a set of LMS points and gives you back the octants where the points lie 
uniquepts = (unique(sign(lmsvals),'rows')+1)/2;
whichoctants = uniquepts*[4; 2; 1];
end