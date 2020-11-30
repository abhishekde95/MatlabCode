function [M,bkgndrgb] = M_from_calibration(cal)

if ischar(cal)
    cal = load(cal);
    cal = cal.cals{end};
elseif iscell(cal)
    cal = cal{end};
end

funds = load('T_cones_smj10');
funds = funds.T_cones_smj10;

bkgndRGB = round(255*cal.bgColor)';
bkgndrgb = [cal.gammaTable(bkgndRGB(1)+1,1)
            cal.gammaTable(bkgndRGB(2)+1,2)
            cal.gammaTable(bkgndRGB(3)+1,3)];

spd = SplineSpd(linspace(380, 780, size(cal.P_device, 1))', cal.P_device, ...
    linspace(380,780,size(funds,2))');
M = funds*spd;
