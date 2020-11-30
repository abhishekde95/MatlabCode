den_lens_smj = [
    380 3.2 % linear extrapolation
    385 2.8 % linear extrapolation
    390 2.4
    395 1.985
    400 1.62
    405 1.308
    410 1.05
    415 0.84
    420 0.675
    425 0.557
    430 0.468
    435 0.393
    440 0.335
    445 0.29
    450 0.26
    455 0.24
    460 0.225
    465 0.215
    470 0.203
    475 0.191
    480 0.18
    485 0.168
    490 0.162
    495 0.151
    500 0.145
    505 0.133
    510 0.128
    515 0.122
    520 0.116
    525 0.11
    530 0.104
    535 0.099
    540 0.093
    545 0.087
    550 0.081
    555 0.075
    560 0.07
    565 0.064
    570 0.058
    575 0.052
    580 0.046
    585 0.041
    590 0.036
    595 0.031
    600 0.028
    605 0.024
    610 0.021
    615 0.017
    620 0.014
    625 0.012
    630 0.009
    635 0.007
    640 0.005
    645 0.003
    650 0.002
    655 0.001
    (660:5:780)' zeros(25,1)];

% valid between Mar 5 2013 - Oct 4 2014 (Dell 4)
rgbmat = [0.838070055625675,0.289365615659916,0.857826812321819,0.269608858963772;
    0.871068406039999,0.0663076893912883,0.373489014042769,0.563887081388518;
    0.452986415118669,0.539144322966528,0.498781563965455,0.493349174119741];
% only care about the chromatic directions
rgbmat = rgbmat(:,3:4);

stro = nex2stro(findfile('N062313002.nex', '~/Dropbox/Horwitz Lab/'));
mon_spd = reshape(stro.sum.exptParams.monspd, [], 3);
orig_funds = reshape(stro.sum.exptParams.fundamentals, [], 3)';

lens = den_lens_smj(:,2)./1.16;  % Dividing by 1.16 to make it "open pupil". Thats what SMJ did.
lenstransmittance = 1./(10.^(lens*1.28));
fund10 = orig_funds./repmat(max(orig_funds,[],2),1,size(orig_funds,2));
absorptance = fund10'./repmat(lenstransmittance,1,3);
absorptance = absorptance./repmat(max(absorptance),81,1);

% ****************
% Here's the value to change the degree of lens density
% 0 = no lens density
% 1.28 = human lens density
lensmults = linspace(.2,2);
lm_diffs = zeros(size(lensmults));
for k = 1:length(lensmults)
    % Now creating a NEW set of fundamentals starting with the action spectra
    % Change parameters here
    lenstransmittance = 1./(10.^(lens*lensmults(k)));
    funds = absorptance .*repmat(lenstransmittance,1,3);
    funds = funds./repmat(max(funds),81,1);
    
    mon_spd_rescaled = SplineRaw(linspace(380,780,size(mon_spd,1))', ...
        mon_spd, linspace(380,780,size(funds,1))');
    Mnew = funds'*mon_spd_rescaled;
    mulms = Mnew*stro.sum.exptParams.bkgndrgb;
    stimlms = (Mnew*rgbmat)';
    new_cc = bsxfun(@rdivide, bsxfun(@minus, stimlms, mulms'), mulms');
    lm_diffs(k) = new_cc(1,1)-new_cc(1,2);
end

figure; plot(lensmults, lm_diffs, '.')
