function drive  = getcumdrive3(Subunit1wts,Subunit2wts,lightS1,lightS2,mode)
N = size(lightS1,2);
w1 = 1;
w2 = 1;
if mode == 1
%     keyboard;
    driveS1 = repmat(Subunit1wts,[1 N]).*lightS1; 
    driveS2 = repmat(Subunit2wts,[1 N]).*lightS2; 
    driveS1 = sum(driveS1,1);
    driveS2 = sum(driveS2,1);
    drive = w1*driveS1 + w2*driveS2; % linear spatial summation and linear cone-signal combination
elseif mode == 2
    driveS1 = repmat(Subunit1wts,[1 N]).*(lightS1.^2);
    driveS2 = repmat(Subunit2wts,[1 N]).*(lightS2.^2);
    driveS1 = sum(driveS1,1);
    driveS2 = sum(driveS2,1);
    drive = w1*max(0,driveS1.^2) + w2*max(0,driveS2.^2); % linear spatial summation and quadratic cone-signal combination
elseif mode == 3
    driveS1 = repmat(Subunit1wts,[1 N]).*lightS1; 
    driveS2 = repmat(Subunit2wts,[1 N]).*lightS2; 
    driveS1 = sum(driveS1,1);
    driveS2 = sum(driveS2,1);
    drive = w1*(driveS1.^2) + w2*(driveS2.^2); % broader-than-linear spatial summation and linear cone-signal combination
elseif mode == 4
    driveS1 = repmat(Subunit1wts,[1 N]).*(lightS1.^2);
    driveS2 = repmat(Subunit2wts,[1 N]).*(lightS2.^2);
    driveS1 = sum(driveS1,1);
    driveS2 = sum(driveS2,1);
    drive = w1*(driveS1.^2) + w2*(driveS2.^2); % broader-than-linear spatial summation and quadratic cone-signal combination
elseif mode == 5
    driveS1 = repmat(Subunit1wts,[1 N]).*lightS1; 
    driveS2 = repmat(Subunit2wts,[1 N]).*lightS2; 
    driveS1 = sum(driveS1,1);
    driveS2 = sum(driveS2,1);
    drive = w1*(driveS1.^2) - w2*(driveS2.^2); % narrower-than-linear spatial summation and linear cone-signal combination
elseif mode == 6
    driveS1 = repmat(Subunit1wts,[1 N]).*(lightS1.^2);
    driveS2 = repmat(Subunit2wts,[1 N]).*(lightS2.^2);
    driveS1 = sum(driveS1,1);
    driveS2 = sum(driveS2,1);
    drive = w1*(driveS1.^2) - w2*(driveS2.^2); % narrower-than-linear spatial summation and quadratic cone-signal combination
end

end

