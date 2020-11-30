function FR = calcresponse(Subunit1wts,Subunit2wts,lightS1,lightS2,mode)
max_FR = 100;
if mode == 1
    driveS1 = Subunit1wts'*lightS1;
    driveS2 = Subunit2wts'*lightS2;
    drive = [1 0.5]*[driveS1 driveS2]'; % linear spatial summation and linear cone-signal combination
elseif mode == 2
    driveS1 = Subunit1wts'*(lightS1.^2);
    driveS2 = Subunit2wts'*(lightS2.^2);
    drive = [1 0.5]*[driveS1 driveS2]'; % linear spatial summation and quadratic cone-signal combination
elseif mode == 3
    driveS1 = Subunit1wts'*lightS1;
    driveS2 = Subunit2wts'*lightS2;
    drive = [1 0.5]*[driveS1.^2 driveS2.^2]'; % broader-than-linear spatial summation and linear cone-signal combination
elseif mode == 4
    driveS1 = Subunit1wts'*(lightS1.^2);
    driveS2 = Subunit2wts'*(lightS2.^2);
    drive = [1 0.5]*[driveS1.^2 driveS2.^2]'; % broader-than-linear spatial summation and quadratic cone-signal combination
elseif mode == 5
    driveS1 = Subunit1wts'*lightS1;
    driveS2 = Subunit2wts'*lightS2;
    drive = [1 -0.5]*[driveS1.^2 driveS2.^2]'; % narrower-than-linear spatial summation and linear cone-signal combination
elseif mode == 6
    driveS1 = Subunit1wts'*(lightS1.^2);
    driveS2 = Subunit2wts'*(lightS2.^2);
    drive = [1 -0.5]*[driveS1.^2 driveS2.^2]'; % narrower-than-linear spatial summation and quadratic cone-signal combination
end
FR = max_FR * 1./(1+exp(-1*(drive-1))); % Spike generation mechanism: sigmoidal non-linearity 
end

