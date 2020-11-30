function [err, errs] = tf_fiterr3(params, data, mode)
% [err, errs] = tf_fiterr3(params, data, mode)
%
% Used to fits a model to an entire LMTF dataset where the two xi parameters (for
% LUM and RG) are allowed to vary across spatial location, but all the
% other parameters are fixed across locations. The theta parameter, which
% controls the rotation of the funnel, is allowed to change in some
% modes.
%
% This function can operate in several modes.
% Mode 0, or 11+2n, fits the two parameters, [xi_LUM, xi_RG], independently at each
%    spatial location.
% Mode 1, or 10+3n, fits the three parameters, [xi_LUM, xi_RG, theta], independently
%    at each spatial location.
% Mode 1.1, or 10+3n, fits the three parameters, [xi_LUM, n_LUM, xi_RG], independently
%    at each spatial location.
% Mode 2, or rampy trough, fits two parameters, [xi_LUM, xi_RG], jointly across
%    spatial locations. Under this model xi_LUM = b1*r+b2*cos(2*phi)+b0
%                                        xi_RG  = a1*r+a2*cos(2*phi)+a0
% Mode 3, or tilted rampy trough, fits the two parameters, [xi_LUM, xi_RG], jointly across
%    spatial locations, and allows xi_LUM to have an asymmetry about the horizontal meridian.
%    Under this model xi_LUM = b1*r+b2*sin(2*phi)+b3*cos(2*phi)+b0
%                     xi_RG =  a1*r+a2*cos(2*phi)+a0
% Mode 4, or double tilted rampy trough, fits the two parameters, [xi_LUM, xi_RG], jointly across
%    spatial locations, and allows xi_LUM and xi_RG to have different asymmetries about the
%    horizontal meridian.
%    Under this model xi_LUM = b1*r+b2*sin(2*phi)+b3*cos(2*phi)+b0
%                     xi_RG =  a1*r+a2*sin(2*phi)+a3*cos(2*phi)+a0
% Mode 5, or yoked double tilted rampy trough, fits the two parameters, [xi_LUM, xi_RG], jointly across
%    spatial locations, and forces xi_LUM and xi_RG to have matched asymmetries about the horizontal meridian.
%    Under this model xi_LUM = b1*r+b2*sin(2*phi)+b3*cos(2*phi)+b0
%                     xi_RG =  a1*r+b2*sin(2*phi)+a3*cos(2*phi)+a0
%
% INPUT arguments: (params, data, mode)
%  params (mode 0): a ridiculously long vector. The first eleven elements are the
%        ones that are common to every retinal location:
%        [zeta_LUM, n_LUM, delta n_LUM, tau_LUM, kappa_LUM, zeta_RG, n_RG, delta n_RG, tau_RG, kappa_RG, theta]
%         After that, arguments come in pairs: [xi_LUM, xi_RG].
%         There is one [xi_LUM, xi_RG] pair for each unique retinal
%         location represented in the second input argument, data.
%  params (mode 1): As above, but the first ten elements are the
%        ones that are common to every retinal location:
%        [zeta_LUM, n_LUM, delta n_LUM, tau_LUM, kappa_LUM, zeta_RG, n_RG, delta n_RG, tau_RG, kappa_RG]
%         After that, arguments come in sets of three: [xi_LUM, xi_RG, theta]
%         and there is one set of [xi_LUM, xi_RG, theta] for each unique retinal
%         location represented in the second input argument, data.
%  params (model 1.1): As above, but the first ten elements are the
%        ones that are common to every retinal location:
%        [zeta_LUM, delta n_LUM, tau_LUM, kappa_LUM, zeta_RG, n_RG, delta n_RG, tau_RG, kappa_RG, theta]
%         After that, arguments come in sets of three: [xi_LUM, n_lum, xi_RG]
%         and there is one set of [xi_LUM, n_lum, xi_RG] for each unique retinal
%         location represented in the second input argument, data.
%  params (mode 2):  The first eleven elements are the ones that are common to every retinal location:
%         [zeta_LUM, n_LUM, delta n_LUM, tau_LUM, kappa_LUM, zeta_RG, n_RG, delta n_RG, tau_RG, kappa_RG, theta]
%         After that, come [b0, b1, b2, a0, a1, a2].
%  params (mode 3):  The first eleven elements are the ones that are common to every retinal location:
%         [zeta_LUM, n_LUM, delta n_LUM, tau_LUM, kappa_LUM, zeta_RG, n_RG, delta n_RG, tau_RG, kappa_RG, theta]
%         After that, come [b0, b1, b2, b3, a0, a1, a2].
%  params (mode 4):  The first eleven elements are the ones that are common to every retinal location:
%         [zeta_LUM, n_LUM, delta n_LUM, tau_LUM, kappa_LUM, zeta_RG, n_RG, delta n_RG, tau_RG, kappa_RG, theta]
%         After that, come [b0, b1, b2, b3, a0, a1, a2, a3].
%  params (mode 5):  The first eleven elements are the ones that are common to every retinal location:
%         [zeta_LUM, n_LUM, delta n_LUM, tau_LUM, kappa_LUM, zeta_RG, n_RG, delta n_RG, tau_RG, kappa_RG, theta]
%         After that, come [b0, b1, b2, b3, a0, a1, a3].
%
%  data: an nx6 matrix, the columns of which are [L, M, TF, Loog, RFX,
%  RFY]. RFX and RFY are assumed to be in 10ths of degrees.
%
% In data(:,[1 2]), out of gamut points should indicate the edge of the
% monitor gamut. These points contribute to the errors of the fit if the
% predicted threshold is lower than the actual threshold.
%
% OUTPUT is some reasonable measure of error
% Optional second argument is error per retinal location.
%
% See tf_fiterr2.m for details of the 2D Watson temporal contrast sensitivity model.

uniqueXYs = unique(data(:,[5 6]),'rows','stable');
nuniqueXYs = size(uniqueXYs,1);
tfs = data(:,3);
Loog = data(:,4);
rfx = data(:,5);
rfy = data(:,6);

% Input error checking
actualnumberofparameters = length(params);
switch mode
    case 0
        expectednumberofparameters = 11+2*nuniqueXYs;
    case {1, 1.1, 1.2}
        expectednumberofparameters = 10+3*nuniqueXYs;
    case 2
        expectednumberofparameters = 17;
    case {3,5}
        expectednumberofparameters = 18;
    case 4
        expectednumberofparameters = 19;
    otherwise
        error('Unknown mode: %d',mode);
end
if (actualnumberofparameters ~= expectednumberofparameters)
    error('Mode %d: expecting %d parameters, got %d.', mode, expectednumberofparameters, actualnumberofparameters);
end

err = 0;
errs = nan*ones(nuniqueXYs,1);
if (mode == 0 || mode == 1 || mode == 1.1 || mode == 1.2)
    % Setting up parameters
    if ~mode % mode 0
        ximat = reshape(params(12:end),2, nuniqueXYs);
    else
        ximat = reshape(params(11:end),3, nuniqueXYs);
    end
    for ecc_counter = 1:nuniqueXYs
        ecc = uniqueXYs(ecc_counter,:);
        Lecc = all(data(:,[5 6]) == repmat(ecc,size(data,1),1),2);
        % Just pulling out the data that's at the right retinal position (ecc).
        data_one_ecc = data(Lecc,:); % From here on out in the loop, using "data_one_ecc" which is the subset of the data at a particular retinal position.
        if ~mode %mode 0
            model = [ximat(1,ecc_counter); params(1:5); ximat(2,ecc_counter); params(6:11)];
        elseif mode == 1
            model = [ximat(1,ecc_counter); params(1:5); ximat(2,ecc_counter); params(6:10); ximat(3,ecc_counter)];
        elseif mode == 1.1
            model = [ximat(1,ecc_counter); params(1); ximat(2,ecc_counter); params(2:4); ximat(3,ecc_counter); params(5:10);];
        elseif mode == 1.2
            model = [ximat(1,ecc_counter); params(1:5); ximat(2, ecc_counter); params(6); ximat(3, ecc_counter); params(7:10);];
        else
            error(['unknown mode ',num2str(mode)]);
        end
        err = err+tf_fiterr2(model, data_one_ecc);
    end
else
    % mode2 = RT: xi_LUM = b0+b1*r+b2*r*cos(2*phi) and  xi_RG = a0+a1*r+a2*r*cos(2*phi)
    % mode3 = TRT: xi_LUM = b0+b1*r+b2*r*sin(2*phi)+b3*r*cos(2*phi)   and  xi_RG = a0+a1*r+a2*r*cos(2*phi)
    % mode4 = DTRT: xi_LUM = b0+b1*r+b2*r*sin(2*phi)+b3*r*cos(2*phi)  and  xi_RG = a0+a1*r+a2*r*sin(2*phi)+a3*r*cos(2*phi)
    % mode5 = YDTRT: xi_LUM = b0+b1*r+b2*r*sin(2*phi)+b3*r*cos(2*phi) and  xi_RG = a0+a1*r+b2*r*sin(2*phi)+a2*r*cos(2*phi)
    
    for ecc_counter = 1:nuniqueXYs
        ecc = uniqueXYs(ecc_counter,:);
        model = LMTF_global_to_local_model(params, abs(ecc(1))/10,ecc(2)/10, mode);
        Lecc = all(data(:,[5 6]) == repmat(ecc,size(data,1),1),2);
        data_one_ecc = data(Lecc,:); % From here on out in the loop, using "data_one_ecc" which is the subset of the data at a particular retinal position.
        %errs(ecc_counter) = UpdateErr(data_one_ecc,model);
        errs(ecc_counter) = tf_fiterr2(model, data_one_ecc);
    end
    err = sum(errs);
end
end
