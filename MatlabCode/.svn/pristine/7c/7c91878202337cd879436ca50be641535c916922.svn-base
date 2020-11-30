function Spectrum=GeneratePhotoreceptorSpectrum(Type,Wavelength)

% Spectrum=GeneratePhotoreceptorSpectrum(Type,Wavelength)
% Type options: 'LCone','MCone','SCone','Rod'
% Generates Normalized Photoreceptor Spectral Sensitivity Curve relying on
% fits derived from Baylor's et al. suction recordings (1985, 1987) on
% macaque photoreceptors and fitted according to the Govadorvskii nomogram
% (2000). Fitting procedure can be checked in FittingGovadorvskii.m
% Created Dec_2010 Angueyra

switch Type
    case 'LCone'
        %values derived from fitting procedure
        LambdaMax=565.2836;
        BetaBandScaling=0.1460;
        Spectrum=GovadorvskiiTemplateA1([LambdaMax,BetaBandScaling],Wavelength);
    case 'MCone'
        %values derived from fitting procedure
        LambdaMax=534.3201;
        BetaBandScaling=0.1797;
        Spectrum=GovadorvskiiTemplateA1([LambdaMax,BetaBandScaling],Wavelength);        
    case 'SCone'
        %values derived from fitting procedure
        LambdaMax=430.4053;
        BetaBandScaling=0.5210;
        Spectrum=GovadorvskiiTemplateA1([LambdaMax,BetaBandScaling],Wavelength);  
    case 'Rod'
        %values derived from fitting procedure
        LambdaMax=493.3022;
        BetaBandScaling=0.1233;
        Spectrum=GovadorvskiiTemplateA1([LambdaMax,BetaBandScaling],Wavelength);
end