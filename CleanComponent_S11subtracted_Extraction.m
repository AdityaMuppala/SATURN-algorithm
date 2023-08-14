function [ShatClean,RClean,S11_prediction,brightnessFunc] = CleanComponent_S11subtracted_Extraction(S11_subtracted, RTargettemp, thetaTarget,phiTarget,ShatMaxdB,fnpts, M, frange, fstart, fstop, k, rho_n, phi_n, z_n)

%% GO model for S11 prediction

% Finding Rtarget inside range bin (Make sure to sweep around 2 lambda
% since RR is 2.5 lambda
lambda_10GGHz = physconst('Lightspeed')/(1e10);
Radjust_range = -10e-3:0.1e-3:10e-3; %Range inside which to look for target
brightnessFunc = zeros(1,length(Radjust_range));
ShatClean1D = zeros(1,length(Radjust_range));
S11_prediction_3D = zeros(fnpts,M,length(Radjust_range));
thetaSpan = 1;
thetaStep = 0.1;
phiStep = 5;
for Radjust_ind = 1:length(Radjust_range)
    Radjust = Radjust_range(Radjust_ind);
    RTarget2 = RTargettemp+Radjust;
    xpm = RTarget2.*sind(thetaTarget).*cosd(phiTarget)-rho_n.*cosd(phi_n);
    ypm = RTarget2.*sind(thetaTarget).*sind(phiTarget)-rho_n.*sind(phi_n);
    zpm = RTarget2.*cosd(thetaTarget)-z_n;
    drpm = 2*sqrt(xpm.^2+ypm.^2+zpm.^2);

    GO_2D_phase = zeros(length(frange),M); GO_2D_phase2 = zeros(length(frange),M);
    for freqpos = 1:length(frange)
        fdes = 1e9*frange(freqpos);
        GO_2D_phase(freqpos,:) = (2*pi*fdes/physconst('LightSpeed')).*(drpm);
    end

    for elementnum = 1:M
        wrappedPhase = mod(GO_2D_phase(:,elementnum),2*pi);
        unwrappedPhase = rad2deg(unwrap(wrappedPhase));
        GO_2D_phase2(:,elementnum) = deg2rad(unwrappedPhase);
    end

    % Cable loss

    Shat_prime = 10^(ShatMaxdB/20);
    ScaleFactor = fnpts*M/((4*pi*RTarget2)^2); 
    ShatClean1D(Radjust_ind) = Shat_prime/ScaleFactor; %The scale factor removes the amplification in signal from fnpts and M

    % S11 Prediction
    rpm_2D = repmat(drpm/2,[fnpts,1]);
    S11_prediction = exp(1i*GO_2D_phase2).*ShatClean1D(Radjust_ind)./((4*pi*rpm_2D).^2);
    S11_prediction_3D(:,:,Radjust_ind) = S11_prediction;
    S11_diff = S11_subtracted - S11_prediction;

    %%
    Shat_2D = Shat3D_from_S11_sectioned(S11_diff,RTarget2,thetaSpan,thetaStep,...
        phiStep,thetaTarget,phiTarget,k,fnpts,rho_n,phi_n,z_n,M);

    brightnessFunc(Radjust_ind) = sum(sum(abs(Shat_2D)));

end

[~,posMin] = min(brightnessFunc);
S11_prediction = S11_prediction_3D(:,:,posMin);
RClean = RTargettemp+Radjust_range(posMin);
ShatClean = ShatClean1D(posMin);

end

