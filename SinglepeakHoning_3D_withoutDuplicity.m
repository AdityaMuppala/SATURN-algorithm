function [ShatMaxdB,R_center,theta_center,phi_center] = SinglepeakHoning_3D_withoutDuplicity(S11_subtemp, CleanMatrix, M, fstart, fstop, fnpts, RBigTarget, RR, k, rho_n, phi_n, z_n)

frange = linspace(fstart,fstop,fnpts);
bw = (fstop-fstart)*1e9;

%% Generating 3D Matrix with gating each bin

theta_span_array = [10 1];
theta_step_array = [1 0.1];
theta_center = 0;
phi_step_array = [1 0.2];
phi_center = 90;
R_span_array = [0.5 0.1];
R_step_array = [RR RR/3];
R_center = RBigTarget;

%Set span in theta and phi around clean components where peak search is
%avoided

for trial = 1:2
%     trial
    Rrange = R_center-R_span_array(trial)/2:R_step_array(trial):R_center+R_span_array(trial)/2;

    thetaDup = CleanMatrix(CleanMatrix(:,1)~=0,3);
    phiDup = CleanMatrix(CleanMatrix(:,1)~=0,4);
    RDup = CleanMatrix(CleanMatrix(:,1)~=0,2);

    % Image generation by phase compensation (digital beamforming)
    %%
    [Shat_3D,phi,theta,rabs] = Shat3D_from_S11_sectioned(S11_subtemp,Rrange,theta_span_array(trial),theta_step_array(trial),...
        phi_step_array(trial),theta_center,phi_center,k,fnpts,rho_n,phi_n,z_n,M);

    %Removing duplicate peaks from CleanMatrix
    mask = imregionalmax(abs(Shat_3D));
    [Shat3D_peaks,loc] = maxk(abs(Shat_3D(mask)),10);
    thetatemp = theta(mask); phitemp = phi(mask); rtemp = rabs(mask);
    thetapeaks = thetatemp(loc); phipeaks = phitemp(loc); rpeaks = rtemp(loc);

    %Removing duplicate peaks
    dupMatrix = [RDup.*sind(thetaDup).*cosd(phiDup), RDup.*sind(thetaDup).*sind(phiDup), RDup.*cosd(thetaDup)];
    thresholdBallRad = 0.1; %Ball radius around duplicate peak where peaks found are discarded

    for peakDupind = 1:length(thetaDup)
        peakMatrix = [rpeaks.*sind(thetapeaks).*cosd(phipeaks), rpeaks.*sind(thetapeaks).*sind(phipeaks), rpeaks.*cosd(thetapeaks)];
        newPeakInd = vecnorm(peakMatrix-dupMatrix(peakDupind,:),2,2)>thresholdBallRad;
        thetapeaks = thetapeaks(newPeakInd);
        phipeaks = phipeaks(newPeakInd);
        rpeaks = rpeaks(newPeakInd);
        Shat3D_peaks = Shat3D_peaks(newPeakInd);
    end

    %}
    [ShatM,pos] = max(Shat3D_peaks); theta_center = thetapeaks(pos); phi_center = phipeaks(pos); R_center = rpeaks(pos);
    ShatMaxdB = 20*log10(ShatM);

end
