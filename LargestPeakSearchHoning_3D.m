function [ShatMaxdB, R_center,theta_center,phi_center] = LargestPeakSearchHoning_3D(S11n, M, fstart, fstop, fnpts, RR, k, rho_n, phi_n, z_n)

frange = linspace(fstart,fstop,fnpts);
bw = (fstop-fstart)*1e9;
%% Generating 3D Matrix with gating each bin

theta_span_array = [10 10 10 5];
theta_step_array = [1 1 0.5 0.1];
theta_center = 0;
phi_span_array = [180 180 60 30];
phi_step_array = [1 1 0.5 0.1];
phi_center = 90;
R_span_array = [1 0.2 0.2 0.2];
R_step_array = [RR RR/3 RR/3 RR/5];
R_center = 5;

for trial = 1:2
    thetarange = theta_center-theta_span_array(trial)/2:theta_step_array(trial):theta_center+theta_span_array(trial)/2;
    phirange = phi_center-phi_span_array(trial)/2:phi_step_array(trial):phi_center+phi_span_array(trial)/2;
    Rrange = R_center-R_span_array(trial)/2:R_step_array(trial):R_center+R_span_array(trial)/2;

    Shat_3D = Shat3D_from_S11(S11n,Rrange,thetarange,phirange,k,fnpts,rho_n,phi_n,z_n,M);
    % Displaying the maximum target location in the imaging region
    [ShatMaxdB,p] = max(20*log10(abs(Shat_3D(:))));
    sz = size(Shat_3D);

    [RMaxPos, thetaMaxPos, phiMaxPos] = ind2sub(sz,p);
    if trial > 1
        
        theta_center = thetarange(thetaMaxPos);
        phi_center = phirange(phiMaxPos);
    end

    R_center = Rrange(RMaxPos);
end

