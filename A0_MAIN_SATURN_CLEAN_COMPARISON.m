% This code implements Back-projection CLEAN and SATURN for an X-band SAR radar.

clc
clear
close all

fnpts = 201; %No. of frequency points
fstart = 8; fstop = 10; %Start and stop frequencies in GHz
f = linspace(fstart*1e9,fstop*1e9,fnpts); % Frequency range in Hz
frange = linspace(fstart,fstop,fnpts); % Frequency range in GHz
omega = 2*180*f;
c = physconst('lightspeed');
k = omega/c;

lambda_10GGHz = c./(1e10); % Center frequency
bw = 2e9; % Bandwidth in Hz
AFR = (fnpts-1)/(bw)*c/2; % Alias free range
RR = c/(2*bw); % Range resolution

%% S11 Nx1 vector generation for M Targets
% TXRX Array
M = 31; %number of antenna elements
a = 15*lambda_10GGHz; %Array radius
rho_n = a*linspace(1,1,M); % Radial coordinate
phi_n = linspace(360/M, 360,M); % Azimuthal coordinate
z_n = linspace(0,0,M); % Z-coordinate

% Target vector
Ntargs_iter = 5; %Total number of targets
iter = 1; %Iteration number

Ntargets = Ntargs_iter; 
rt_arr = 5+(rand(1,Ntargets)-0.5)*RR; % Radial coordinate of targets
thetat_arr = 7*ones(1,Ntargets); % Elevation coordinate of targets
phit_arr = linspace(0,360-360/Ntargets,Ntargets); % Azimuth coordinate of targets
Shat_arr_dB = -10*rand(1,Ntargets); % Target reflectivity in dB

Shat_arr = 10.^(Shat_arr_dB/20);
xarr = rt_arr.*sind(thetat_arr).*cosd(phit_arr); yarr = rt_arr.*sind(thetat_arr).*sind(phit_arr); zarr = rt_arr.*cosd(thetat_arr);

PlotSTAR(Shat_arr,rt_arr,thetat_arr,phit_arr,'go') %Plotting original target locations

S11_original = S11_generation(M,rho_n,phi_n,z_n,Ntargets,rt_arr,thetat_arr,phit_arr,Shat_arr,fnpts,k); % Generating S-parameter data at the antennas

%% Initializing the algorithm
%-- Finding the imaging window where the targets lie. The maximum target in the scene is
%found and the window is set to be around that target.
[ShatMaxdB,RBigTarget,thetaTarget,phiTarget] = LargestPeakSearchHoning_3D(S11_original, M, fstart, fstop, fnpts, RR, k, rho_n, phi_n, z_n);

Ntargets_total = Ntargets; %Number of targets to be decoupled in total

%Magnitude and location of each clean component
CleanMatrix = zeros(Ntargets_total,4); % Coloumn 1 is ShatClean, Coloumn 2 is RClean, Coloumn 3 is thetaClean and Coloumn 4 is phiClean
CleanMatrixOld = zeros(Ntargets_total,4); % Coloumn 1 is ShatClean, Coloumn 2 is RClean, Coloumn 3 is thetaClean and Coloumn 4 is phiClean
CleanMatrix(1,1) = 0.1; %Add fake value to prevent error later from CleanMatrix being all zeros

S11_predictionMatrix = zeros(fnpts,M,Ntargets_total);


%% CLEAN Algorithm
tic
for Ntargets = 1:Ntargets_total
    Ntargets
    S11_subtracted = S11_original-sum(S11_predictionMatrix,3); %Add all predicted S11's and subtract from measurement - This generates the residue S11

    %Step 1: Peak Search: This function scans the residual image and
    %finds the largest peak.
    [ShatMaxdB,RTarget,thetaTarget,phiTarget] = SinglepeakHoning_3D_withoutDuplicity(S11_subtracted, CleanMatrix, M, fstart, fstop, fnpts, RBigTarget, RR, k, rho_n, phi_n, z_n);

    %Step 2: Sweep and Extinguish: Although this step does not exist in the
    %original CLEAN algorithm, it is added here to extend the CLEAN for
    %near-field imaging. In this step, the location of the target inside
    %the range bin is found by sweeping it's location and subtracting it's
    %PSF. It also returns the predicted S11 of the target after it is
    %found.
    [ShatClean,RClean,S11_prediction,brightnessFunc] = CleanComponent_S11subtracted_Extraction(S11_subtracted, RTarget,thetaTarget,phiTarget,ShatMaxdB,fnpts, M, frange, fstart, fstop, k, rho_n, phi_n, z_n);
    S11_predictionMatrix(:,:,Ntargets) = S11_prediction;
    
    % Updating CLEAN components
    CleanMatrix(Ntargets,1) = ShatClean; CleanMatrix(Ntargets,2) = RClean; CleanMatrix(Ntargets,3) = thetaTarget; CleanMatrix(Ntargets,4) = phiTarget;

end

PlotSTAR(CleanMatrix(:,1),CleanMatrix(:,2),CleanMatrix(:,3),CleanMatrix(:,4),'m*')

CleanMatrixRand = CleanMatrix;
CleanMatrixRandOld = CleanMatrixOld;

time_CLEAN = toc;

%% Reconstructing De-convolved CLEAN image

S11_residual = S11_original-sum(S11_predictionMatrix,3); % Residual S11

%Imaging domain
R_range = 5;
theta_range = -15:0.1:15;
phi_range = 0:0.2:180;

% Convolving the clean components with a Gaussian
GaussianCleanImage = gaussian_CleanMatrix_convolution(R_range,theta_range,phi_range,CleanMatrixRand,fnpts,M);

% Residual image
Shat_3D_residual = Shat3D_from_S11(S11_residual,R_range,theta_range,phi_range,k,fnpts,rho_n,phi_n,z_n,M);

% Final CLEAN Image
Shat_CLEAN_abs = GaussianCleanImage+abs(Shat_3D_residual);

% Backprojection image of the same scene
Shat_3D_original = Shat3D_from_S11(S11_original,R_range,theta_range,phi_range,k,fnpts,rho_n,phi_n,z_n,M);
Shat_3D_original_abs = abs(Shat_3D_original);

disp('CLEAN done')
%}

%% SATURN Algorithm

SaturnMatrix2 = zeros(Ntargets_total,4); % Coloumn 1 is ShatClean, Coloumn 2 is RClean, Coloumn 3 is thetaClean and Coloumn 4 is phiClean
SaturnMatrixOld = zeros(Ntargets_total,4); % Coloumn 1 is ShatClean, Coloumn 2 is RClean, Coloumn 3 is thetaClean and Coloumn 4 is phiClean
SaturnMatrix2(1,1) = 0.1; %Add fake value to prevent error later from CleanMatrix being all zeros

S11_predictionMatrix = zeros(fnpts,M,Ntargets_total);
RSaturn = RBigTarget;

tic
for Ntargets = 1:Ntargets_total %Start by uncoupling two targets and then move up to uncouple one more at a time
    Ntargets
    S11_subtracted = S11_original-sum(S11_predictionMatrix,3); %Add all predicted S11's and subtract from measurement

    %Finding next peak from subtracted image in a 1m range around main peak
    [ShatMaxdB,RTarget,thetaTarget,phiTarget] = SinglepeakHoning_3D_withoutDuplicity(S11_subtracted, SaturnMatrix2, M, fstart, fstop, fnpts, RBigTarget, RR, k, rho_n, phi_n, z_n);

    %Finding target location inside range bin and returning Clean
    %component and S11_prediction
    [ShatSaturn,RSaturn,S11_prediction,brightnessFunc] = CleanComponent_S11subtracted_Extraction(S11_subtracted, RTarget,thetaTarget,phiTarget,ShatMaxdB,fnpts, M, frange, fstart, fstop, k, rho_n, phi_n, z_n);

    S11_predictionMatrix(:,:,Ntargets) = S11_prediction;
    SaturnMatrix2(Ntargets,1) = ShatSaturn; SaturnMatrix2(Ntargets,2) = RSaturn; SaturnMatrix2(Ntargets,3) = thetaTarget; SaturnMatrix2(Ntargets,4) = phiTarget;

    while_counter = 1;
    StoppingFactor = 20; % Increase this number to allow for a smaller convergence error i.e. more rounds of Recursive Decoupling

    while any(abs(SaturnMatrix2(:,1) - SaturnMatrixOld(:,1)) > SaturnMatrix2(:,1)/StoppingFactor)

        SaturnMatrixOld = SaturnMatrix2;
        if Ntargets == 1
            break
        end

        for targetNum = 1:Ntargets

            KnownPeak = SaturnMatrix2(targetNum,:);
            SaturnMatrix2(targetNum,:) = zeros(1,4); %Nulling the clean component that must be revised

            S11_predictionMatrix(:,:,targetNum) = zeros(fnpts,M);

            S11_subtracted = S11_original-sum(S11_predictionMatrix,3); %Add all predicted S11's and subtract from measurement

            %Finding next peak from subtracted image in a 1m range around
            %main peak
            [ShatMaxdB,RTarget,thetaTarget,phiTarget] = KnownpeakHoning_3D_withoutDuplicity(S11_subtracted, SaturnMatrix2, KnownPeak, M, fstart, fstop, fnpts, RR, k, rho_n, phi_n, z_n);

            %Finding target location inside range bin and returning Clean
            %component and S11_prediction
            [ShatSaturn,RSaturn,S11_prediction,brightnessFunc] = KnownComponent_S11subtracted_Extraction(S11_subtracted, RTarget,thetaTarget,phiTarget,ShatMaxdB,fnpts, M, frange, fstart, fstop, k, rho_n, phi_n, z_n);
            S11_predictionMatrix(:,:,targetNum) = S11_prediction;
            SaturnMatrix2(targetNum,1) = ShatSaturn; SaturnMatrix2(targetNum,2) = RSaturn; SaturnMatrix2(targetNum,3) = thetaTarget; SaturnMatrix2(targetNum,4) = phiTarget;

        end

        while_counter = while_counter+1;

        if while_counter > 10 % Increase this number to allow for more Recursive decoupling
            break
        end

    end

end
timeSATURN = toc;

PlotSTAR(SaturnMatrix2(:,1),SaturnMatrix2(:,2),SaturnMatrix2(:,3),SaturnMatrix2(:,4),'r*')

CleanMatrixRand = SaturnMatrix2;
CleanMatrixRandOld = SaturnMatrixOld;
disp('SATURN done')

%% Reconstructing De-convolved SATURN image

S11_residual = S11_original-sum(S11_predictionMatrix,3);

R_range = 5;
theta_range = -15:0.1:15;
phi_range = 0:0.2:180;

GaussianCleanImage = gaussian_CleanMatrix_convolution(R_range,theta_range,phi_range,CleanMatrixRand,fnpts,M);
Shat_3D_residual = Shat3D_from_S11(S11_residual,R_range,theta_range,phi_range,k,fnpts,rho_n,phi_n,z_n,M);
Shat_SATURN_abs = GaussianCleanImage+abs(Shat_3D_residual);

%% Peak_finding
%Finding the top N peaks in the image and comparing with actual target locations.

numPeaks = Ntargets_total; threshold = -35;

R_ind = 1;
r = R_range(R_ind);
[phi,theta] = meshgrid(phi_range,theta_range);
x = r.*sind(theta).*cosd(phi);
y = r.*sind(theta).*sind(phi);

Shat_2D_dB_original = 20*log10(Shat_3D_original_abs(R_ind,:,:));
Shat_2D_norm_original = abs(Shat_3D_original_abs(R_ind,:,:))./max(max(abs(Shat_3D_original_abs(R_ind,:,:))));
mask = imregionalmax(Shat_2D_norm_original);
[Shat2D_peaks_original,loc] = maxk(Shat_2D_dB_original(mask),numPeaks);
xtemp = x(mask); ytemp = y(mask);
xpeaks_original = xtemp(loc); ypeaks_original = ytemp(loc);
threshold_mask = Shat2D_peaks_original<threshold;
Shat2D_peaks_original(threshold_mask) = []; xpeaks_original(threshold_mask) = []; ypeaks_original(threshold_mask) = [];

Shat_2D_dB_CLEAN = 20*log10(Shat_CLEAN_abs(R_ind,:,:));
Shat_2D_norm_CLEAN = abs(Shat_CLEAN_abs(R_ind,:,:))./max(max(abs(Shat_CLEAN_abs(R_ind,:,:))));
mask = imregionalmax(Shat_2D_norm_CLEAN);
[Shat2D_peaks_CLEAN,loc] = maxk(Shat_2D_dB_CLEAN(mask),numPeaks);
xtemp = x(mask); ytemp = y(mask);
xpeaks_CLEAN = xtemp(loc); ypeaks_CLEAN = ytemp(loc);
threshold_mask = Shat2D_peaks_CLEAN<threshold;
Shat2D_peaks_CLEAN(threshold_mask) = []; xpeaks_CLEAN(threshold_mask) = []; ypeaks_CLEAN(threshold_mask) = [];

Shat_2D_dB_SATURN = 20*log10(Shat_SATURN_abs(R_ind,:,:));
Shat_2D_norm_SATURN = abs(Shat_SATURN_abs(R_ind,:,:))./max(max(abs(Shat_SATURN_abs(R_ind,:,:))));
mask = imregionalmax(Shat_2D_norm_SATURN);
[Shat2D_peaks_SATURN,loc] = maxk(Shat_2D_dB_SATURN(mask),numPeaks);
xtemp = x(mask); ytemp = y(mask);
xpeaks_SATURN = xtemp(loc); ypeaks_SATURN = ytemp(loc);
threshold_mask = Shat2D_peaks_SATURN<threshold;
Shat2D_peaks_SATURN(threshold_mask) = []; xpeaks_SATURN(threshold_mask) = []; ypeaks_SATURN(threshold_mask) = [];

%% Plotting

figure(2)
clf
surf(x,y,squeeze(Shat_2D_dB_original))
xlabel('X')
ylabel('Y')
set(gca,'FontSize',20)
shading interp
view(180,270); axis equal; axis tight
colormap('jet');
colorbar
set(gca, 'clim', [-35 0]);

figure(3)
clf
surf(x,y,squeeze(Shat_2D_dB_CLEAN))
hold on
xlabel('X')
ylabel('Y')
set(gca,'FontSize',20)
shading interp
view(180,270); axis equal; axis tight
colormap('jet');
colorbar
set(gca, 'clim', [-35 0]);

figure(4)
clf
surf(x,y,squeeze(Shat_2D_dB_SATURN))
hold on
xlabel('X')
ylabel('Y')
set(gca,'FontSize',20)
shading interp
view(180,270); axis equal; axis tight
colormap('jet');
colorbar
set(gca, 'clim', [-35 0]);
%}
