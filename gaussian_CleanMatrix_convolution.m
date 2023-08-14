function GaussianCleanImage = gaussian_CleanMatrix_convolution(R_range,theta_range,phi_range,CleanMatrixTemp,fnpts,M)

CleanMatrixTemp( ~any(CleanMatrixTemp,2), : ) = [];
R_target_range = CleanMatrixTemp(:,2);
theta_target_range = CleanMatrixTemp(:,3);
phi_target_range = CleanMatrixTemp(:,4);
Shat_target_range = CleanMatrixTemp(:,1);

[R,theta,phi] = ndgrid(R_range,theta_range,phi_range);
GaussianCleanImage = zeros(size(R));

for Cleanind = 1:length(Shat_target_range)

    R_target = R_target_range(Cleanind);
    Scale_Factor = fnpts*M/((4*pi*R_target)^2);
    theta_target = theta_target_range(Cleanind);
    phi_target = phi_target_range(Cleanind);
    Shat_target = Shat_target_range(Cleanind);

    distX = (R.*sind(theta).*cosd(phi)-R_target.*sind(theta_target).*cosd(phi_target)).^2;
    distY = (R.*sind(theta).*sind(phi)-R_target.*sind(theta_target).*sind(phi_target)).^2;
    distZ = (R.*cosd(theta)-R_target.*cosd(theta_target)).^2;
    distsq = distX+distY+distZ;

    varInv  = 500; %This matches the range and angular resolution of the main beam from X-band array
    GaussianCleanImage = GaussianCleanImage+abs(Scale_Factor*Shat_target.*exp(-varInv*(abs(distsq))));
end
