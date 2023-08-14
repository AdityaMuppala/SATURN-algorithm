function [Shat_3D,phirot,thetarot,Rrot] = Shat3D_from_S11_sectioned(S11n,Rsweep,thetaSpan,thetaStep,phiStep,theta_center,phi_center,k,fnpts,rho_n,phi_n,z_n,M)

Ry = [cosd(theta_center), 0, sind(theta_center); 0, 1, 0; -sind(theta_center), 0, cosd(theta_center)];
Rz = [cosd(phi_center), -sind(phi_center), 0; sind(phi_center), cosd(phi_center), 0; 0, 0, 1;];

thetasweep = -thetaSpan:thetaStep:thetaSpan; phisweep = 0:phiStep:180;

[R1,theta1,phi1] = ndgrid(Rsweep,thetasweep,phisweep);
x1 = R1.*sind(theta1).*cosd(phi1); y1 = R1.*sind(theta1).*sind(phi1); z1 = R1.*cosd(theta1);
r1 = [x1(:).'; y1(:).'; z1(:).'];
r1new = Rz*Ry*r1;
x1rot = reshape(r1new(1,:),size(x1)); y1rot = reshape(r1new(2,:),size(y1)); z1rot = reshape(r1new(3,:),size(z1));
[phirot,elrot,Rrot] = cart2sph(x1rot,y1rot,z1rot);
thetarot = 90-rad2deg(elrot); phirot = rad2deg(phirot);

Shat_3D = zeros(size(x1rot));
for n = 1:M
    %     n
    xxn = x1rot-rho_n(n)*cosd(phi_n(n)); yyn = y1rot-rho_n(n)*sind(phi_n(n)); zzn = z1rot-z_n(n);
    mod_n = sqrt(xxn.^2+yyn.^2+zzn.^2);

    for f = 1:fnpts
        U = cosd(2*k(f)*(mod_n))-1i*sind(2*k(f)*(mod_n));
        Shat_3D = Shat_3D+U*S11n(f,n);
    end
end

% Displaying the maximum target location in the imaging region. Uncomment the remaining code to plot each step of the algorithm.
%{
for R_ind = 1:length(Rsweep)
    %
    figure(2)
    clf
    Shat_dB_norm = 20*log10(squeeze(abs(Shat_3D(R_ind,:,:))));%-max(max( 20*log10(squeeze(abs(Shat_3D(R_ind,:,:))))));
    surf(squeeze(x1rot(R_ind,:,:)),squeeze(y1rot(R_ind,:,:)),Shat_dB_norm) %%%%%%%%%
    xlabel('X (m)')
    ylabel('Y (m)')
    set(gca,'FontSize',20)
    shading interp
    view(180,270); axis tight; axis equal
    %     view(90,0)
    %     zlim([-100 0])
    colormap('jet');
    colorbar
    set(gca, 'clim', [-25 0]);

    w = waitforbuttonpress;

    %
end
%}

end

