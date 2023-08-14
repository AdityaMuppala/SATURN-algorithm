function Shat_3D = Shat3D_from_S11(S11n,Rrange,thetarange,phirange,k,fnpts,rho_n,phi_n,z_n,M,count)

Shat_3D = zeros(length(Rrange),length(thetarange),length(phirange));
[R,theta,phi] = ndgrid(Rrange,thetarange,phirange);
x = R.*sind(theta).*cosd(phi); y = R.*sind(theta).*sind(phi); z = R.*cosd(theta);

for n = 1:M
    %     n
    xxn = x-rho_n(n)*cosd(phi_n(n)); yyn = y-rho_n(n)*sind(phi_n(n)); zzn = z-z_n(n);
    mod_n = sqrt(xxn.^2+yyn.^2+zzn.^2);

    for f = 1:fnpts
        U = cosd(2*k(f)*(mod_n))-1i*sind(2*k(f)*(mod_n));
        Shat_3D = Shat_3D+U*S11n(f,n);
    end
end

end

