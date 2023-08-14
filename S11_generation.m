function S11n = S11_generation(M,rho_n,phi_n,z_n,Ntargets,rt_arr,thetat_arr,phit_arr,Shat_arr,fnpts,k)

%% S11 Nx1 vector generation for M Targets

S11n = zeros(fnpts,M); S11n_temp = zeros(fnpts,M);

for m = 1:Ntargets
    thetat = thetat_arr(m);
    phit = phit_arr(m);
    rt = rt_arr(m);
    for n = 1:M
        rn_vec = [rho_n(n)*cosd(phi_n(n)),rho_n(n)*sind(phi_n(n)),z_n(n)];
        rt_vec = [rt*sind(thetat)*cosd(phit), rt*sind(thetat)*sind(phit), rt*cosd(thetat)];
        mod_tn = norm(rn_vec-rt_vec);
        S11n_temp(:,n) = Shat_arr(m).*(1/(4*pi*mod_tn)^2)*(cosd(2*k*(mod_tn))+1i*sind(2*k*(mod_tn)));
    end
    S11n = S11n+S11n_temp;
end

end

