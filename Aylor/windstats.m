function [sigu, sigw, dsigu2, dsigw2, uw, duw, tau, Ubar] ...
    = windstats(z, L, ustar, z0, k, v_s, h, d, zw, ...
                Ubar_h, sigu_h, sigw_h, T_Lh, T_L0, beta, sigu_zw, sigw_zw)
     
        % at canopy boundary z = h
        
        zoverh = z/h;
        
        if z >= zw    %above roughness sublayer
            ex = (1-15*(z-d)/L)^0.25;                                       %Aylor 2001, just below eq. A2 
            psi = -2*log((1+ex)/2) - log((1+ex^2)/2) + 2*atan(ex) - pi/2;   %eq. A2 for L < 0
            sigu = sigu_zw;                                                 %eq. A8, calculated in param.m because it remains constant above zw for z >= zw, L < 0
            sigw = 1.25*ustar*(1-3*((z-d)/L))^(1/3);                        %eq. A9 for z >= zw, L < 0
            dsigu2 = 0;                                                     %derivative of sigu^2 wrt z 
            dsigw2 = -3.125*ustar^2*(1-3*(z-d)/L)^(-1/3)/L;                 %derivative of sigw^2 wrt z 
            uw = -ustar^2;                                                  %eq. A3 for z >= h
            duw = 0;                                                        %derivative of uw wrt z 
            T_L = (0.5*z/sigw)*(1-6*(z)/L)^0.25;                          %eq. A5 for z >= h
            tau = T_L/(1+(beta*v_s/sigw)^2)^0.5;                            %eq. 7 
            Ubar = ustar/k*(log((z-d)/z0) + psi);                           %eq. A1
            if tau < T_L0
                tau = T_L0;                                                 %from code "to keep high particles under control"
            end
            
            
            
            
        elseif zoverh < 1    %inside canopy
            gamma1 = 4;         
            gamma2 = 2.4;
            gamma3 = gamma2;    %gammas are attenuation factors from appendix part d Aylor 2001
            gamma4 = 3.9;
            Ubar = Ubar_h*exp(-gamma1*(1-zoverh)); %eq. A10, for 0 < z <= h
            sigu = sigu_h*exp(-gamma2*(1-zoverh)); %eq. A11, for 0 < z < h
            dsigu2 = 2*gamma2/h*sigu_h^2*exp(2*gamma2*(zoverh-1)); %derivative of sigu^2 wrt z 
            sigw = sigw_h*exp(-gamma3*(1-zoverh)); %eq. A12, for 0 < z < h
            dsigw2 = 2*gamma3/h*sigw_h^2*exp(2*gamma3*(zoverh-1)); %derivative of sigw^2 wrt z 
            uw = -ustar^2*exp(-gamma4*(1-zoverh)); %eq. A13, for 0 < z < h
            duw = -ustar^2*gamma4*exp(gamma4*(zoverh-1))/h; %derivative of uw wrt z 
            tau = T_Lh; %eq. A14 for 0.25h < z < h
            if zoverh < 0.25
                tau = T_Lh*(.1+3.6*zoverh); %eq. A15 for z < 0.25*h
            end
            
            
        else %in roughness sublayer
            ex = (1-15*(z-d)/L)^0.25; %Aylor 2001, just below eq. A2 
            psi = -2*log((1+ex)/2) - log((1+ex^2)/2) + 2*atan(ex) - pi/2; %eq. A2 for L < 0
 
            sigu = sigu_zw + 0.15*zw*sigu_zw*(zw-z)/(h-zw); %sigu decreases linearly by 15% from top to bottom of rough sublayer
            dsigu2 = 2*(sigu)*(-0.15*zw*sigu_zw/(h-zw)); %derivative of sigu^2 wrt z

            sigw = sigw_zw + 0.15*zw*sigw_zw*(zw-z)/(h-zw); %sigw decreases linearly by 15% from top to bottom of rough sublayer
            dsigw2 = 2*(sigw)*(-0.15*zw*sigw_zw/(h-zw)); %derivative of sigw^2 wrt z

            uw = -ustar^2;  %eq. A3 for z >= h
            duw = 0; %derivative of uw wrt z
            T_L = (0.5*z/sigw)*(1-6*(z)/L)^0.25; %eq. A5 for z >= h
            tau = T_L/(1+(beta*v_s/sigw)^2)^0.5; %eq. 7 
            Ubar = ustar/k*(log((z-d)/z0) + psi); %eq. A1
            if tau < T_L0
                tau = T_L0;
            end
            
        end
    end





