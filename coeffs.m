    function [a_u, a_w, b_u, b_w] = coeffs(sigu, sigw, dsigu2, dsigw2, uw, duw, tau, up, wp)
        A = 2*(sigu^2*sigw^2 - uw^2);
        b_u = (2*sigw^2/tau)^0.5;
        b_w = b_u;
        a_u = 1/A*b_u^2*(uw*wp - sigw^2*up) + 0.5*duw + 1/A*(sigw^2*dsigu2*wp*up - uw*dsigu2*wp^2 - uw*duw*up*wp + sigu^2*duw*wp^2);
        a_w = 1/A*b_w^2*(uw*up - sigu^2*wp) + 0.5*dsigw2 + 1/A*(sigw^2*duw*up*wp -uw*duw*wp^2 - uw*dsigw2*up*wp + sigu^2*dsigw2*wp^2);
    end

