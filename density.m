function [] = density()
% Purpose: To calculate the density in the fluid as a function of pressure and temperature

% constants
global NPI NPJ R
% variables
global rho p T pc P_ATM relax_rho p_abs

for I = 2:NPI+1
    for J = 2:NPJ+1
         rho(I,J) = (1-relax_rho)*rho(I,J) + relax_rho*(p_abs(I,J)) ./(R(I,J).*T(I,J));
    end
end

end