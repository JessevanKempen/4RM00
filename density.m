function [] = density()
% Purpose: To calculate the density in the fluid as a function of pressure and temperature

% constants
global NPI NPJ R
% variables
global rho p T pc

for I = 2:NPI+1
    for J = 2:NPJ+1
         rho(I,J) = rho(I,J) + p(I,J) ./(R(I,J).*T(I,J));
    end
end

end