function [] = epscoeff()
% Purpose: To calculate the epsilon.

% constants
global NPI NPJ Dt Cmu LARGE SMALL sigmaeps kappa C1eps C2eps
% variables
global x x_u y y_v SP Su F_u F_v mut rho Istart Iend ...
    Jstart Jend b aE aW aN aS aP k eps E2 eps_old uplus wood_I wood_J inlet_J


Istart = 2;
Iend = NPI+1;
Jstart = 2;
Jend = NPJ+1;

convect();
viscosity();

for I = Istart:Iend
    i = I;
    for J = Jstart:Jend
        j = J;
        
        % Geometrical parameters
        % Areas of the cell faces
        AREAw = y_v(j+1) - y_v(j); % = A(i,J) See fig. 6.3
        AREAe = AREAw;
        AREAs = x_u(i+1) - x_u(i); % = A(I,j)
        AREAn = AREAs;
        
        % eq. 6.9a-6.9d - the convective mass flux defined in eq. 5.8a
        % note:  F = rho*u but Fw = (rho*u)w = rho*u*AREAw per definition.
        Fw = F_u(i,J)*AREAw;
        Fe = F_u(i+1,J)*AREAe;
        Fs = F_v(I,j)*AREAs;
        Fn = F_v(I,j+1)*AREAn;
        
        % eq. 6.9e-6.9h - the transport by diffusion defined in eq. 5.8b
        % note: D = mu/Dx but Dw = (mu/Dx)*AREAw per definition
        % The conductivity, Gamma, at the interface is calculated with the use of a harmonic mean.
        Dw = mut(I-1,J)*mut(I,J)/sigmaeps/(mut(I-1,J)*(x(I) - x_u(i)) + ...
            mut(I,J)*(x_u(i)-x(I-1)))*AREAw;
        De = mut(I,J)*mut(I+1,J)/sigmaeps/(mut(I,J)*(x(I+1) - x_u(i+1)) + ...
            mut(I+1,J)*(x_u(i+1)-x(I)))*AREAe;
        Ds = mut(I,J-1)*mut(I,J)/sigmaeps/(mut(I,J-1)*(y(J) - y_v(j)) + ...
            mut(I,J)*(y_v(j)-y(J-1)))*AREAs;
        Dn = mut(I,J)*mut(I,J+1)/sigmaeps/(mut(I,J)*(y(J+1) - y_v(j+1)) + ...
            mut(I,J+1)*(y_v(j+1)-y(J)))*AREAn;
        
        % The source terms
        if J == 2 || I == NPJ+1 || (I == 2 && J > wood_J)
           SP(I,J) = -LARGE;
           Su(I,J) = Cmu^0.75*k(I,J)^1.5/(kappa*0.5*AREAw)*LARGE;
        
        
        % Hot wood source is source term for turbulence (eq. 9.23)
        elseif I < wood_I && J > inlet_J && J < wood_J
            SP(I,J) = -LARGE;
            Su(I,J) = Cmu^0.75*k(I,J)^1.5/(kappa*0.5*AREAw)*LARGE;

         else
            SP(I,J) = -C2eps*rho(I,J)*eps(I,J)/(k(I,J) + SMALL);
            Su(I,J) = C1eps*eps(I,J)/k(I,J)*2.*mut(I,J)*E2(I,J);
         end
             
        Su(I,J) =  Su(I,J)*AREAw*AREAs;
        SP(I,J) =  SP(I,J)*AREAw*AREAs;
        
        % The coefficients (hybrid differencing scheme)
        aN(I,J) = max([-Fn, Dn - Fn/2, 0.]);
        
        if J==2 %|| (I > 1+ceil(NPI/4) && I < 1+3*ceil(NPI/4) && J == 1+4*ceil(NPJ/5)) %Bottom or bottom of pan                                 
            aS(I,J) = 0.;
        else
            aS(I,J) = max([ Fs, Ds + Fs/2, 0.]);
        end
        
        if I == 2 && J > ceil(NPJ/3+1)              %Left
            aW(I,J) = 0.;
        else
            aW(I,J) = max([ Fw, Dw + Fw/2, 0.]);
        end

        if I ==NPJ+1                                %Right
            aE(I,J) = 0.;
        else
            aE(I,J) = max([-Fe, De - Fe/2, 0.]);
        end
        
        aPold   =  rho(I,J)*AREAe*AREAn/Dt;
        
        % eq. 8.31 with time dependent terms (see also eq. 5.14):
        aP(I,J) = aW(I,J) + aE(I,J) + aS(I,J) + aN(I,J) + Fe - Fw + Fn - Fs - SP(I,J) + aPold;
        
        % setting the source term equal to b
        b(I,J) = Su(I,J) + aPold*eps_old(I,J);
        
        % now the TDMA algorithm can be called to solve the equation.
        % This is done in the next step of the main program. */
        
    end
end

end
