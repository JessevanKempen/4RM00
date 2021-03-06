function [] = Tcoeff()
% Purpose: To calculate the coefficients for the T equation.


% constants
global NPI NPJ LARGE Dt Cmu
% variables
global x x_u y y_v T Gamma SP Su F_u F_v relax_T T_old rho Istart Iend ...
    Jstart Jend b aE aW aN aS aP k Cp Tplus wood_I wood_J inlet_J Thot


Istart = 2;
Iend = NPI+1;
Jstart = 2;
Jend = NPJ+1;

convect();
calculateuplus();

for I = Istart:Iend
    i = I;
    for J = Jstart:Jend
        j = J;
        % Geometrical parameters
        % Areas of the cell faces
        AREAw = y_v(j+1) - y_v(j); % = A(i,J) See fig. 6.2 or fig. 6.5
        AREAe = AREAw;
        AREAs = x_u(i+1) - x_u(i); % = A(I,j)
        AREAn = AREAs;
        Acell = AREAw*AREAs;
        
        % The convective mass flux defined in eq. 5.8a
        % note:  F = rho*u but Fw = (rho*u)w = rho*u*AREAw per definition.    
        Fw = F_u(i,J)*AREAw;
        Fe = F_u(i+1,J)*AREAe;
        Fs = F_v(I,j)*AREAs;
        Fn = F_v(I,j+1)*AREAn;
        
        % The transport by diffusion defined in eq. 5.8b
        % note: D = mu/Dx but Dw = (mu/Dx)*AREAw per definition
        
        % The conductivity, Gamma, at the interface is calculated
        % with the use of a harmonic mean.       
        Dw = ((Gamma(I-1,J)*Gamma(I,J))/(Gamma(I-1,J)*(x(I) - x_u(i)) ...
            + Gamma(I,J)*(x_u(i) - x(I-1))))*AREAw;
        De = ((Gamma(I,J)*Gamma(I+1,J))/(Gamma(I,J)*(x(I+1) - x_u(i+1)) ...
            + Gamma(I+1,J)*(x_u(i+1) - x(I))))*AREAe;
        Ds = ((Gamma(I,J-1)*Gamma(I,J))/(Gamma(I,J-1)*(y(J) - y_v(j)) ...
            + Gamma(I,J)*(y_v(j) - y(J-1))))*AREAs;
        Dn = ((Gamma(I,J)*Gamma(I,J+1))/(Gamma(I,J)*(y(J+1) - y_v(j+1)) ...
            + Gamma(I,J+1)*(y_v(j+1) - y(J))))*AREAn;
            
        %% The source terms
        % Heat flux from wood block       
        if I < wood_I && J > inlet_J && J < wood_J+1
            SP(I,J) = -LARGE;
            Su(I,J) = LARGE*Thot;
        
        % Horizontal wall bottom, Vertical wall right, Vertical wall left,
        % wood block right
        elseif J == 2 || I == NPJ+1 || (I == 2 && J > wood_J) || (I == wood_I+1 && J > inlet_J && J < wood_J+1)          
            SP(I,j) = -rho(I,J) * Cmu^0.25 * k(I,J)^0.5 * Cp(I,J) / Tplus(I,J) *Acell;
            Su(I,J) = rho(I,J) * Cmu^0.25 * k(I,J)^0.5 * Cp(I,J) * T(I,J) / Tplus(I,J) *Acell;
            
        else    
            SP(I,J) = 0.;
            Su(I,J) = 0.;
        end
        
        %% The coefficients (hybrid differencing scheme)
        aN(I,J) = max([-Fn, Dn - Fn/2, 0.]);
        
        % Suppressing link at bottom wall
        if J==2                                 
            aS(I,J) = 0.;
        else
            aS(I,J) = max([ Fs, Ds + Fs/2, 0.]);
        end
        
        % Suppressing link at left wall
        if I == 2 && J > wood_J                   
            aW(I,J) = 0.;
        else
            aW(I,J) = max([ Fw, Dw + Fw/2, 0.]);
        end

        % Suppressing link at right wall
        if I ==NPJ+1                               
            aE(I,J) = 0.;
        else
            aE(I,J) = max([-Fe, De - Fe/2, 0.]);
        end
        
        aPold   = rho(I,J)*AREAe*AREAn/Dt;
         
        
        % eq. 8.31 without time dependent terms (see also eq. 5.14):
        aP(I,J) = aW(I,J) + aE(I,J) + aS(I,J) + aN(I,J) + Fe - Fw + Fn - Fs - SP(I,J) + aPold;
        
        % Setting the source term equal to b        
        b(I,J) = Su(I,J) + aPold*T_old(I,J);
        
        % Introducing relaxation by eq. 6.36 . and putting also the last
        % term on the right side into the source term b(i,J)
        
        aP(I,J) = aP(I,J)/relax_T;
        b(I,J)  = b(I,J) + (1.0 - relax_T)*aP(I,J)*T(I,J);
        
        % now the TDMA algorithm can be called to solve the equation.
        % This is done in the next step of the main program.
        
    end
end

end

