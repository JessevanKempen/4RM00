function [] = calculateuplus()
% Purpose: To calculate uplus vplus tplus and tw.

% constants
global NPI NPJ Cmu kappa ERough sigmaturb
% variables
global x y rho k tw yplus yplus1 yplus2 yplus3 uplus T Tplus mu u v Pee wood_J inlet_J

viscosity();

for I = 1:NPI+1
    i = I;
    if yplus1(I,2) < 11.63
        tw(I,2)      = mu(I,2)*0.5*(u(i,2)+u(i+1,2))/(y(2)-y(1));
        yplus1(I,2)  = sqrt(rho(I,2)*abs(tw(I,2)))*(y(2)-y(1))/mu(I,2);
        yplus(I,2)   = yplus1(I,2);
        uplus(I,2)   = yplus(I,2);
        %Tplus(I,2)   = T(I,2);
    else
        tw(I,2)      = rho(I,2)*Cmu^0.25*sqrt(k(I,2))*0.5*(u(i,2)+u(i+1,2))/uplus(I,2);
        yplus1(I,2)  = sqrt(rho(I,2)*abs(tw(I,2)))*(y(2)-y(1))/mu(I,2);
        yplus(I,2)   = yplus1(I,2);
        uplus(I,2)   = log(ERough*yplus(I,2))/kappa;
        %Tplus(I,2)   = sigmaturb(I,2) * (uplus(I,2) + Pee(I,2));
    end
end

for I = 1:NPI+2
    i = I;
    if yplus1(I,2) < 11.63
       Tplus(I,2)   = T(I,2);
    else
       Tplus(I,2)   = sigmaturb(I,2) * (uplus(I,2) + Pee(I,2));
    end
end

for J = wood_J+1:NPJ+1
    j = J;
    if yplus2(2,J) < 11.63
        tw(2,J)      = mu(2,J)*0.5*(v(2,j)+v(2,j+1))/(x(2)-x(1));
        yplus2(2,J)  = sqrt(rho(2,J)*abs(tw(2,J)))*(x(2)-x(1))/mu(2,J);
        yplus(2,J)   = yplus2(2,J);
        uplus(2,J)   = yplus(2,J);
        %Tplus(2,J)   = T(2,J);
    else
        tw(2,J)      = rho(2,J)*Cmu^0.25*sqrt(k(2,J))*0.5*(v(2,j)+v(2,j+1))/uplus(2,J);
        yplus2(2,J)  = sqrt(rho(2,J)*abs(tw(2,J)))*(x(2)-x(1))/mu(2,J);
        yplus(2,J)   = yplus2(2,J);
        uplus(2,J)   = log(ERough*yplus(2,J))/kappa;
        %Tplus(2,J)  = sigmaturb(2,J) * (uplus(2,J) + Pee(2,J));
    end
end

for J = wood_J+1:NPJ+2
    j = J;
    if yplus2(2,J) < 11.63
       Tplus(2,J)   = T(2,J); 
    else   
       Tplus(2,J)  = sigmaturb(2,J) * (uplus(2,J) + Pee(2,J)); 
    end    
end

for J = 1:NPJ+1
    j = J;
    if yplus3(NPI+1,J) < 11.63          
        tw(NPI+1,J)      = mu(NPI+1,J)*0.5*(v(NPI+1,j)+v(NPI+1,j+1))/(x(NPI+2)-x(NPI+1));
        yplus3(NPI+1,J)  = sqrt(rho(NPI+1,J)*abs(tw(NPI+1,J)))*(x(NPI+2)-x(NPI+1))/mu(NPI+1,J);
        yplus(NPI+1,J)   = yplus3(NPI+1,J);
        uplus(NPI+1,J)   = yplus(NPI+1,J);
        %Tplus(NPI+1,J)   = T(NPI+1,J);
    else
        tw(NPI+1,J)      = rho(NPI+1,J)*Cmu^0.25*sqrt(k(NPI+1,J))*0.5*(v(NPI+1,j)+v(NPI+1,j+1))/uplus(NPI+1,J);
        yplus3(NPI+1,J)  = sqrt(rho(NPI+1,J)*abs(tw(NPI+1,J)))*(x(NPI+2)-x(NPI+1))/mu(NPI+1,J);
        yplus(NPI+1,J)   = yplus3(NPI+1,J);
        uplus(NPI+1,J)   = log(ERough*yplus(NPI+1,J))/kappa;
        %Tplus(NPI+1,J)   = sigmaturb(NPI+1,J) * (uplus(NPI+1,J) + Pee(NPI+1,J));
    end 
end 

for J = 1:NPJ+2
    j = J;
    if yplus3(NPI+1,J) < 11.63
       Tplus(NPI+1,J)   = T(NPI+1,J);
    else
       Tplus(NPI+1,J)   = sigmaturb(NPI+1,J) * (uplus(NPI+1,J) + Pee(NPI+1,J)); 
    end
end
        
end

