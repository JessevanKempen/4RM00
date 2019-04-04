function [] = calculateuplus()
% Purpose: To calculate uplus vplus and tw.

% constants
global NPI NPJ Cmu kappa ERough
% variables
global y rho k tw yplus yplus1 yplus2 uplus mu u 

viscosity();

for I = 1:NPI+1
    i = I;
    
%     if J == 2 || J == NPJ+1  (want ander model, snelheid is niet overal
%     parallel aan muur)

%         if yplus(I,J) < 11.63
%             tw(I,J)      = mu(I,J)*0.5*(u(i,J)+u(i+1,J))/(y(J)-y(J-1));
%             yplus(I,J)  = sqrt(rho(I,J)*abs(tw(I,J)))*(y(J)-y(J-1))/mu(I,J);
%             %yplus(I,2)   = yplus1(I,2);
%             uplus(I,2)   = yplus(I,2);
%     else
%         tw(I,J)      = rho(I,J)*Cmu^0.25*sqrt(k(I,J))*0.5*(u(i,J)+u(i+1,J))/uplus(I,J);
%         yplus(I,J)  = sqrt(rho(I,J)*abs(tw(I,J)))*(y(J)-y(J-1))/mu(I,J);
%         %yplus(I,2)   = yplus1(I,2);
%         uplus(I,J)   = log(ERough*yplus(I,J))/kappa;
%     end
%     
%     if yplus2(I,NPJ+1) < 11.63
%         tw(I,NPJ+1)      = mu(I,NPJ+1)*0.5*(u(i,NPJ+1)+u(i+1,NPJ+1))/(y(NPJ+2)-y(NPJ+1));
%         yplus1(I,NPJ+1)  = sqrt(rho(I,NPJ+1)*abs(tw(I,NPJ+1)))*(y(NPJ+2)-y(NPJ+1))/mu(I,NPJ+1);
%         yplus(I,NPJ+1)   = yplus1(I,NPJ+1);
%         uplus(I,NPJ+1)   = yplus(I,NPJ+1);
%     else
%         tw(I,NPJ+1)      = rho(I,NPJ+1)*Cmu^0.25*sqrt(k(I,NPJ+1))*0.5*(u(i,NPJ+1)+u(i+1,NPJ+1))/uplus(I,NPJ+1);
%         yplus1(I,NPJ+1)  = sqrt(rho(I,NPJ+1)*abs(tw(I,NPJ+1)))*(y(NPJ+2)-y(NPJ+1))/mu(I,NPJ+1);
%         yplus(I,NPJ+1)   = yplus1(I,NPJ+1);
%         uplus(I,NPJ+1)   = log(ERough*yplus(I,NPJ+1))/kappa;
%     end   
        
end
