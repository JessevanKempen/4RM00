function [] = bound()
% Purpose: Specify boundary conditions for a calculation

% constants
global NPI NPJ U_IN XMAX YMAX Cmu Ti
% variables
global y x u v T m_in m_out y_v x_u F_u F_v k eps Tinlet Thot inlet_J wood_J wood_I P_ATM p

% Fixed temperature/velocity at the upper and lower wall
% for J = 1:NPJ+2
%           u(2,1:NPJ+2) = U_IN; % inlet (zijkanten extra op 0 zetten!)
%     u(2,J) = U_IN*1.5*(1.-(2.*(y(J)-YMAX/2.)/YMAX)^2); % inlet
% end

%% left wall conditions
u(2,2:inlet_J) = U_IN;                 % inlet uniform x velocity
u(2, inlet_J+1:NPJ+2) = 0;             % left boundary has x velocity 0
v(1,2:NPJ+2) = 0;                      % left boundary has y velocity 0

%% bottom wall conditions
u(2:NPI+2, 1) = 0;                     % bottom wall has x velocity 0
v(1:NPI+2, 2) = 0;                     % bottom wall has y velocity 0

%% right wall conditions
u(NPI+2, 1:NPJ+2) = 0;                 % right wall has x velocity 0
v(NPI+2, 2:NPJ+2) = 0;                 % right wall has y velocity 0

%% Temperature at the walls in Kelvin
T(1:NPI+2,1) = T(1:NPI+2,2);                             % gradient 0 over bottom wall
T(NPI+2, 1:NPJ+2) = T(NPI+1, 1:NPJ+2);                   % gradient 0 over right wall
T(1, wood_J+1:NPJ+2) = T(2, wood_J+1:NPJ+2);             % gradient 0 over left wall

T(1, 2:inlet_J) = Tinlet;                                % Inlet T Left wall
T(1, inlet_J+1:wood_J) = Thot;                           % Wood continues left on boundary

%% Pressure conditions
%p(1, 2:inlet_J) = P_ATM;                              % Inlet P is ambient pressure
%p(2, 2) = P_ATM;                                      % Inlet P is ambient pressure
%p(2:NPI+1,NPJ+2) = P_ATM;



% begin: globcont();
% Purpose: Calculate mass in and out of the calculation domain to
%          correct for the continuity at outlet.
convect();

m_in = 0.;
m_out = 0.;

% for J = 2:13
%     j = J;
%     AREAw = y_v(j+1) - y_v(j); % See fig. 6.3
%     m_in  = m_in  + F_u(2,J)*AREAw;
% end
% 
% for I = 4:NPI+1
%     i = I;
%     AREAs = x_u(i+1) - x_u(i); % See fig. 6.3
%     m_out = m_out + F_u(I,NPJ+1)*AREAs;   
% end

for J = 2:NPJ+1
    j = J;
    AREAw = y_v(j+1) - y_v(j); % See fig. 6.3
    m_in  = m_in  + F_u(2,J)*AREAw;
end

for I = 2:NPI+1
    i = I;
    AREAs = x_u(i+1) - x_u(i); % See fig. 6.3
    m_out = m_out + F_v(I,NPJ+1)*AREAs;
end

% end: globcont()

% Velocity and temperature gradient at outlet = zero:
% Correction factor m_in/m_out is used to satisfy global continuity
u(3:NPI+1,NPJ+2) = u(3:NPI+1,NPJ+1)*m_in/m_out;
v(2:NPI+1,NPJ+2) = v(2:NPI+1,NPJ+1)*m_in/m_out;

k(2:NPI+1,NPJ+2) = k(2:NPI+1,NPJ+1);
eps(2:NPI+1,NPJ+2) = eps(2:NPI+1,NPJ+1);
p(2:NPI+1,NPJ+2) = p(2:NPI+1,NPJ+1);                     % pressure gradient 0 at outlet

T(2:NPI+2,NPJ+2) = T(2:NPI+2,NPJ+1);

k(1,1:NPJ+2)     = 1.5*(U_IN*Ti)^2; % at inlet
eps(1,1:NPJ+2)   = Cmu^0.75 *k(1,1:NPJ+2).^1.5/(0.07*YMAX*0.5); % at inlet



end
