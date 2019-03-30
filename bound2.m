function [] = bound()
% Purpose: Specify boundary conditions for a calculation

% constants
global NPI NPJ U_IN YMAX Cmu Ti
% variables
global y u v T m_in m_out y_v F_u k eps

% Fixed temperature/velocity at the upper and lower wall
% for J = 1:NPJ+2
%     
%     %u(2,J) = U_IN*1.5*(1.-(2.*(y(J)-YMAX/2.)/YMAX)^2); % inlet
% end


u(1:NPI+2, 1) = 0;              % bottom wall
u(NPI+2, 1:NPJ+2) = 0;          % right wall
u(1, ceil((NPJ+2)/3):NPJ+2) = 5;% left wall

v = u;                          % vertical boundary velocities

u(1, 2:ceil((NPJ+2)/3)-1) = U_IN*1 ; % inlet
v(1, 2:ceil((NPJ+2)/3)-1) = 0;       % inlet

% Temperature at the walls in Kelvin
T(1:NPI+2,1) = 100.; % bottom wall
T(1:NPI+2,NPJ+2) = 100.; % top wall

% begin: globcont();
% Purpose: Calculate mass in and out of the calculation domain to
%          correct for the continuity at outlet.
convect();

m_in = 0.;
m_out = 0.;

for J = 2:NPJ+1
    j = J;
    AREAw = y_v(j+1) - y_v(j); % See fig. 6.3
    m_in  = m_in  + F_u(2,J)*AREAw;
    m_out = m_out + F_u(NPI+1,J)*AREAw;
end
% end: globcont()

% Velocity and temperature gradient at outlet = zero:
% Correction factor m_in/m_out is used to satisfy global continuity
u(2:NPI+1,NPJ+2) = u(2:NPI+1,NPJ+1)*m_in/m_out;
v(2:NPI+1,NPJ+2) = v(2:NPI+1,NPJ+1);
k(2:NPI+1,NPJ+2) = k(2:NPI+1,NPJ+1);
eps(2:NPI+1,NPJ+2) = eps(2:NPI+1,NPJ+1);
k(1, 1:NPJ+2) = 2./3.*(U_IN*Ti)^2; % at inlet
eps(1, 1:NPJ+2) = Cmu^0.75 *k(1,1:NPJ+2).^1.5/(0.07*YMAX*0.5); % at inlet
T(2:NPI+1,NPJ+2) = T(2:NPI+1,NPJ+1);

end
