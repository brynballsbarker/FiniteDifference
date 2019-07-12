
function system = CreateSystem1(uinc)

global N m;
global r;
global r0 rmax;
global k;
global delta_r delta_theta;

delta_r = (rmax-r0)/(N-1);
delta_theta = 2*pi/m;
A = sparse(N*m, N*m);

% put it all together
A(1:m, 1:m) = eye(m);
start = m;

for j=2:N-1
    % create T
    d = 1/(r(j)^2*delta_theta^2);
    C = k^2 - 2/(delta_r^2) - 2/(r(j)^2 * delta_theta^2);
    T = sparse(m, m);
    T(1, 1:2) = [C d];
    T(1, m) = d;
    for h=2:m-1
        T(h, h-1:h+1) = [d C d];
    end
    T(m, 1) = d;
    T(m, m-1:m) = [d C];
    
    % create D+ and D-
    coefmin = 1/(delta_r^2) - 1/(2*r(j)*delta_r);
    coefplus = 1/(delta_r^2) + 1/(2*r(j)*delta_r);
    Dmin = coefmin*eye(m);
    Dplus = coefplus*eye(m);
    
    A(start+1:start+m, start-m+1:start) = Dmin;
    A(start+1:start+m, start+m+1:start+2*m) = Dplus;
    A(start+1:start+m, start+1:start+m) = T;
    start = start + m;
end

% create Ttilde
d = 1/(rmax^2*delta_theta^2);
c = 2*delta_r*(1i*k - 1/(2*rmax))*(1/(delta_r^2) + 1/(2*rmax*delta_r));
C = k^2 - 2/(delta_r^2) - 2/(rmax^2 * delta_theta^2) + c;
T = sparse(m, m);
T(1, 1:2) = [C d];
T(1, m) = d;
for h=2:m-1
    T(h, h-1:h+1) = [d C d];
end
T(m, 1) = d;
T(m, m-1:m) = [d C];
Ttilde = T;


% create Dtilde
d = 2/(delta_r^2);
Dtilde = d*eye(m);

% set last row of A
A((N-1)*m + 1:N*m, (N-2)*m + 1:N*m) = [Dtilde Ttilde]; 

% create F
F = sparse(N*m,1);
for j=1:m
    F(j,1) = -uinc(1,j);
end

% return full system
system = [A F];
    