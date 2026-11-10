function shallow_water_sediment_HLL
% =========================================================================
% Given System
%
%   d_t H + V d_x H + H d_x V              = 0                         (1a)
%   d_t V + g d_x H + V d_x V + g d_x Z    = g Sb - Cf V^2 / H         (1b)
%   d_t Z + a V^2 d_x V                    = vs/(1-p') (D - E)         (1c)
%   d_t C + V d_x C                        = vs/H      (E - D)         (1d)
%     H depth, V velocity, Z bed elevation, C suspended-sediment conc.
%     E = A1 |V|^5 / (1 + A2 |V|^5)            
%     D = ( r1 + r2 (Cf |V| / vs)^r3 ) C     
%
% CONSERVATIVE FORM (almost, -g H d_x Z is not conservative)
%   q = [ H , m , Z , s ]          m = H V ,  s = H C
%   F = [ m , m^2/H + g H^2/2 , (a/3)(m/H)^3 , m s / H ]
%   S = [ 0 ; -g H d_x Z + g H Sb - Cf V^2 ; vs/(1-p')(D-E) ; vs(E-D) ]
%
% BOUNDARY CONDITIONS  
%   inlet  x=0 : m = Q0 (4a), Z = Z0 (4b), C = c1 + c2 Q0^c3 
%   outlet x=L : m = kG max(Ueq-Z,0) sqrt(2 g max(H+Z-HL,0))   
%
% INITIAL CONDITIONS ARE ASSUMED AS EQUILIBRIUM BUT CAN BE CHANGED
% =========================================================================


%% ------------------------- parameters --------------------------------

g       = 9.81;       
Cf      = 2.0e-3;     
A       = 1.3e-7;     
c1      = 5.0e-3;      
c2      = 1.3e-11;     
c3      = 2.75;       
r1      = 1.0;        
r2      = 31.5;      
r3      = -1.46;       
nu      = 1.0e-6;      
p_prime = 0.30;        
Ds      = 2.5e-4;     
R       = 1.65;        
kG      = 0.60;       

% Adjustable parameters
Q0  = 2.0;             
Z0  = 0.0;             
HL  = 1.0;             


% Equilibrium

Heq = 2.0;             % also adjustable
Veq = Q0 / Heq;                                     
Zeq = Z0;                                            
Sb  = Cf * Veq^2 / (g * Heq);                       
Ueq = Q0 / (kG * sqrt(2*g*(Heq + Zeq - HL))) + Zeq; 

% Derived constants
Rep = Ds * sqrt(R*g*Ds) / nu;         
vs  = g * R * Ds^2 / (18*nu);          
a   = 39 * Cf^1.5 / ((1 - p_prime)*g*R);    
A1  = A * (sqrt(Cf) * Rep^0.6 / vs)^5;     
A2  = A1 / 0.3;                           

%% ------------------------- grid ------------------------------------
L   = 1000;            % domain length  
Nx  = 400;             % number of cells 
T   = 200;             % final time    
CFL = 0.9;             % CFL number     

dx = L / Nx;                        % cell width
xc = (0.5:1:Nx-0.5).' * dx;        % cell center

%% ------------------------- initial condition(adjustable) ---------------------
H = Heq * ones(Nx,1);
V = Veq * ones(Nx,1);
Z = Z0  * ones(Nx,1);


E0 = A1 * Veq^5 / (1 + A2 * Veq^5);
K0 = r1 + r2 * (Cf * Veq / vs)^r3;
C  = (E0 / K0) * ones(Nx,1);

C_in = c1 + c2 * Q0^c3;

% variables changed to conserved variables for scheme
q = [H, H.*V, Z, H.*C];

%% ------------------------- HLL -------------------------
t = 0;  step = 0;
while t < T 

    % ghost cells 
    qg = zeros(Nx+2,4);
    qg(2:end-1,:) = q;

    % inlet  (x=0)
    Hin      = q(1,1);
    qg(1,:)  = [Hin, Q0, Z0, Hin*C_in];

    % outlet (x=L)
    Hout     = q(end,1);
    Zout     = q(end,3);
    Cout     = q(end,4) / q(end,1);
    head     = max(Hout + Zout - HL, 0);      
    opening  = max(Ueq - Zout, 0);            
    qg(end,:) = [Hout, kG*opening*sqrt(2*g*head), Zout, Hout*Cout];

    % flux
    Hg = qg(:,1);  mg = qg(:,2);  Zg = qg(:,3);  sg = qg(:,4);
    Vg = mg ./ Hg;
    Fg = [ mg, mg.^2./Hg + 0.5*g*Hg.^2, (a/3)*Vg.^3, mg.*sg./Hg ];

    % HLL with Davis speed estimates
    hL = Hg(1:end-1);  hR = Hg(2:end);
    vL = Vg(1:end-1);  vR = Vg(2:end);
    cL = sqrt(g.*hL + abs(a)*g.*vL.^2);   
    cR = sqrt(g.*hR + abs(a)*g.*vR.^2);
    SL = min(vL - cL, vR - cR);
    SR = max(vL + cL, vR + cR);

    % CFL step
    dt = CFL * dx / max(max(abs(SL), abs(SR)));
    if t + dt > T, dt = T - t; end        

    % HLL numerical flux 
    FL = Fg(1:end-1,:);  FR = Fg(2:end,:);   % flux left / right
    qL = qg(1:end-1,:);  qR = qg(2:end,:);   % left / right conserved state
    % flux of q "star" state between SL and SR:
    Fhll = (SR.*FL - SL.*FR + (SL.*SR).*(qR - qL)) ./ (SR - SL);
    Fhat        = FL;                    % SL >= 0 
    mid         = (SL < 0) & (SR > 0);
    Fhat(mid,:) = Fhll(mid,:);           % SL < 0 < SR 
    rgt         = (SR <= 0);
    Fhat(rgt,:) = FR(rgt,:);             % SR <= 0 

    % conservative update
    qn = q - (dt/dx) * (Fhat(2:end,:) - Fhat(1:end-1,:));

    % original variables
    H = qn(:,1);
    V = qn(:,2) ./ qn(:,1);

    % momentum sources: -g H d_x Z (non-conservative, centered)
    dZdx    = (Zg(3:end) - Zg(1:end-2)) / (2*dx);
    qn(:,2) = qn(:,2) + dt * ( -g .* H .* dZdx + g .* H .* Sb );

    qn(:,2) = qn(:,2) ./ (1 + dt * Cf .* abs(V) ./ H);

    % sediment exchange with backward euler for physical correctness
    aV    = abs(V);
    Eero  = A1 .* aV.^5 ./ (1 + A2 .* aV.^5);         
    Kdep  = r1 + r2 .* (Cf .* aV ./ vs).^r3;      
    Cstar = qn(:,4) ./ H;                             
    rC    = dt * vs ./ H;                             
    Cnew  = (Cstar + rC .* Eero) ./ (1 + rC .* Kdep); 
    Dnew  = Kdep .* Cnew;                            
    qn(:,4) = H .* Cnew;                              
    qn(:,3) = qn(:,3) + dt * (vs / (1 - p_prime)) .* (Dnew - Eero);

    % step update
    q = qn;
    t = t + dt;  step = step + 1;

    % Negative depth error
    if min(q(:,1)) <= 0
        error('Negative depth at step %d, t = %.4g.', step, t);
    end
end

%% ------------------------- output ------------------------------
% output in original variables
H = q(:,1);
V = q(:,2) ./ q(:,1);
Z = q(:,3);
C = q(:,4) ./ q(:,1);
assignin('base','result', struct('x',xc,'H',H,'V',V,'Z',Z,'C',C,'t',t));
