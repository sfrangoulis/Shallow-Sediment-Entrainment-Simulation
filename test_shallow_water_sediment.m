%% test_shallow_water_sediment.m
% =========================================================================
% Verification suite for  shallow_water_sediment_HLL.m
%
% Run:   matlab -batch test_shallow_water_sediment
% Exit code is non-zero if any test fails (CI-friendly).
%
% Most tests run against REFSOLVE (bottom of this file): a parametrised,
% headless twin of the solver.  It is numerically identical to
% shallow_water_sediment_HLL.m -- same HLL interface flux, same operator
% split, same semi-implicit friction and implicit erosion/deposition -- but
% takes a parameter struct and returns the state instead of plotting, so a
% single configuration can be swept.  Test 6 additionally runs the real
% shallow_water_sediment_HLL.m and checks the struct it leaves in the base
% workspace.
%
% What each test pins down
%   1  Uniform moving equilibrium is preserved to round-off
%      -> flux signs, source signs, friction/gravity balance, indexing.
%   2  Periodic conservation: total water mass, and total sediment
%      HC + (1-p')Z, are conserved to machine precision
%      -> telescoping fluxes + consistency of the coupled source step.
%   3  Still water over a flat bed stays at rest
%      -> hydrostatic pressure-flux balance, sources do not leak into V.
%   4  HLL flux consistency  Fhat(q,q) = F(q), and correct upwinding when
%      the fan is entirely supersonic.
%   5  Self-convergence under grid refinement is ~first order.
%   6  The real shallow_water_sediment_HLL.m runs, stays finite/positive,
%      honours its inlet-discharge and downstream-depth BCs.
%   7  Dam break: no NaN, depth stays positive, momentum conserved when
%      Cf = 0, discrete energy does not increase.
% =========================================================================

clear; clc; close all force;
RES = {};   % rows: {name, pass(logical), detail(string)}

%% ---------------- Test 1: uniform equilibrium is preserved --------------
try
    P = derive(base_params());  o = base_opts();
    Heq = o.Heq;  Veq = o.Q0/Heq;  Zeq = o.Z0;
    P.Sb  = P.Cf * Veq^2 / (P.g * Heq);                                % g Heq Sb = Cf Veq^2
    o.Ueq = o.Q0 / (P.kG * sqrt(2*P.g*(Heq + Zeq - o.HL))) + Zeq;      % eq (58)
    [H0,V0,C0] = equil_state(P,o);
    P.c1 = C0;  P.c2 = 0;            % inlet C_in = c1 + c2 Q0^c3 = C0 exactly
    o.T  = 50;
    out = refsolve(P,o);

    dH = max(abs(out.H - H0)) / H0;
    dV = max(abs(out.V - V0)) / abs(V0);
    dZ = max(abs(out.Z - o.Z0));
    dC = max(abs(out.C - C0)) / C0;

    RES(end+1,:) = {'equilibrium preserved: H', dH < 1e-9,  sprintf('rel err %.2e', dH)};
    RES(end+1,:) = {'equilibrium preserved: V', dV < 1e-9,  sprintf('rel err %.2e', dV)};
    RES(end+1,:) = {'equilibrium preserved: Z', dZ < 1e-10, sprintf('abs err %.2e', dZ)};
    RES(end+1,:) = {'equilibrium preserved: C', dC < 1e-9,  sprintf('rel err %.2e', dC)};
catch ME
    RES(end+1,:) = {'Test 1 (equilibrium)', false, ME.message};
end

%% ---------------- Test 2: periodic conservation ------------------------
try
    P = base_params();  P.Sb = 0;
    o = base_opts();  o.bc = 'periodic';  o.Nx = 64;  o.L = 200;  o.T = 30;
    x = ((0:o.Nx-1)+0.5).' * (o.L/o.Nx);
    o.U0 = [ 1.0  + 0.10*sin(2*pi*x/o.L), ...
             0.50 + 0.10*cos(4*pi*x/o.L), ...
             0.05*sin(2*pi*x/o.L), ...
             1e-3*(1 + 0.3*sin(2*pi*x/o.L)) ];
    out = refsolve(P,o);
    H = out.hist;                                  % [t sumH sumM sumS sumSed E]
    relH   = max(abs(H(:,2) - H(1,2))) / abs(H(1,2));
    relSed = max(abs(H(:,5) - H(1,5))) / abs(H(1,5));

    RES(end+1,:) = {'periodic: water mass conserved',            relH   < 1e-12, sprintf('rel err %.2e', relH)};
    RES(end+1,:) = {'periodic: total sediment HC+(1-p'')Z kept',  relSed < 1e-11, sprintf('rel err %.2e', relSed)};
    RES(end+1,:) = {'periodic: run stayed finite/positive',      ~out.crashed,   sprintf('min H %.3g', out.Hmin)};
catch ME
    RES(end+1,:) = {'Test 2 (conservation)', false, ME.message};
end

%% ---------------- Test 3: still water, flat bed, stays at rest ---------
try
    P = base_params();
    P.a = 0;  P.Sb = 0;  P.A = 0;  P.r2 = 0;  P.r3 = 1;   % kill erosion & V-dependent deposition
    o = base_opts();  o.bc = 'periodic';  o.Nx = 50;  o.L = 100;  o.T = 20;
    o.U0 = [ 2*ones(o.Nx,1), zeros(o.Nx,1), zeros(o.Nx,1), 1e-3*ones(o.Nx,1) ];
    out = refsolve(P,o);
    dV = max(abs(out.V));
    dH = max(abs(out.H - 2));

    RES(end+1,:) = {'lake at rest: V stays 0',        dV < 1e-12, sprintf('max|V| %.2e', dV)};
    RES(end+1,:) = {'lake at rest: H stays constant', dH < 1e-12, sprintf('max|H-2| %.2e', dH)};
catch ME
    RES(end+1,:) = {'Test 3 (lake at rest)', false, ME.message};
end

%% ---------------- Test 4: HLL flux consistency -----------------------
try
    P = derive(base_params());  rng('default');
    emax = 0;
    for k = 1:300
        q = [0.2 + 3*rand, (-2 + 4*rand), (-1 + 2*rand), 1e-3*rand];   % [H m Z s]
        F = phys_flux(q,P);
        [Fh,~] = hll_flux(q, q, F, F, P);
        emax = max(emax, max(abs(Fh - F) ./ (abs(F) + 1e-30)));
    end
    % entirely right-going fan (V > c on both sides)  ->  Fhat = F(qL)
    qL = [1.0, 9.0, 0.0, 1e-3];  qR = [1.3, 12.0, 0.1, 1.3e-3];
    [FhR,~] = hll_flux(qL,qR, phys_flux(qL,P), phys_flux(qR,P), P);
    okR = max(abs(FhR - phys_flux(qL,P))) < 1e-10;
    % entirely left-going fan  ->  Fhat = F(qR)
    qL = [1.0,-9.0,0.0,1e-3];  qR = [1.3,-12.0,0.1,1.3e-3];
    [FhL,~] = hll_flux(qL,qR, phys_flux(qL,P), phys_flux(qR,P), P);
    okL = max(abs(FhL - phys_flux(qR,P))) < 1e-10;

    RES(end+1,:) = {'HLL: Fhat(q,q) = F(q)',        emax < 1e-9, sprintf('max rel err %.2e', emax)};
    RES(end+1,:) = {'HLL: supersonic-right upwind', okR,         'Fhat = F(qL)'};
    RES(end+1,:) = {'HLL: supersonic-left upwind',  okL,         'Fhat = F(qR)'};
catch ME
    RES(end+1,:) = {'Test 4 (HLL flux)', false, ME.message};
end

%% ---------------- Test 5: self-convergence (smooth, periodic) --------
try
    P = base_params();  P.a = 0;  P.Sb = 0;
    o = base_opts();  o.bc = 'periodic';  o.L = 200;  o.T = 8;
    Ns = [40 80 160 320];  sols = cell(numel(Ns),1);
    for k = 1:numel(Ns)
        ok = o;  ok.Nx = Ns(k);
        xk = ((0:Ns(k)-1)+0.5).' * (o.L/Ns(k));
        ok.U0 = [ 1 + 0.15*sin(2*pi*xk/o.L), ...
                  0.5 + 0.10*sin(2*pi*xk/o.L + 1), ...
                  zeros(Ns(k),1), 1e-3*ones(Ns(k),1) ];
        rk = refsolve(P,ok);  sols{k} = rk.H;
    end
    e = zeros(numel(Ns)-1,1);
    for k = 1:numel(Ns)-1
        fine2coarse = mean(reshape(sols{k+1},2,[]),1).';   % restrict fine grid
        e(k) = mean(abs(sols{k} - fine2coarse));           % L1 self-error
    end
    ord = log2(e(1:end-1) ./ e(2:end));

    RES(end+1,:) = {'convergence: L1 self-error decreasing', all(diff(e) < 0), ...
                    sprintf('e = [%.2e %.2e %.2e]', e)};
    RES(end+1,:) = {'convergence: observed order ~1',        all(ord > 0.6) && mean(ord) < 1.8, ...
                    sprintf('orders = [%.2f %.2f]', ord)};
catch ME
    RES(end+1,:) = {'Test 5 (convergence)', false, ME.message};
end

%% ---------------- Test 6: run the real shallow_water_sediment_HLL ------
try
    if exist('shallow_water_sediment_HLL','file') ~= 2
        error('shallow_water_sediment_HLL.m not found on the path');
    end
    clear result;  close all force;
    shallow_water_sediment_HLL;               % assigns 'result' into base ws
    R = result;                              % this script IS the base ws
    close all force;

    fin  = all(isfinite([R.H;R.V;R.Z;R.C]));
    pos  = all(R.H > 0);
    tOK  = abs(R.t - 200) < 1e-6;            % T   in file-1 header
    nOK  = numel(R.x) == 400;                % Nx  in file-1 header
    qin  = R.H(1) * R.V(1);
    qOK  = abs(qin - 2) / 2 < 0.05;          % Q0  in file-1 header  (4a)
    hLOK = abs(R.H(end) - 2) < 0.02;         % equilibrium depth Heq = 2 m
    cOK  = all(R.C > -1e-10) && all(R.C < 1);

    RES(end+1,:) = {'file1: runs & output finite',     fin,        ''};
    RES(end+1,:) = {'file1: depth positive',           pos,        sprintf('min H %.3g', min(R.H))};
    RES(end+1,:) = {'file1: reached T with Nx cells',  tOK && nOK, sprintf('t=%.4g  n=%d', R.t, numel(R.x))};
    RES(end+1,:) = {'file1: inlet discharge = Q0',     qOK,        sprintf('VH|in = %.4g  (Q0 = 2)', qin)};
    RES(end+1,:) = {'file1: downstream H ~ Heq',       hLOK,       sprintf('H(end) = %.4g  (Heq = 2)', R.H(end))};
    RES(end+1,:) = {'file1: C within [0,1)',           cOK,        sprintf('C in [%.2e, %.2e]', min(R.C), max(R.C))};
catch ME
    RES(end+1,:) = {'Test 6 (real file 1)', false, ME.message};
end

%% ---------------- Test 7: dam break robustness / energy --------------
try
    P = base_params();
    P.a = 0;  P.Sb = 0;  P.Cf = 0;  P.A = 0;  P.r2 = 0;  P.r3 = 1;   % pure hydrodynamics
    o = base_opts();  o.bc = 'periodic';  o.Nx = 200;  o.L = 100;  o.T = 3;
    x  = ((0:o.Nx-1)+0.5).' * (o.L/o.Nx);
    H0 = 2*ones(o.Nx,1);  H0(x > o.L/2) = 1;
    o.U0 = [ H0, zeros(o.Nx,1), zeros(o.Nx,1), 1e-3*ones(o.Nx,1) ];
    out = refsolve(P,o);
    H = out.hist;
    momMax = max(abs(H(:,3)));                       % sum(m); starts at 0
    dE     = max(H(:,6)) - H(1,6);                   % energy rise (want <= 0)

    RES(end+1,:) = {'dam break: no NaN, no crash',           ~out.crashed && all(isfinite(out.H)), ...
                    sprintf('min H %.3g', out.Hmin)};
    RES(end+1,:) = {'dam break: depth stays positive',       out.Hmin > 0, sprintf('min H %.3g', out.Hmin)};
    RES(end+1,:) = {'dam break: momentum conserved (Cf=0)',  momMax < 1e-8, sprintf('max|sum m| %.2e', momMax)};
    RES(end+1,:) = {'dam break: discrete energy non-increasing', dE <= 1e-9*abs(H(1,6)) + 1e-12, ...
                    sprintf('dE = %.2e', dE)};
catch ME
    RES(end+1,:) = {'Test 7 (dam break)', false, ME.message};
end

%% ------------------------------ report ------------------------------
fprintf('\n=================== VERIFICATION RESULTS ===================\n');
nfail = 0;
for i = 1:size(RES,1)
    p = RES{i,2};
    p = ~isempty(p) && all(logical(p));
    nfail = nfail + ~p;
    tag = 'PASS';  if ~p, tag = 'FAIL'; end
    fprintf('  %s  %-44s  %s\n', tag, RES{i,1}, RES{i,3});
end
fprintf('-----------------------------------------------------------\n');
fprintf('  %d passed, %d failed  (of %d)\n\n', size(RES,1)-nfail, nfail, size(RES,1));
if nfail > 0
    error('test_shallow_water_sediment: %d test(s) failed', nfail);
end

% =========================================================================
%                          LOCAL FUNCTIONS
% =========================================================================
function P = base_params()
P.g       = 9.81;
P.Cf      = 2.0e-3;
P.A       = 1.3e-7;
P.c1      = 5.0e-3;
P.c2      = 1.3e-11;
P.c3      = 2.75;
P.r1      = 1.0;
P.r2      = 31.5;
P.r3      = -1.46;
P.nu      = 1.0e-6;
P.p_prime = 0.30;
P.Ds      = 2.5e-4;
P.R       = 1.65;
P.kG      = 0.60;
P.Sb      = 1.0e-3;
end
% -----------------------------------------------------------------------
function o = base_opts()
o.L   = 1000;
o.Nx  = 200;
o.T   = 50;
o.CFL = 0.9;
o.Q0  = 2.0;
o.Z0  = 0.0;
o.bc  = 'gate';
o.HL  = 1.0;
o.Heq = 2.0;
end
% -----------------------------------------------------------------------
function P = derive(P)
% coefficients from the given constants (paper eqs 47-51); an explicitly
% preset field is kept (lets a test force a = 0 etc.).
if ~isfield(P,'Rep'), P.Rep = P.Ds * sqrt(P.R*P.g*P.Ds) / P.nu;              end
if ~isfield(P,'vs'),  P.vs  = P.g * P.R * P.Ds^2 / (18*P.nu);                end   % Stokes
if ~isfield(P,'a'),   P.a   = 39 * P.Cf^1.5 / ((1 - P.p_prime)*P.g*P.R);     end
if ~isfield(P,'A1'),  P.A1  = P.A * (P.Cf / (P.vs * P.Rep^0.6))^5;           end
if ~isfield(P,'A2'),  P.A2  = P.A1 / 0.3;                                    end
end
% -----------------------------------------------------------------------
function [H0,V0,C0] = equil_state(P,o)
% uniform equilibrium: depth Heq carrying Q0, Sb chosen to match
P  = derive(P);
H0 = o.Heq;
V0 = o.Q0 / o.Heq;
E0 = P.A1 * V0^5 / (1 + P.A2 * V0^5);
K0 = P.r1 + P.r2 * (P.Cf * V0 / P.vs)^P.r3;
C0 = E0 / K0;
end
% -----------------------------------------------------------------------
function F = phys_flux(q,P)
H = q(:,1);  m = q(:,2);  s = q(:,4);  V = m ./ H;
F = [ m, ...
      m.^2 ./ H + 0.5 * P.g * H.^2, ...
      (P.a / 3) * V.^3, ...
      m .* s ./ H ];
end
% -----------------------------------------------------------------------
function [E,Kdep] = erosion_deposition(V,P)
aV   = abs(V);
V5   = aV.^5;
E    = P.A1 .* V5 ./ (1 + P.A2 .* V5);
Kdep = P.r1 + P.r2 .* (P.Cf .* aV ./ P.vs).^P.r3;
end
% -----------------------------------------------------------------------
function [Fhat,smax] = hll_flux(qL,qR,FL,FR,P)
% HLL interface flux (Davis speeds, a-coupling in the gravity speed).
% Identical formula to shallow_water_sediment_HLL.m.
hL = qL(:,1);  hR = qR(:,1);
vL = qL(:,2) ./ qL(:,1);  vR = qR(:,2) ./ qR(:,1);
cL = sqrt(P.g.*hL + abs(P.a)*P.g.*vL.^2);
cR = sqrt(P.g.*hR + abs(P.a)*P.g.*vR.^2);
SL = min(vL - cL, vR - cR);
SR = max(vL + cL, vR + cR);
SR = max(SR, SL + 1e-12);

Fhll = (SR.*FL - SL.*FR + (SL.*SR).*(qR - qL)) ./ (SR - SL);
Fhat = FL;
mid  = (SL < 0) & (SR > 0);   Fhat(mid,:) = Fhll(mid,:);
rgt  = (SR <= 0);             Fhat(rgt,:) = FR(rgt,:);
smax = max(max(abs(SL),abs(SR)));
end
% -----------------------------------------------------------------------
function [qg,Zg] = pad(q,o,P)
% ghost cells: 'periodic' wrap, or 'gate' = inlet (Q0,Z0,C_in) + sluice-gate
% outlet  m = kG (Ueq - Z) sqrt(2 g (H + Z - HL))  (eq 58 form)
Nx = size(q,1);
if strcmp(o.bc,'periodic')
    qg = [q(end,:); q; q(1,:)];
else
    qg = zeros(Nx+2,4);  qg(2:end-1,:) = q;
    Hin       = q(1,1);
    qg(1,1)   = Hin;
    qg(1,2)   = o.Q0;
    qg(1,3)   = o.Z0;
    qg(1,4)   = Hin * o.C_in;
    Hout      = q(end,1);
    Zout      = q(end,3);
    Cout      = q(end,4) / q(end,1);
    head      = max(Hout + Zout - o.HL, 0);
    opening   = max(o.Ueq - Zout, 0);
    qg(end,1) = Hout;
    qg(end,2) = P.kG * opening * sqrt(2 * P.g * head);
    qg(end,3) = Zout;
    qg(end,4) = Hout * Cout;
end
Zg = qg(:,3);
end
% -----------------------------------------------------------------------
function out = refsolve(P,o)
% Headless parametrised twin of shallow_water_sediment_HLL.m.
% On negative/non-finite depth it records out.crashed and stops (the real
% file calls error() there); everything else is byte-for-byte the same
% numerics.
P = derive(P);
Nx = o.Nx;  dx = o.L / Nx;
x  = ((0:Nx-1)+0.5).' * dx;

if isfield(o,'U0') && ~isempty(o.U0)
    H = o.U0(:,1);  V = o.U0(:,2);  Z = o.U0(:,3);  C = o.U0(:,4);
else
    [H0,V0,C0] = equil_state(P,o);
    H = H0*ones(Nx,1);  V = V0*ones(Nx,1);  Z = o.Z0*ones(Nx,1);  C = C0*ones(Nx,1);
end
if ~isfield(o,'C_in') || isempty(o.C_in)
    o.C_in = P.c1 + P.c2 * o.Q0^P.c3;
end

q = [H, H.*V, Z, H.*C];
rec = @(tt,w)[ tt, sum(w(:,1))*dx, sum(w(:,2))*dx, sum(w(:,4))*dx, ...
               sum(w(:,4) + (1-P.p_prime)*w(:,3))*dx, ...
               sum(0.5*w(:,2).^2./w(:,1) + 0.5*P.g*w(:,1).^2)*dx ];

t = 0;  step = 0;  crashed = false;
Hmin = min(q(:,1));
hist = rec(t,q);
maxsteps = 200000;

while t < o.T - 1e-12 && step < maxsteps
    [qg,Zg] = pad(q,o,P);
    Fg = phys_flux(qg,P);
    qL = qg(1:end-1,:);  qR = qg(2:end,:);
    FL = Fg(1:end-1,:);  FR = Fg(2:end,:);

    [Fhat,smax] = hll_flux(qL,qR,FL,FR,P);
    dt = o.CFL * dx / smax;
    if t + dt > o.T, dt = o.T - t; end

    qn = q - (dt/dx) * (Fhat(2:end,:) - Fhat(1:end-1,:));

    % --- source sub-step ---
    H = qn(:,1);  V = qn(:,2) ./ qn(:,1);
    dZdx    = (Zg(3:end) - Zg(1:end-2)) / (2*dx);
    qn(:,2) = qn(:,2) + dt * (-P.g .* H .* dZdx + P.g .* H .* P.Sb);
    qn(:,2) = qn(:,2) ./ (1 + dt * P.Cf .* abs(V) ./ H);

    [Eero,Kdep] = erosion_deposition(V,P);
    Cstar   = qn(:,4) ./ H;
    rC      = dt * P.vs ./ H;
    Cnew    = (Cstar + rC .* Eero) ./ (1 + rC .* Kdep);
    Dnew    = Kdep .* Cnew;
    qn(:,4) = H .* Cnew;
    qn(:,3) = qn(:,3) + dt * (P.vs / (1 - P.p_prime)) .* (Dnew - Eero);

    q = qn;  t = t + dt;  step = step + 1;
    Hmin = min(Hmin, min(q(:,1)));
    hist(end+1,:) = rec(t,q); %#ok<AGROW>

    if min(q(:,1)) <= 0 || any(~isfinite(q(:)))
        crashed = true;  break;
    end
end

out.x = x;
out.H = q(:,1);  out.V = q(:,2) ./ q(:,1);  out.Z = q(:,3);  out.C = q(:,4) ./ q(:,1);
out.t = t;  out.step = step;  out.hist = hist;  out.crashed = crashed;  out.Hmin = Hmin;
end
