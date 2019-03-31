clear all; 
clc; 
dbstop if error;
close('all');

H.test2=0;
H.test3=0;
Tend = 250;         % [s]

%% Model parameters
N = 20;
H.N = N;            % Number of discrete points in model
H.Nele = N-1;       % Number of elements in model
sigma_int_max_total = 0;
H.t_n = 0;

%% Turbine
H.D       = 2.5;
H.t       = 0.03;
H.A_st    = pi/4*(H.D + 0.5*H.t)^2 - pi/4*(H.D-0.5*H.t)^2;
H.I_st    = pi/64*(H.D + 0.5*H.t)^4 - pi/64*(H.D-0.5*H.t)^4;
H.Lsub    = 30;     % [m]   Submerged length
H.Lsoil   = 30;     % [m]   Penetration depth
H.Lbeam   = 96+H.Lsub;
H.L_tot = ...
    H.Lbeam + H.Lsoil;
H.dL = ...
    H.L_tot/(N-1);  % Distance between nodes
H.vol     = ...     % Volume of water displaced by 1 element of beam 
    (pi/4*(H.D + 0.5*H.t)^2*H.dL);
H.E       = 210E9;
H.rho_st  = 8050;
H.c_drag  = 0.47;   
H.c_addm  = 1;      % added mass from water 
H.c_m     = 1 + H.c_addm;
H.sigma_y = 250E6;
H.W = H.I_st/((H.D + H.t)/2);

%% Nacelle & Loading at top
H.m_top   = 3.5E5;
H.rho_air = 1.225;
H.c_drag_top = 0.05;
H.A_top   = pi*60^2;
H.v_wind  = load('wind_signal.mat');

%% Environmental
H.H_sign  = 7;      % [m]   Significant wave height
H.T_p     = 12;     % [s]   Peak Period
H.H_depth = 30;     % [m]   Water depth
H.rho_wat = 1025;   % [kg/m^3] Water density
H.g       = 9.81;   % [m/s^2] Gravitational constant
H.mu      = 1E-3;   % [Pa*s] Dynamic viscosity
H.c_soil = 50000;   % [N*s/m] soil damping parameter
H.k_soil = 50000;   % [N/m] soil stiffness parameter
H.L_tot = ...       % total beam length
    H.Lbeam + H.Lsoil;
H.w_c0    = 0.5;    
H.alpha   = 1/7;
H.ddw     = 0.1*ones(N,1);
H.dw      = 0.5*ones(N,1);

%% Creation of location vector: 1 = soil, 2 = water, 3 = air
loc = zeros(N,1);
for i = 1:N
    if (i-1)*H.L_tot/(N-1) <= H.Lsoil
        loc(i) = 1;
    elseif ( (i-1)*H.L_tot/(N-1) > H.Lsoil ) && ( (i-1)*H.L_tot/(N-1) <= H.Lsoil+H.Lsub)
        loc(i) = 2;
    else
        loc(i) = 3;
    end
end
H.loc = loc;


H = calculate_wave_loads_JONSWAP(H);



%% Stiffness matrix of beam
d = ...             % Specification of diagonals in stiffness matrix
    ones(N,1)*[1 -4 6 -4 1];
d(end-2,1) = 2;
d(end-1,2) = -4;
d(1,3) = 7;
d(end-1,3) = 5;
d(end,3) = 2;
d(end,4) = -2;

Kbeam = ...             % Construction of stiffness matrix
    spdiags((H.E*H.I_st/H.dL^4)*d,[-2 -1 0 1 2],N,N);
H.Kbeam = Kbeam;

%% Mass matrix (from beam and added mass from water)
% Diagonals from beam mass
m = ...             % mass of all nodes except node N = (rho A l^4)/EI
    ones(N,1)* (H.rho_st*H.A_st);
m(end) = ...        % mass of node N = (rho A l^4 + 2ml^3)/EI
    (H.rho_st*H.A_st + 2*H.m_top/H.dL);

% Diagonals from added mass
m_vec = ...         % Specify which entries are non-zero
    zeros(N,1) + (loc == 2);
rhoVc = ...         % Magnitude of non-zero values
    H.rho_wat * H.vol * H.c_m;

% Diagonal from beam mass + added mass
mtot = m + (m_vec*(rhoVc));
M = spdiags(mtot,0,N,N);  % Generation of sparse mass matrix
H.Minv = inv(M);

%% Stiffness and damping matrix of soil
d = ...             % Specification of diagonals in soil stiffness and damping matrices
    zeros(N,1) + (loc == 1);
Csoil = ...         % Formulation of damping matrix
    spdiags(H.c_soil*d,0,N,N);
Ksoil = ...         % Formulation of stiffness matrix
    spdiags(H.k_soil*d,0,N,N);
H.Ksoil = Ksoil;
H.Csoil = Csoil;

%% Damping matrix of added mass and of drag
d = ...             % Specification of diagonals in damping matrix
    zeros(N,1) + (loc == 2);
C_m = ...           % Generation of added mass damping matrix
    spdiags(H.c_m*d,0,N,N);
H.C_m = C_m;

d = ...             % Specification of diagonals in drag damping matrix
    zeros(N,1) + (loc == 2);
C_drag = ...        % Generation of drag damping matrix
    spdiags(H.c_drag*d,0,N,N);
H.C_drag = C_drag;

%% Static check
H.Fstatic = 10000;          % [N]
u_static = ...            % [m]
    H.Fstatic*H.L_tot^3/(3*H.E*H.I_st);

%% Calculation of initial conditions
% Force related to water
f_water = H.rho_wat*H.vol* (H.ddw + H.C_m*H.ddw + ...
    0.5*((H.D+0.5*H.t)/H.vol)*(H.C_drag*(H.dw).*abs(H.dw)));

% Force related to wind
v_wind  = ...       % [m/s] Velocity of the wind at beam's end
    H.v_wind.wind(2,1);
F_top   = ...       % [N]   Contribution of wind load
    0.01*0.5*H.rho_air*H.c_drag_top*H.A_top*v_wind*abs(v_wind);
f_wind = zeros(H.N,1);
f_wind(end) = 2*F_top/H.dL;
fext = f_wind + f_water;
Ktot = Kbeam + Ksoil;
u = inv(Ktot)*fext;


%% Running the model
plotfullbeam = 1;         % 1 for plot of entire beam, 0 for just the top


H.Tend = Tend;
q0 = [u' zeros(1,N)]';   % initial conditions
tvec=[H.v_wind.wind(1,1) H.v_wind.wind(1,Tend*100)];

if plotfullbeam == 0
    tic
    [T,Q] = ode45(@(t_n, q_n) solve_statespace_vector(t_n, q_n, H), tvec, q0);
    toc
    plot(T,Q(:,N))

elseif plotfullbeam == 1
    [T,W_beam] = ode45(@(t_n,q_n) solve_statespace_vector(t_n,q_n,H), tvec, q0);
    % Post
    W_beam_orig = W_beam;
    W_beam = [fliplr(W_beam(:,1:N-1)), fliplr(W_beam(:,N:end))];
    % Plot results
    % Open figure
    figure('units','normalized','outerposition',[0 0.1 0.45 0.85],...
        'PaperPositionMode','auto');
    % Plot
    title('Cantilever Beam Subjected To Dynamic Wind, Wave And Current Loading')
    x = linspace(0,H.L_tot,H.N);
    surface(T, x(2:end), rot90(W_beam(:,1:H.N-1)),'edgecolor','none')
    xlabel('Time [s]')
    ylabel('Beam Height [m]')
    zlabel('Displacement [m]')
end

%% Stress testing

deflections = W_beam_orig(:,1:H.N-1);
u_star = [zeros(numel(deflections(:,1)),1) deflections deflections(:,end)];

for qq = 1:numel(u_star(:,1))
    
    for zz = 2:(numel(u_star(1,:))-1)
        
        u_starrow = u_star(qq,:);
        
        udotdot = (u_starrow(zz-1) - 2*u_starrow(zz) + u_starrow(zz+1))/((H.dL)^2);        
        M_int = H.E*H.I_st*udotdot;
        sigma_int_max(qq,zz) = M_int/H.W/1E6;
        
    end
        
end

figure()
title('Bending Stresses Due To Environmental Loading')
   
surface(sigma_int_max(:,1:end-1),'edgecolor','none')
    ylabel('Time [s]')
    xlabel('Node number [-]')
    zlabel('Bending stress [N/m^2]')
hold on
limiter=H.sigma_y/1E6*ones(numel(sigma_int_max(:,1)),numel(sigma_int_max(1,:)));
surface(limiter,'FaceColor','r','EdgeColor','none');
    



