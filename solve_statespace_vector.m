function q_n_dot = solve_statespace_vector(t_n, q_n, H)

% Create vector u, containing displacements at time n
u = q_n(1:length(q_n)/2);

% Create vector v, containing velocities at time n
v = q_n(length(q_n)/2+1:end);


%% Compute wind load scalar at time n
v_top   = v(end);   % [m/s] Velocity of discretized mass at beam's end
v_wind  = ...       % [m/s] Velocity of the wind at beam's end
    extract_wind_speed(t_n,H.v_wind.wind);
v_rel   = ...       % [m/s] Relative velocity at the Nacelle
    v_top-v_wind;
F_top   = ...       % [N]   Contribution of wind load
    0.01*0.5*H.rho_air*H.c_drag_top*H.A_top*v_rel*abs(v_rel);
% F_top = H.Fstatic;    % Static load


%% Compute wave load vector at time n

wdot_addt = sin(H.wave_omega*t_n+H.wave_eps);
wdot_waves = bsxfun(@times,wdot_addt,H.wdot_wav_vec);
wdot_waves_sum = sum(wdot_waves);

wdot_current = H.wdot_cur_vec;

wdot = bsxfun(@plus,wdot_waves_sum,wdot_current);
wdot = wdot.';

wdotdot_addt = cos(H.wave_omega*t_n+H.wave_eps); % omega product already performed in A3_waves
wdotdot_waves = bsxfun(@times,wdotdot_addt,H.wdotdot_wav_vec);
wdotdot_waves_sum = sum(wdotdot_waves);

wdotdot = wdotdot_waves_sum.';

% A = D for circular
% H.ddw     = 0.1*cos(t_n)*ones(H.N,1);
% H.dw      = 0.1*sin(t_n)*ones(H.N,1);
H.ddw       = wdotdot;
H.dw        = wdot;

% Checking sizes
sizewdotdot=size(H.ddw);
sizewdot=size(H.dw);

% F_top = H.Fstatic;
F_wind = zeros(H.N,1);
F_wind(end) = 2*F_top/H.dL;
F_Kbeam = H.Kbeam*u;
F_Csoil = H.Csoil*v;
F_Ksoil = H.Ksoil*u;
F_water = H.rho_wat*H.vol* (H.ddw + H.C_m*H.ddw + ...
    0.5*((H.D+0.5*H.t)/H.vol)*(H.C_drag*(H.dw-v).*abs(H.dw-v)));

a = H.Minv * (F_wind - F_Kbeam -F_Csoil -F_water - F_Ksoil);

% Combine velocity and acceleration vectors in derivative of state vector
q_n_dot = [v; a];

%% Runtime information
% clc
% disp('t_n equals')
% t_n
% disp('total time:')
% H.Tend
% disp('kinematics: wdot wdotdot u v')
% 
% % [wdot wdotdot u v]
% 
% disp('Forces: F_wind F_water F_Csoil F_Ksoil F_Kbeam')

% [F_wind F_water F_Csoil F_Ksoil F_Kbeam]



end