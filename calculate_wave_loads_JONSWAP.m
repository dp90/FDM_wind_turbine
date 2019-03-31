function H = calculate_wave_loads_JONSWAP(H)

%% Constructing JONSWAP spectrum
H.H_sign= 7;        % [m]   Significant wave height
H.T_p   = 12;       % [s]   Peak period 

gamma   = 3.3;      % [-]   Peak of the wave
omega_p = ...       % [rad] Wave frequency
    2*pi / H.T_p; 
sigma   = ...       %
    @(omega) (omega <= omega_p)*0.07 + (omega > omega_p)*0.09;
A       = ...       % 
    @(omega) exp(-( (omega/omega_p - 1) ./ (sigma(omega)*sqrt(2) ) ).^2);
S_omega = ...       %
    @(omega) 320 * H.H_sign^2 / H.T_p^4 ./omega.^5 .* ...
    exp(-1950/(H.T_p).^4 ./omega.^4) .* gamma.^A(omega);

%% Computing kinematics

H.Nele = H.N;

omega   = linspace(0, 10, 1E3); % Generating omega
domega  = 1/(numel(omega)-1);
omega   = [0.0000001 omega(2:end)]; % testing

k_vec   = zeros(1, length(omega)); % Building vector for k

wdot_cur_vec = zeros(1, H.Nele);
wdot_wav_vec = zeros(length(omega),H.Nele);
wdotdot_wav_vec = zeros(length(omega),H.Nele);

for kk = 1:H.Nele % varying number of elements

    if (H.loc(kk)==2)
        % assign current velocity
        zloc    =(kk-0.5)*H.dL;
        zloc_sub=zloc-H.Lsoil;
        wdot_cur_vec(kk) =  ...
            H.w_c0*((zloc_sub+H.Lsub)/H.Lsub)^(H.alpha);        
    end

end

H.wdot_cur_vec = wdot_cur_vec;

figure()
plot(wdot_cur_vec)

for ii = 1:length(omega)
    wave_number = ...
        @(k) omega(ii).^2/H.g - k*tanh(k*H.Lsub);
    k_vec(ii)    = ...
        fzero(wave_number, omega(ii).^2/H.g); %results in k-vector, each wavenumber corresponding to the wave frequencies

    for jj = 1:H.Nele 
        
        if (H.loc(jj)==2)        
        zloc    =(jj-0.5)*H.dL;
        zloc_sub=zloc-H.Lsoil;

        wdot_wav_vec(ii,jj) = ...
            sqrt(2*S_omega(omega(ii))*domega)*omega(ii) ...
            *cosh(k_vec(ii)*(H.Lsub-zloc_sub)) ...
            /(sinh(k_vec(ii)*H.Lsub)); %velocities^ in x-direction
    
        wdotdot_wav_vec(ii,jj) = ...
            sqrt(2*S_omega(omega(ii))*domega)*(omega(ii)^2) ...
            *cosh(k_vec(ii)*(H.Lsub-zloc_sub)) ...
            /(sinh(k_vec(ii)*H.Lsub)); %accelerations^ in x-direction

        
        end
        
    end

end

H.wave_k    = k_vec.';
H.wave_eps  = 2*pi.*rand(length(omega),1);
H.wave_omega= omega.';

rowswitch = find(H.loc==2);
row1      = rowswitch(1);
rowend    = rowswitch(end);

wdot_wav_vec(:,row1:1:rowend) = fliplr(wdot_wav_vec(:,row1:1:rowend));
H.wdot_wav_vec=wdot_wav_vec;

wdotdot_wav_vec(:,row1:1:rowend) = fliplr(wdotdot_wav_vec(:,row1:1:rowend));
H.wdotdot_wav_vec=wdotdot_wav_vec;

%% Plotting wave kinematics

    figure()
    surf(H.dL*(1:H.Nele),omega,wdot_wav_vec(:,:))
    ylabel('Omega [rad/s]')
    xlabel('Length [-]')

%% Plotting JONSWAP spectrum
if H.test3 == 1
    figure()    
    plot(omega, S_omega(omega));
    xlabel('Omega [rad/s]')
    ylabel('S(omega) [m^2/s]')
end

end

