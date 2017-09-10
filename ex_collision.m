lambda=1083e-9;     % Raman beam wavelength [m]
theta_raman=90;     % Raman beam angle [deg]
m_He=6.6465e-27;    % mass of helium atom [kg]
hbar=1.0546e-34;    % reduced Planck constant [m2kg/s]

k_photon=2*pi/lambda;       % k-vector for photon [m^-1]
K_c=hbar*k_photon*sind(theta_raman/2)         % BEC COM momentum (i.e. atom collision momentum) [kgm/s]
v_c=K_c/m_He        % helium collision velocity [m/s]

