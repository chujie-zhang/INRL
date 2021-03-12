function headingResult = gyroMagnetometerIntegration(time,heading,gyro_heading)
Define_Constants

%apply two state kalman filter for Gyro-Magnetometer Integration

%initalize data
result = zeros(size(time,1),1);
h_minus = [0;0];
sigma_bias = 1;
sigma_heading = 4*deg_to_rad;
tau_s=0.5;
P_minus = [sigma_heading^2,0;0,sigma_bias^2];
%The heading error is the integral of the gyro bias
phi = [1,tau_s;
      0,1];

%Gyro random noise with power spectral density (PSD) 
S_rg =  1*10^-4;
%Gyro bias variation with PSD
S_bgd = 3*10^-6;

Q = [S_rg*tau_s+1/3*S_bgd*tau_s^3,1/2*S_bgd*tau_s^2;
              1/2*S_bgd*tau_s^2,S_bgd*tau_s];

for i=1:size(time,1)
    %step 1 use transition matrix to propogate state estimate
    x = phi*h_minus;
    
    %step 2 propagate error covariance matrix
    P = phi*P_minus*phi' + Q;
    
    %step 3 compute measurement matrix
    H_k = [-1 0];
    
    %step 4 Formulate measurement innovation vector
    d_z = heading(i) - gyro_heading(i) - H_k*x;
    
    %step 5 compute measurement noise covaraince matrix
    sigma_m = 4*deg_to_rad;
    R = diag([sigma_m^2]);
    
    %step 6 Compute Kalman Gain matrix
    K = P*H_k'/(H_k*P*H_k' + R);
    
    % step 7 update state estimates
    x_plus = x + K*d_z;
    P_plus = (eye(size(P,1)) - K*H_k)*P;
    
    %store results
    result(i,:) = (gyro_heading(i) - x_plus(1))';
    
    %update variables
    h_minus = x_plus;
    P_minus = P_plus;
end
headingResult = result';
end
