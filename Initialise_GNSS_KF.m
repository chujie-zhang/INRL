function [x_est,P_matrix] = Initialise_GNSS_KF(r_eb_e,v_eb_e,d_rho_c,dd_rho_c)
%Initialise_Integration_KF - Initializes the Integration KF state estimates
% and error covariance matrix for Workshop 2
%
% This function created 30/11/2016 by Paul Groves
%
% Outputs:
%   x_est                 Kalman filter estimates:

%   P_matrix              state estimation error covariance matrix

% Copyright 2016, Paul Groves
% License: BSD; see license.txt for details

% Begins

% Initialise state estimates
x_est = [r_eb_e;v_eb_e;d_rho_c;dd_rho_c];

% Initialise error covariance matrix
P_matrix =  zeros(8);
for i=1:3
    P_matrix(i,i) = 10^2; 
end
for i=4:6
    P_matrix(i,i) = 0.05^2; 
end
P_matrix(7,7) = 100000^2; 
P_matrix(8,8) = 200^2; 

% Ends