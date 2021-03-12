function GNSS_results=GNSS
clear variables;
Define_Constants
%load data from csv file
pseudoRanges = csvread('Pseudo_ranges.csv');
pseudoRangeRates = csvread('Pseudo_range_rates.csv');

%store data into separate variables
time = pseudoRanges(2:end,1);
id = pseudoRanges(1,2:end);
pseudo_ranges = pseudoRanges(2:end,2:end);
pseudo_range_rates = pseudoRangeRates(2:end,2:end);

%step 1 because we dont know the initial position, so I use least-square to
%estimate the initial position until it dont converge
startingPosition = initialPositioning(time,id,pseudo_ranges);
%assume initial velocity is zero
startingVelocity = [0;0;0]; 

%step 2 implement the same method for all epochs
outlier_list = [0 0];
[positions,velocities,d_rho_c,dd_rho_c,outlier_solutions] = multipleEpochs(time,id,...
    pseudo_ranges,pseudo_range_rates,startingPosition,startingVelocity);


% find the outlier satellites list
for k=1:size(time,1)
    [index,sat_id_number] = max(abs(outlier_solutions(k,:)));
    if index > 0
       outlier_list = [outlier_list;k,sat_id_number]; 
    end
end
%disp(outlier_list);

%step 3 implement kalman filter 
GNSS_results = gnssKalmanFilter(time,id,pseudo_ranges,pseudo_range_rates, ...
    positions(:,1),velocities(:,1),d_rho_c,dd_rho_c,outlier_list);

%save result and write to a csv.file
gnss_result=zeros(size(time,1),5);
gnss_result(:,1)=time;
gnss_result(:,2:3)=GNSS_results(1:2,:)'*rad_to_deg;
gnss_result(:,4:5)=GNSS_results(4:5,:)';
writematrix(gnss_result,'GNSS.csv');

%draw
figure
plot(GNSS_results(1,:)'*rad_to_deg,GNSS_results(2,:)'*rad_to_deg);
xlabel('Latitude')
ylabel('Longitude')
title('Position-GNSS-only')

figure
plot(gnss_result(:,4))
xlabel('Time')
ylabel('velocity')
title('North velocity-GNSS-only')

figure
plot(gnss_result(:,5))
xlabel('Time')
ylabel('velocity')
title('East velocity-GNSS-only')
end


function [positions,velocities,clockOffset,clockOffset2,outlier_list] = ...
    multipleEpochs(time,sat_id,pseudo_ranges,pseudo_range_rates,...
    initial_positions,initial_velocirt)

Define_Constants
%step a. convert latitude, longititude and height to cartesian ECEF position
latitude = initial_positions(1,1); 
longitude = initial_positions(2,1);
height = initial_positions(3,1);

%change ned to ecef using the function given by workshop
[r_eb_e,v_eb_e] = pv_NED_to_ECEF(latitude,longitude,height,initial_velocirt);

%inital clock offset estimation
clockOffset = 0;
clockOffset2 = 0;

%skew symmetric matrix
omegaE = [0,-omega_ie,0;
         omega_ie,0,0;
         0,0,0];

%define variables for using later
positions = zeros(3,size(time,1));
velocities = zeros(3,size(time,1));
outlier_list = zeros(size(time,1),size(pseudo_ranges,2));

for i=1:size(time,1)
    
    %define variables for using later
    total_sat_r_es_e = zeros(3,size(sat_id,2));
    total_sat_v_es_e = zeros(3,size(sat_id,2));
    r_aj = zeros(1,size(sat_id,2));
    u_aj = zeros(3,size(sat_id,2));
    v_aj = zeros(1,size(sat_id,2));
    dz = zeros(size(sat_id,2),1);
    d_z = zeros(size(sat_id,2),1);
    H = zeros(size(sat_id,2),4);
    
    for j=1:size(sat_id,2)
        %get value for satellite
        %step b cartesian ecef positions of satellites at time 0
        [sat_r_es_e,sat_v_es_e] = ...
            Satellite_position_and_velocity(time(i),sat_id(j));
        total_sat_r_es_e(:,j) = sat_r_es_e';
        total_sat_v_es_e(:,j) = sat_v_es_e';
        
        %step c Predict range from the approximate user position
        temp=eye(3,3)*total_sat_r_es_e(:,j) - r_eb_e;
        r_a=sqrt(temp'*temp);
        
        %Sagnac effect compensation matrix
        C_e = [1,omega_ie*r_a/c,0;
            -omega_ie*r_a/c,1,0;
            0,0,1];
        
        %recalcuate r_aj
        [C_e,r_a] = Raj(r_eb_e,total_sat_r_es_e(:,j)');
        r_aj(:,j) = r_a;
        
        %step d compute line-of-sight unit vector for satellite
        u_a = (C_e*total_sat_r_es_e(:,j) - r_eb_e) / r_aj(:,j);
        u_aj(:,j) = u_a;
        
        %get velocity
        v_a = u_a'*(C_e*(total_sat_v_es_e(:,j) + ...
            omegaE*total_sat_r_es_e(:,j)) - (v_eb_e + omegaE*r_eb_e));
        v_aj(:,j) = v_a;
        
        %step e Formulate the predicted state vector, 
        %measurement innovation vector and
        %measurement matrix
        x_minus = [r_eb_e;clockOffset];
        %measurement innovation vector
        dz(j,1) = pseudo_ranges(i,j) - r_aj(1,j) - clockOffset;
        
        %measurement matrix
        H(j,:) = [-u_aj(:,j)' 1];
        
        %using the same method for velocity
        %predicted state vector
        x_minus_v = [v_eb_e;clockOffset2];
        %measurement innovation vector
        d_z(j,1) = pseudo_range_rates(i,j) - v_aj(1,j) - clockOffset2;
    end
    
    %check for outliers
    outlier_list(i,:) = findOutliers(H,dz);
    
    %f. Compute position and reciever clock offset using unweighted
    %least-squares
    x_new = x_minus + pinv(H'*H)*H'*dz;
    r_eb_e = x_new(1:3,1);
    clockOffset = x_new(4,1);
    
    %f. Compute velocity and reciever clock offset using unweighted
    %least-squares
    x_plus = x_minus_v + (H'*H)\H'*d_z;
    v_eb_e = x_plus(1:3,1);
    clockOffset2 = x_plus(4,1);
    
    %return in ecef format
    positions(:,i) = r_eb_e;
    velocities(:,i) = v_eb_e;
end
end

function gnss_solutions = gnssKalmanFilter(time,sat_id,pseudo_ranges,pseudo_range_rates, ...
    r_eb_e,v_eb_e,clockOffset,clockOffset2,outlier_list)

%workshop2 task 2b:  GNSS Kalman Filter Multiple Epochs with 8-state kalman
%filter

Define_Constants

%initalize matrices
solutions = zeros(6,size(time,1));
total_sat_r_es_e = zeros(3,size(sat_id,2));
total_sat_v_es_e = zeros(3,size(sat_id,2));
r_aj = zeros(1,size(sat_id,2));
r_aj_dot = zeros(1,size(sat_id,2));
d_z = zeros(2*size(sat_id,2),1);

%initalize kalman filter state vector
[x_est,P_matrix] = Initialise_GNSS_KF(r_eb_e,v_eb_e,clockOffset,clockOffset2);

%compute transition matrix
tau_s = 0.5;
phi = ...
    [eye(3,3),tau_s*eye(3,3),zeros(3,1),zeros(3,1);
    zeros(3,3),eye(3,3),zeros(3,1),zeros(3,1);
    zeros(1,3),zeros(1,3),1,tau_s;
    zeros(1,3),zeros(1,3),0,1];


%compute system noise covariance matrix
Sa = 5; 
Scphi = 0.01; 
Scf = 0.04;
Q = ...
    [1/3*Sa*tau_s^3*eye(3,3) 1/2*Sa*tau_s^2*eye(3,3) zeros(3,1) zeros(3,1);
    1/2*Sa*tau_s^2*eye(3,3) Sa*tau_s*eye(3,3) zeros(3,1) zeros(3,1);
    zeros(1,3) zeros(1,3) Scphi*tau_s + 1/3*Scf*tau_s^3 1/2*Scf*tau_s^2;
    zeros(1,3) zeros(1,3) 1/2*Scf*tau_s^2 Scf*tau_s];
tempPseudo_ranges=pseudo_ranges;
tempPseudo_ranges_rates=pseudo_range_rates;
tempId=sat_id;
for epoch=1:size(time,1)
    %use transition matrix to propogate state estimate
    x_k = phi*x_est;
    for i=2:size(outlier_list,1)
        if epoch==outlier_list(i,1)
            pseudo_ranges(:,outlier_list(i,2)) = [];
            pseudo_range_rates(:,outlier_list(i,2)) = [];
            sat_id(outlier_list(i,2)) = [];
        end
    end
    %propogate state covariance matrix
    P = phi*P_matrix*phi' + Q;
    
    %compute line of sight vectors
    clear u_a_all;
    clear d_z;
    for j=1:size(sat_id,2)
        
        r_eb_e = x_k(1:3,1);
        v_eb_e = x_k(4:6,1);
        %step b get value for satellite
        [sat_r_es_e,sat_v_es_e] = ...
            Satellite_position_and_velocity(time(epoch),sat_id(j));
        total_sat_r_es_e(:,j) = sat_r_es_e';
        total_sat_v_es_e(:,j) = sat_v_es_e';
        
        %step c. Predict range from the approximate user position
        r_a = sqrt((eye(3,3)*total_sat_r_es_e(:,j) - r_eb_e)' * ...
            (eye(3,3)*total_sat_r_es_e(:,j) - r_eb_e));
        %Sagnac effect compensation matrix
        C_e = [1,omega_ie*r_a/c,0;
            -omega_ie*r_a/c,1,0;
            0,0,1];
        %recalcuate r_aj
        [C_e,r_a] = Raj(r_eb_e,total_sat_r_es_e(:,j)');
        r_aj(:,j) = r_a;
        
        %compute line of sight vector
        u_a = (C_e*total_sat_r_es_e(:,j) - r_eb_e) / r_aj(:,j);
        u_a_all(:,j) = u_a;
        
        %calculate range rates for each satellite
        r_a_dot = u_a'*(C_e*(total_sat_v_es_e(:,j) + ...
            Omega_ie*total_sat_r_es_e(:,j)) - (v_eb_e + Omega_ie*r_eb_e));
        r_aj_dot(:,j) = r_a_dot;
        
        %formulate measurement innovation vector
        d_z(j,1) = pseudo_ranges(epoch,j) - r_aj(1,j) - x_k(7,1);
        d_z(j+size(sat_id,2),1) = pseudo_range_rates(epoch,j) ...
            - r_aj_dot(1,j) - x_k(8,1);
        
    end
    
    %compute measurement matrix
    R_k = zeros(2*size(sat_id,2),2*size(sat_id,2));
    for r = 1:size(sat_id,2)
        R_k(r,r) = 10^2;
        R_k(r+size(sat_id,2),r+size(sat_id,2)) = 0.05^2;
    end
    %H_k = zeros(2*size(sat_id,2),2*size(u_a_all,1)+2);
    clear H_k
    for k=1:size(u_a_all,2)
        H_k(k,:) = [-u_a_all(:,k).',zeros(1,3),1,0];
        H_k(k+size(sat_id,2),:) = [zeros(1,3),-u_a_all(:,k).',0,1];
        
    end

    %disp(size(R_k));
    %disp(size(H_k));
    %Compute Kalman Gain matrix
    K = P*H_k'/(H_k*P*H_k' + R_k);
    
    %update state estimates
    x_plus = x_k + K*d_z;
    P_plus = (eye(size(P,1)) - K*H_k)*P;
    
    %append solutions
    [L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(x_plus(1:3),x_plus(4:6));
    solutions(:,epoch) = [L_b;lambda_b;h_b;v_eb_n];
    
    %update variables
    x_est = x_plus;
    P_matrix = P_plus;
    pseudo_ranges=tempPseudo_ranges;
    pseudo_range_rates=tempPseudo_ranges_rates;
    sat_id=tempId;
end
gnss_solutions = solutions;
end

function outlier_index = findOutliers(H,d_z)
%define variables
listss = size(H,1);
%step a compute residuals vector
v = (H*pinv(H'*H)*H' - eye(listss))*d_z;

%step b compute residual covariance
C_v = (eye(listss) - H*inv(H'*H)*H')*5^2;

%step c compute normalized residuals and compare to threshold
outlier_index = zeros(1,listss);
for i=1:size(H,1)
    if norm(v(i)) > sqrt(C_v(i,i))*6
        outlier_index(i) = v(i);
    end
end
end

function [C,r_aj] = Raj(rea,rej)
Define_Constants;
r_aj = 0;
temp = inf;
%recursion to find r_aj when it converges
while r_aj~=temp
temp = r_aj;
C = [1,omega_ie*r_aj/c,0;
    -omega_ie*r_aj/c,1,0;
    0,0,1];
temp2=C*rej'-rea;
r_aj=sqrt(temp2'*temp2);
end
end
