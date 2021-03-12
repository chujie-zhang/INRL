clear variables;
Define_Constants
%load data from another two methods, GNSS and dead reckoning.
dr_result=deadReckoning;
gnss_result=GNSS;
%store data in separate variables
time2=dr_result(:,1);
position=dr_result(:,2:3);
ins_dr_velocity=dr_result(:,4:5);

geodetic_position=gnss_result(1:3,:)';
referenced_velocity=gnss_result(4:6,:)';
%---------------------------------------------------------------------------------------%
%---------------------------------------task 2-------------------------------------------%
%---------------------------------------------------------------------------------------%
%define a 4 state kalman filter estimating north and east DR
%velocity error, DR latitude error and DR longitude error
%the state vector is thus
x=zeros(4,1);
newPosition = geodetic_position(1,1:2);
newVelocity = referenced_velocity(1,1:2);
sigma_v=0.1;
sigma_r=10;
[R_N,R_E] = Radii_of_curvature(geodetic_position(1,1));
%The state estimation error covariance matrix
%is therefore initialised at
P_plus=eye(4,4);
P_plus(1,1)=sigma_v^2;
P_plus(2,2)=sigma_v^2;
P_plus(3,3)=(sigma_r^2)/((geodetic_position(1,3)+R_N)^2);
P_plus(4,4)=(sigma_r^2)/((R_E+geodetic_position(1,3))...
^2*cos(geodetic_position(1,1))^2);

%ten steps of Kalman filter
for i=2:size(time2,1)
    [R_N,R_E] = Radii_of_curvature(geodetic_position(i,1));
    %step 1 Compute the transition matrix
    tau_s=0.5;
    phi=eye(4,4);
    phi(3,1)=tau_s/(R_N+geodetic_position(i-1,3));
    phi(4,2)=tau_s/((R_E+geodetic_position(i-1,3))*...
        cos(geodetic_position(i-1,1)));
    
    %step 2 Compute the system noise covariance matrix
    Q=zeros(4,4);
    S_DR=0.2;
    Q(1,1)=S_DR*tau_s;
    Q(1,3)=0.5*((S_DR*tau_s^2)/(R_N+geodetic_position(i-1,3)));
    Q(2,2)=S_DR*tau_s;
    Q(2,4)=0.5*((S_DR*tau_s^2)/((R_E+geodetic_position(i-1,3))*...
        cos(geodetic_position(i-1,1))));
    Q(3,1)=0.5*((S_DR*tau_s^2)/(R_N+geodetic_position(i-1,3)));
    Q(3,3)=(1/3)*((S_DR*tau_s^3)/(R_N+geodetic_position(i-1,3))^2);
    Q(4,2)=0.5*((S_DR*tau_s^2)/((R_E+geodetic_position(i-1,3))*...
        cos(geodetic_position(i-1,1))));
    Q(4,4)=(1/3)*((S_DR*tau_s^3)/((R_E+geodetic_position(i-1,3))^...
        2*cos(geodetic_position(i-1,1))^2));
    
    %step 3 Propagate the state estimates:
    x_minus=phi*x;
    
    %step 4 Propagate the error covariance matrix:
    P_minus=phi*P_plus*phi'+Q;
    
    %step 5 Compute the measurement matrix
    H=[0,0,-1,0;
        0,0,0,-1;
        -1,0,0,0;
        0,-1,0,0];
    
    %step 6 Compute the measurement noise covariance matrix
    positionError_std=5;
    velocityError_std=0.02;
    R=zeros(4,4);
    R(1,1)=positionError_std^2/(R_N+geodetic_position(i,3))^2;
    R(2,2)=positionError_std^2/((R_E+geodetic_position(i,3))^2....
    *cos(geodetic_position(i,1))^2);
    R(3,3)=velocityError_std^2;
    R(4,4)=velocityError_std^2;
    
    %step 7 Compute the Kalman gain matrix
    K=P_minus*H'* pinv(H*P_minus*H'+R);
    
    %step 8 Formulate the measurement innovation vector
    dz=[geodetic_position(i,1)-position(i,1)*deg_to_rad;
        geodetic_position(i,2)-position(i,2)*deg_to_rad;
        referenced_velocity(i,1)-ins_dr_velocity(i,1);
        referenced_velocity(i,2)-ins_dr_velocity(i,2)]-H*x_minus;

    %step 9 Update the state estimates
    x_plus=x_minus+K*dz;
    
    %step 10 Update the error covariance matrix
    P_plus=(eye(4,4)-K*H)*P_minus;

    %Use the Kalman filter estimates to correct the DR solution at each epoch
    newPosition(i,1)=position(i,1)*deg_to_rad-x_plus(3);
    newPosition(i,2)=position(i,2)*deg_to_rad-x_plus(4);
    newVelocity(i,1)=ins_dr_velocity(i,1)-x_plus(1);
    newVelocity(i,2)=ins_dr_velocity(i,2)-x_plus(2);
    
    %update for next epoch
    x=x_plus;
end
newPosition=[newPosition,geodetic_position(:,3)];
newVelocity=[newVelocity,referenced_velocity(:,3)];
%disp(newPosition*rad_to_deg);
%disp(roundn(newVelocity,-2));

%save result and write to a csv.file
integration_result=zeros(size(time2,1),6);
integration_result(:,1)=time2;
integration_result(:,2:3)=newPosition(:,1:2)*rad_to_deg;
integration_result(:,4:5)=newVelocity(:,1:2);
integration_result(:,6)=dr_result(:,6);
writematrix(integration_result,'integration_DR_and_GNSS.csv');

%draw 
figure
plot(integration_result(:,2),integration_result(:,3));
xlabel('Latitude')
ylabel('Longitude')
title('Position-integration-DR-and-GNSS')

figure
plot(integration_result(:,6))
xlabel('Time')
ylabel('Degrees')
title('Heading-integration-DR-and-GNSS')

figure
plot(integration_result(:,4))
xlabel('Time')
ylabel('velocity')
title('North velocity-integration-DR-and-GNSS')

figure
plot(integration_result(:,5))
xlabel('Time')
ylabel('velocity')
title('East velocity-integration-DR-and-GNSS')
