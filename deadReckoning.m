function dr_result=deadReckoning
clear variables
Define_Constants

%load data from csv.file
pseudoRanges = csvread('Pseudo_ranges.csv');
%save into separate variables
time = pseudoRanges(2:end,1);
id = pseudoRanges(1,2:end);
pseudo_ranges = pseudoRanges(2:end,2:end);
initialPosition = initialPositioning(time,id,pseudo_ranges);

%load data from csv.file and save into separate variables
file = csvread('Dead_reckoning.csv');
time = file(:,1);
left_front = file(:,2);
right_front = file(:,3);
left_back = file(:,4);
right_back = file(:,5);
gyro = file(:,6);
heading = file(:,7)*deg_to_rad;
forward_speed=((right_front+right_back)/2+(left_front+left_back)/2)/2;

%forward_speed = file(:,2);
%heading = file(:,3)*deg_to_rad;
h=initialPosition(3);
position = zeros(size(time,1),2);
position(1,:) = initialPosition(1:2);
%set initial data
average_velocity = zeros(size(time,1)-1,2);
ins_dr_velocity = zeros(size(time,1)-1,2);
ins_dr_velocity(1,1) = forward_speed(1)*cos(heading(1));
ins_dr_velocity(1,2) = forward_speed(1)*sin(heading(1));


%-------------------------------------------------------------------------%
%---------------------------------------task 1----------------------------%
%-------------------------------------------------------------------------%
for i=2:size(time,1)
    %compute the average velocity in north and east
    average_velocity(i,1) = 0.5*(cos(heading(i))+cos(heading(i-1)))...
        *forward_speed(i); %average_V_N
    average_velocity(i,2) = 0.5*(sin(heading(i))+sin(heading(i-1)))...
        *forward_speed(i);%averge_V_E
    
    %RN is the meridian radius of curvature and 
    %RE is the transverse radius of curvature
    [R_N,R_E] = Radii_of_curvature(position(i-1,1));
    
    %compute latitude and longitude from their counterparts
    position(i,1)=position(i-1,1)+(average_velocity(i,1)...
        *(time(i)-time(i-1)))/(R_N+h);
    position(i,2)=position(i-1,2)+(average_velocity(i,2)...
        *(time(i)-time(i-1)))/((R_E+h)*cos(position(i,1)));
    
    % compute the damped instantaneous DR velocity at each epoch
    ins_dr_velocity(i,1)=1.7*average_velocity(i,1)-0.7...
        *ins_dr_velocity(i-1,1);%V_N
    ins_dr_velocity(i,2)=1.7*average_velocity(i,2)-0.7...
        *ins_dr_velocity(i-1,2);%V_E
end
position=position*rad_to_deg;
ins_dr_velocity=roundn(ins_dr_velocity,-2);
%disp(position);
%disp(ins_dr_velocity);

%implement Gyro-Magnetometer Integration
gyro_heading = zeros(1,size(time,1));
gyro_heading(1) = heading(1);
for epoch=2:size(time,1)
   gyro_heading(epoch) = gyro_heading(epoch-1) + gyro(epoch)*0.5; 
end
heading_solutions = gyroMagnetometerIntegration(time,heading,gyro_heading);
heading=heading_solutions';


%save result and write to a csv.file
dr_result=zeros(size(time,1),6);
dr_result(:,1)=time;
dr_result(:,2:3)=position;
dr_result(:,4:5)=ins_dr_velocity;
dr_result(:,6)=heading*rad_to_deg;
writematrix(dr_result,'dead_reckon.csv');

%draw
figure
plot(position(:,1),position(:,2));
xlabel('Latitude')
ylabel('Longitude')
title('Position-DR-only')


figure
plot(time,ins_dr_velocity(:,1))
xlabel('Time')
ylabel('velocity')
title('North velocity-DR-only')

figure
plot(time,ins_dr_velocity(:,2))
xlabel('Time')
ylabel('velocity')
title('East velocity-DR-only')
end