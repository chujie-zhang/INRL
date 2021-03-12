function [initial_pos] = initialPositioning(time,id,pseudo_ranges)
Define_Constants

%step b Compute the Cartesian ECEF positions of the satellites at time 0
total_sat_r_es_e = zeros(size(id,2),3);
for i=1:size(id,2)
    [sat_r_es_e,sat_v_es_e] = Satellite_position_and_velocity(time(1),id(i));
    total_sat_r_es_e(i,:) = sat_r_es_e;
end

%set initial data
r_ea = [0;0;0];

clockOffset = 0;
last_r_ea = [0;0;0]; 
thre = 0.10; 
error = inf;
while (error > thre)
    %step c Predict the ranges from the approximate user position 
    %to each satellite
    r_aj = zeros(size(id,2),1);
    for i=1:size(id,2)
        %implement recursion
        %initial range computation
        temp=eye(3,3)*total_sat_r_es_e(i,:)' - r_ea;
        r_a=sqrt( temp.' * temp);
        
        %Sagnac effect compensation matrix
        C_e = [1,omega_ie*r_a/c,0;
              -omega_ie*r_a/c,1,0;
               0,0,1];
        %recalcuate r_aj using C_e
        r_a = sqrt((C_e*total_sat_r_es_e(i,:)' - r_ea)' * ...
            (C_e*total_sat_r_es_e(i,:)' - r_ea));
        
        r_aj(i,:) = r_a;
    end
    
    %step d Compute the line-of-sight unit vector from the 
    %approximate user position to each satellite
    u_aj = zeros(3,size(id,2));
    for i=1:size(id,2)
        u_a = (C_e*total_sat_r_es_e(i,:)' - r_ea) / r_aj(i);
        u_aj(:,i) = u_a;
    end
    
    %step e Formulate the predicted state vector,
    %measurement innovation vector
    % and measurement matrix
    
    % predicted state vector
    x_minus = [r_ea;clockOffset];
    
    %measurement innovation vector
    dz = zeros(size(id,2),1);
    for i=1:size(id,2)
        dz(i,1) = pseudo_ranges(1,i) - r_aj(i,1) - clockOffset;
    end
    
    %measurement matrix
    H = zeros(size(id,2),4);
    for i=1:size(id,2)
        H(i,:) = [-u_aj(:,i)' 1];
    end
    
    %step f Compute the position and receiver clock offset
    %using unweighted least-squares
    x_plus = x_minus + pinv((H'*H))*H'*dz;
    %set current result
    r_ea = x_plus(1:3,1);
    clockOffset = x_plus(4,1);
    error = abs(norm(r_ea) - norm(last_r_ea));
    last_r_ea = r_ea;
end
%step g Convert this Cartesian ECEF position solution 
%to latitude, longitude and height
[latitude,longitude,height,~] = pv_ECEF_to_NED(r_ea,clockOffset);
initial_pos = [latitude;longitude;height];
end