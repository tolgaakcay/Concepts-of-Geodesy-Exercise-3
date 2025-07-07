clc;
clear;
format long;
% Loading the dataset 'Orbit'
load Orbit
% Defining the orbit components and time vector
t = Orbit(:,1); x = Orbit(:,2); y = Orbit(:,3); z = Orbit(:,4);
% Defining the symbolic variables for derivation
syms ti x_v_poly y_v_poly z_v_poly
% Initial plot of location vs time
figure(1);
plot(t,x);
title('Initial Coordinates of the Orbit Components vs Time')
hold on
plot(t,y)
plot(t,z)
legend('X Coordinates', 'Y Coordinates', 'Z Coordinates');
xlabel('Time (Seconds)');
ylabel('Coordinate Values (m)');
grid on

% Creating the design matrix to find the coefficients
for i = 1:length(t)
   for j = 1:5
       A(i,j) = t(i)^(j-1);
   end 
end
% Least squares adjustment used to find the coefficinents matrix for each coordinate
A_x = inv(A' *A)*A'* x;
A_y = inv(A' *A)*A'* y;
A_z = inv(A' *A)*A'* z;
% Interpolation for finding the orbital parameter for each time value from 0s to 60s with 1s spacing
x_poly = []; y_poly = []; z_poly = []; t_full = [];
for i = 1:61
    x_poly(i) = A_x(1)*i^0 + A_x(2)*i^1 + A_x(3)*i^2 + A_x(4)*i^3 + A_x(5)*i^4;
    y_poly(i) = A_y(1)*i^0 + A_y(2)*i^1 + A_y(3)*i^2 + A_y(4)*i^3 + A_y(5)*i^4;
    z_poly(i) = A_z(1)*i^0 + A_z(2)*i^1 + A_z(3)*i^2 + A_z(4)*i^3 + A_z(5)*i^4;
    t_full(i) = i - 1;
end

% Plot of location vs time after the interpolation
figure(2);
plot(t_full,x_poly);
title('Coordinates of the Orbit Components vs Time after Operation')
hold on
plot(t_full,y_poly)
plot(t_full,z_poly)
legend('Xf Coordinates', 'Yf Coordinates', 'Zf Coordinates');
xlabel('Time (Seconds)');
ylabel('Coordinate Values (m)');
grid on

% Plot of initial values vs after operation values
figure(3);
plot(t_full,x_poly);
title('X Initial vs X Final Orbital Value Plot')
hold on
plot(t,x)
legend('Xf Coordinates', 'Xi Coordinates');
xlabel('Time (Seconds)');
ylabel('Coordinate Values (m)');
grid on 

figure(4);
plot(t_full,y_poly)
title('Y Initial vs Y Final Orbital Value Plot')
hold on
plot(t,y)
legend('Yf Coordinates', 'Yi Coordinates');
xlabel('Time (Seconds)');
ylabel('Coordinate Values (m)');
grid on 

figure(5);
plot(t_full,z_poly)
title('Z Initial vs Z Final Orbital Value Plot')
hold on
plot(t,z)
legend('Zf Coordinates', 'Zi Coordinates')
xlabel('Time (Seconds)')
ylabel('Coordinate Values(m)')
grid on 

% The Polynomials Defined with the Symbolic ti
for i = 1:61
    x_v_poly(i) = A_x(1)*ti^0 + A_x(2)*ti^1 + A_x(3)*ti^2 + A_x(4)*ti^3 + A_x(5)*ti^4;
    y_v_poly(i) = A_y(1)*ti^0 + A_y(2)*ti^1 + A_y(3)*ti^2 + A_y(4)*ti^3 + A_y(5)*ti^4;
    z_v_poly(i) = A_z(1)*ti^0 + A_z(2)*ti^1 + A_z(3)*ti^2 + A_z(4)*ti^3 + A_z(5)*ti^4;
end
% Taking the derivative of the polynomials to find velocity
x_velocity = []; y_velocity = []; z_velocity = []; 
for i = 1:61
       % Taking derivative in regards to ti
       x_vel(i) = diff(x_v_poly(i),ti);
       y_vel(i) = diff(y_v_poly(i),ti);
       z_vel(i) = diff(z_v_poly(i),ti);
end
for i = 1:61      
       % Substituting ti with the time value 
       x_velocity(i) = subs(x_vel(i),ti,t_full(i));
       y_velocity(i) = subs(y_vel(i),ti,t_full(i));
       z_velocity(i) = subs(z_vel(i),ti,t_full(i));
end
% Calculation of the Normal Velocity
normalv = [];
for i = 1:61
    % Calculation of the displacement through the distance formula, using
    % the orbital parameters
    normalv(i) = sqrt(x_velocity(i)^2 + y_velocity(i)^2 + z_velocity(i)^2)/1000;
end
% Plot of Normal Velocity
figure(6);
plot(t_full,normalv, 'LineWidth', 1.5)
title('Normal Velocity of the Satellite over Time')
legend('Normal Velocity (km/s)')
xlabel('Time (Seconds)');
ylabel('Displacement(km)');
grid on 
