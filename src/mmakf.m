%% Simulation of Multiple Model Adaptive Kalman Filter

clear;
close all;

% Generate reference trajectory: 
% - Start at (0,0) and move in positive X direction with 1 m/s.
N = 20;
t = 1:N; % time
ref = [1:N; zeros(1,20)];

% Generate noisy range measurements from three stations:
S = [ 0, 20,  10; 
     20, 20, -10];

z = zeros(3,20);
for i=1:20,
  z(1,i) = norm(ref(:,i)-S(:,1));
  z(2,i) = norm(ref(:,i)-S(:,2));
  z(3,i) = norm(ref(:,i)-S(:,3));
end

range_std = 1; % m^2
z = z + range_std * randn(3,20);
R = range_std^2 * eye(3);

% Three different velocity models (m/s): 
alpha = [0.8; 1; 1.2];
p_alpha = [1/3; 1/3; 1/3];

% system noise
sys_std = 0.1;
Q = sys_std^2 * eye(2);

% Initialize state and covariance of 3 Kalman filters
x_bank = zeros(6,1);
P_bank = 0.01 * [...
  eye(2);
  eye(2);
  eye(2)];
  
weight = zeros(3,N);

for i=1:N,
  
  % system noise
  w = sys_std * randn(2,1);
  
  for k=1:3,
    
    ks = 2*k-1;
    ke = 2*k;
    
    % prediction
    
    x_bank(ks:ke,1) = x_bank(ks:ke,1) + alpha(k) * [1;0] + w;
    P_bank(ks:ke,:) = P_bank(ks:ke,:) + Q;
    
    % Compute predicted range measurements and linearize for EKF.
    z_hat = zeros(3,1); H = zeros(3,2);
    
    for j=1:3,
      H(j,:) = (x_bank(ks:ke,1) - S(:,j))';
      z_hat(j) = norm(H(j,:));
      H(j,:) = -H(j,:) ./ z_hat(j);
    end
    
    % measurement update
    inno = z(:,i) - z_hat;
    inno_var = H *  P_bank(ks:ke,:) * H' + R;
    inv_inno_var = inv(inno_var);
    
    likelihood = exp(-0.5 * inno'* inv_inno_var * inno)/sqrt(2*pi * det(inno_var));
    p_alpha(k) = likelihood * p_alpha(k);
    
    K = P_bank(ks:ke,:) * H' * inv_inno_var;
    x_bank(ks:ke,1) = x_bank(ks:ke,1) + K * inno;
    P_bank(ks:ke,:) = (eye(2) - K * H) * P_bank(ks:ke,:);
    
  end
  
  % Normalize weight and store
  p_alpha = p_alpha ./ sum(p_alpha);
  weight(:,i) = p_alpha;
  

end  

% Plot Weight Adaptation History
% Notice that the weight of the 2nd KF approaches to 1 while others to zero.

figure
plot(t, weight(1,:),'r')
hold on;
plot(t, weight(2,:),'g')
plot(t, weight(3,:),'b')
hold off
xlabel('Time (s)')
title('MMAKF Weight Adaptation History')
legend('KF 1', 'KF 2', 'KF 3')
