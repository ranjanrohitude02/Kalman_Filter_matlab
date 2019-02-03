% KalmanFilterParabola
% name : Rohit Ranjan
% student number : 3029204
%% clean up
close all;
clear all;
clc;

a1 = -1*0.34;
kalman_filter(a1)

% estimation with aplha = 0.5

a2 = -1*0.05;
kalman_filter(a2)

function kalman_filter(alpha)
a = alpha;
b = 22;
c = 150;
x = linspace(1,40,400);
y = a*(x-b).^2+c;
figure();
plot(x,y,'b-');
hold on;

% state transition matrix
phi = [1 0 1 0; 0 1 0 1; 0 0 1 0; 0 0 0 1];
% measurement model matrix
H = [1 0 0 0; 0 1 0 0];
% system noise covariance matrix
Q = 10.^(-1*6)*[100 1 1 1; 1 100 1 1; 1 1 100 1; 1 1 1 100];
% measurement noise covariance matrix
R = 10.^(-1*3)*[500 1; 1 500];
% estimated process noise covariance matrix
C_bar = 10.^(-1*9)* eye(4);
% initial state estimation 
s_bar = transpose([0 0 0 0]);
% applying noise
L_trans = chol(R)'*randn(2,1);
L = transpose(L_trans);
for l = 1:length(x)
    chol_norm = chol(R)'*randn(2,1);
    x_norm(l) = x(l) + chol_norm(1);
    y_norm(l) = y(l) + chol_norm(2);
end   
plot(x_norm,y_norm,'r.');

% applying kalman filter
% Initialinsing vectors
x_final = zeros(1,length(x));
y_final = zeros(1,length(x));
for i = 1:length(x)
state_pred = phi* s_bar;
state_cov_pred = phi*C_bar*phi'+Q;
Gain = state_cov_pred*H'*pinv(H*state_cov_pred*H'+R);
m_k = [x_norm(i);y_norm(i)];
state_correct = state_pred + Gain*(m_k - H*phi*s_bar);
s_bar = state_correct;
I = eye(4);
state_cov_correct = (I - Gain*H)*state_cov_pred*(I-Gain*H)'+Gain*R*Gain';
C_bar = state_cov_correct
x_final(i) = state_correct(1);
y_final(i) = state_correct(2);
end
plot(x_final,y_final,'g.');

legend('true position','measured position','estimated position');
hold off;
