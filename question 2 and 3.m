clc;
clear;
close all;
A = importdata("pHdata.dat");
time_steps = A(:,1);
u1 = A(:,2);
u2 = A(:,3);
y = A(:,4);
u1 = (u1 - min(u1))/(max(u1) - min(u1));
u2 = (u2 - min(u2))/(max(u2) - min(u2));
y= (y - min(y))/(max(y) - min(y));

U = [ones(size(u1)) u1.^2 u2.^2];
figure(1);
scatter(u1,y);
xlabel("Acid solution flow in liters");
ylabel('pH of the solution in the tank');
figure(2);
scatter(u2,y);
xlabel('Base solution flow in liters');
ylabel('pH of the solution in the tank');
figure(3);
scatter3(u1,u2,y,'filled');
xlabel("Acid solution flow in liters");
ylabel('Base solution flow in liters');
zlabel('pH of the solution in the tank');
%% least square
theta = lsqr(U,y);
intercept = theta(1);
slope1 = theta(2);
slope2 = theta(3);
y_fit = intercept + slope1 * u1 + slope2 * u2;
figure;
scatter3(u1,u2,y,'filled');
hold on;
scatter3(u1,u2,y_fit,'yellow','filled');
xlabel("Acid solution flow in liters");
ylabel('Base solution flow in liters');
zlabel('pH of the solution in the tank');
legend('actual data','fitted line');
title("fitted line using least sqare method");
grid on;
hold off;
%Error
error = y - y_fit;
figure;
plot(error , 'k', 'LineWidth',2 );
xlabel("data point index");
ylabel("error");
title('error plot');
grid on;
%% forgetting factor
landa = 0.99;
theta = zeros(size(U,2),1);
Q = eye(size(U,2))/landa ;
for i = 1:length(y)
    u_i = U(i,:)';
    new_y = u_i'*theta;
    e = y(i) - new_y;
    k = Q * u_i/(landa + u_i' * Q * u_i);
    theta = theta + k*e;
    Q = (Q - k * u_i' * Q) / landa;
end
% error of forgetting factor 
y_fit = intercept + slope1 * u1 + slope2 * u2;
error = y - y_fit;
figure;
plot(error , 'k' , 'LineWidth',2);
xlabel("data point index");
ylabel("error");
title('error plot of forgetting factor');
grid on;
%% sliding window
window_size = 20;
step = 1;
num_points = length(y);
num_windows = floor((num_points - window_size) / step) + 1;
interpcets = zeros(num_windows,1);
slope1 = zeros(num_windows,1);
slope2 = zeros(num_windows,1);
error = zeros(window_size , window_size);
for i = 1:num_windows
    strat_index = (i-1) * step +1;
    end_index = strat_index + window_size - 1;
    u1_window = u1(strat_index:end_index);
    u2_window = u2(strat_index:end_index);
    y_window = y(strat_index:end_index);
    U = [ones(size(u1_window)) u1_window u2_window];
    theta = pinv(U) * y_window;
    interpcets(i) = theta(1);
    slope1(i) = theta(2);
    slope2(i) = theta(3);
    y_fit = U *theta;
    error(i,:) = y_window - y_fit;
end
% plot
figure;

subplot(2,1,1);
hold on;
for i = 1:num_windows
    
    u1_window = u1((i - 1) * step + 1:(i - 1) * step + window_size);
    u2_window = u2((i - 1) * step + 1:(i - 1) * step + window_size);
    y_fit = interpcets(i) + slope1(i) * u1_window + slope2(i) * u2_window;
    plot3(u1_window , u2_window , y_fit , 'yellow');
end

scatter3(u1_window ,u2_window ,y_fit ,'black');
xlabel("Acid solution flow in liters");
ylabel('Base solution flow in liters');
zlabel('pH of the solution in the tank');
title("fitted line using sliding window least sqare method");
grid on;
hold off;

subplot(2 ,1 ,2);
plot(error','k','LineWidth',2);
xlabel('data point index');
ylabel('error');
title('error plot of the sliding window');
grid on;
%% RLS Method
window_size = 20;
step = 1;
interpcets = zeros(length(time_steps) - window_size + 1 , 1);
slope1 = zeros(length(time_steps) - window_size + 1 , 1);
slope2 = zeros(length(time_steps) - window_size + 1 , 1);
error = zeros(length(time_steps) - window_size + 1 , window_size);

for i=1:length(time_steps) - window_size + 1
    u1_window = u1(i:i+window_size-1);
    u2_window = u2(i:i+window_size-1);
    y_window = y(i:i+window_size-1);

    Q = eye(3);
    theta = zeros(3,1);
    error_window = zeros(1,window_size);
    for j = 1:window_size

        U = [1; u1_window(j); u2_window(j)];
        e = y_window(j) - U' *theta;
        K = (Q * U) / (1 + U' * Q * U);
        theta = theta + K * e;
        Q = (Q - K * U' * Q);
        error_window(j) = e;
    end
    error(i,:) = error_window;
end
figure;
subplot(2,1,1);
hold on;
for i= 1:length(time_steps) -window_size +1

    u1_window = u1(i:i+window_size-1);
    u2_window = u2(i:i+window_size-1);
    y_fit = interpcets(i) + slope1(i) * u1_window + slope2(i) * u2_window;
    plot3(u1_window, u2_window, y_fit, 'red');
end
scatter3(u1, u2, y,'filled');
xlabel("Acid solution flow in liters");
ylabel('Base solution flow in liters');
zlabel('pH of the solution in the tank');
title("fitted line using RLS method");
grid on;
hold off;
subplot(2 ,1 ,2);
plot(error','k','LineWidth',2);
xlabel('data point index');
ylabel('error');
title('error plot of the RLS');
grid on;
