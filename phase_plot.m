
clc
clear

fs = 512;
time = 1000;

%% sin
t = linspace(0, time, fs);
x = sin(t);

delay = 3;
plot(x)

%% Lorenz
% 定义 Lorenz 系统参数
sigma = 10;
rho = 28;
beta = 8/3;

% 定义 Lorenz 系统的微分方程
lorenz_eq = @(t, X) [sigma * (X(2) - X(1)); 
                      X(1) * (rho - X(3)) - X(2); 
                      X(1) * X(2) - beta * X(3)];

% 设置时间范围
t_span = [0 100];
t_eval = linspace(0, 100, fs);  % 采样点

% 初始条件
initial_state = [-13, -14, 47];

% 求解微分方程
[t, sol] = ode45(lorenz_eq, t_eval, initial_state);
x = sol(:, 1)';

delay = 3;
plot(x)

%% Rossler
% 定义 Rossler 系统参数
alpha = 0.2;
rho = 5.7;
beta = 0.2;

% 定义 Rossler 系统的微分方程
rossler_eq = @(t, X) [-X(2) - X(3); 
                      X(1) + alpha * X(2); 
                      beta + X(3) * (X(1) - rho)];

% 设置时间范围
t_span = [0 100];
t_eval = linspace(0, 100, fs);  % 采样点

% 初始条件
initial_state = [-2, -10, 0.2];

% 求解微分方程
[t, sol] = ode45(rossler_eq, t_eval, initial_state);
x = sol(:, 1)';

delay = 18;
plot(x)

%% Logistic
r = 3.8;
x = zeros(1, fs);
x(1) = 0.5;
for i = 2:fs
    x(i) = r * x(i - 1) * (1 - x(i - 1));
end

delay = 1;
plot(x)

%% Henon
a = 1.4;
b = 0.3;
x = zeros(1, fs);
y = zeros(1, fs);
x(1) = -0.0064;
y(1) = -0.4735;
for i = 2:fs
    x(i) = 1 - a * x(i - 1)^2 + y(i - 1);
    y(i) = b * x(i - 1);
end

delay = 1;
plot(x)

%% test
noise = randn(fs, 1) * 0.2;

tau = 1;
n_cluster = 1;
threshold_corr = 0.85;
threshold_nmse = 0.01; 
threshold_svd = 0.96;
periodic_extension = false;

sig = normalize(x(1, :), 'range').' + noise;
SGCs = SpSSA(sig, fs, tau, n_cluster, ...
    threshold_corr, threshold_nmse, threshold_svd, ...
    'kmeans', periodic_extension);

% SGSs = FMD(fs, sig, 30, 1000, 5, 1000);
% [SGCs,residual] = emd(sig, MaxNumIMF=2);
% SGCs = SGCs';

relative_fun = @(x, y) norm(y-x, 2) / norm(x, 2);
fprintf("%.6f\n", relative_fun(x', SGCs(1, :)));
fprintf("%.6f\n", relative_fun(noise, SGCs(2, :)));

figure();
subplot(1, 2, 1)
sig = normalize(x(1, :), 'range');
scatter(sig(1:end-delay), sig(1+delay:end))
subplot(1, 2, 2)
sig = normalize(SGCs(:, 1), 'range');
scatter(sig(1:end-delay), sig(1+delay:end))
