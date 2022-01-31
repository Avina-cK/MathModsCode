%% Exercise 4 : Tikhonov regularization and parameter selection
%{
A: L^2 (0,1) -> L^2 (0,1)
Af(x) = int(0,x)(f(t)dt)

f1(x) = sign(x-0.5)
f2(x) = sin(pi*x)

%}

%% Part a
 % Run this section to view the plot of f, g and noisy g
 
N = 300;
h = 1/N;
i = 1.0:1.0:N;
% creating grid as per Exercise Sheet 6
x = (i-0.5)*h;

% Create A as the sum of a low triangular matrix and a diagonal matrix
A = h*(tril(ones(N,N), -1) + diag(0.5*ones(N, 1)));
% De-comment the next line to view the matrix
% heatmap(A, 'GridVisible','off'); title('Matrix A')

f2_x = f2(x);
subplot(1,2,1)
% Plot f2 to see the function
plot(x,f2_x, 'DisplayName','f2')
xlabel('x')
ylabel('f2')
title('Original function')
legend

% Create g = Af
g2 = A*f2_x';

% seed to make sure noise is the same every time the code is run
rng(314159)
i = 1:1:length(g2);
d1 = datasample(i,length(g2)/2);
d3 = ones(length(g2),1);
d3(d1) = -1;

% Create noisy g
g2_d = g2 + (0.05 .* d3 .*g2);

% Checking if there is exactly 5% noise in g
if norm(g2_d-g2)/norm(g2)==0.05
    disp('5% noise in g')
end

% Plotting g2 and noisy g2
subplot(1,2,2)
hold off
plot(x,g2, 'DisplayName', 'g_2 = A*f2')
hold on
plot(x,g2_d, 'DisplayName','g_2^d = g2 + noise')
xlabel('x')
ylabel('g2')
legend('Location','southeast')

%% Part b : implement the Tikhonov-regularization
% Run this section to view plots that show implemented Tikhonov-regularization
% The regularization method is implemented as the function Tk_reg below

% Implementing the Tikhonov-regularization on g2 for various alpha
hold off
subplot(1,2,1)
scatter(x,f2_x, 5, ...
    'MarkerEdgeColor',[0,0,0], ...
    'MarkerEdgeAlpha', 0.2, ...
    'DisplayName', 'Original function')
hold on
plot(x, Tk_reg(A, 0, g2), ...
    'DisplayName', 'alpha=0')
hold on
plot(x, Tk_reg(A, 1, g2), ...
    'DisplayName', 'alpha=1')
hold on
plot(x, Tk_reg(A, 5, g2), ...
    'DisplayName', 'alpha=5')
title('A^+g2 (no noise)')
legend('Location','northwest')
xlabel('x')
ylim([0.0 1.22])

% Implementing the Tikhonov-regularization on noisy g2 for various alpha
hold off
subplot(1,2,2)
scatter(x,f2_x, 5, ...
    'MarkerEdgeColor',[0,0,0], ...
    'MarkerEdgeAlpha', 0.2, ...
    'DisplayName', 'Original function')
hold on
plot(x, Tk_reg(A, 0.0001, g2_d), ...
    'DisplayName', 'alpha=0.0001', ...
    color=[1,0.75,0])
hold on
plot(x, Tk_reg(A, 0.002, g2_d), ...
    'DisplayName', 'alpha=0.002', ...
    color='blue')
hold on
plot(x, Tk_reg(A, 0.01, g2_d), ...
    'DisplayName', 'alpha=0.01', ...
    color=[1,0,0])
title('A^+g^d (5% noise)')
ylim([0.0 1.22])
xlabel('x')
legend('Location','south')

sgtitle('Tikhonov-regularization')

%% Part c : Morozov 


% slowly monotonically decreasing function : 1/ln(x)
y = 1.5:0.01:2.5;
y = (log(y).^(-1));
y = y ./200;

%{
Let {t_k}_(k in N) be a strictly monotone decreasing sequence and tau>1
fixed. Determine k* such that
        norm(Afd_(t_k*) - gd) <= tau*delta < norm(Afd_(t_i) - gd)
gamma = t_k*
%}

t = y;
tau = 1.1;
delta = 0.4;
barrier = tau*delta*ones(length(t), 1);
% g_delta  = g2_d is defined above
% Tk_reg(A, t_k*, g_delta) is the function that gives fd_(t_k*)

F = zeros(length(f2_x) , length(t));
for i=1:length(t)
     f_temp = Tk_reg(A, t(i), g2_d);
     F(:,i) = f_temp;
end

Af = zeros(length(g2_d), length(t));
for i=1:length(t)
    F_temp = F(:,i);
    Af_temp = A*F_temp;
    Af(:,i) = Af_temp;
end

norms = zeros(length(t), 1);
for i=1:length(t)
    norms(i, 1) = norm(Af(:,i) - g2_d);
end

%Plot norms to find optimal alpha
hold off
plot(barrier, 'DisplayName', 'tau*delta')
hold on
plot(norms, 'DisplayName', 'norm(Afd_(t_k) - gd)')
legend
xlabel('k')
 
 % From the curve, we see that the optimal alpha is t(52)~ 0.007
 op_alpha = t(52);
 hold off
 plot(x, Tk_reg(A, op_alpha, g2_d), ...
     'DisplayName', 'alpha(Morozov) ~ 0.007', ...
     color=[1,0,0])
 hold on 
 plot(x,f2_x,  ...
     'DisplayName', 'Original function', ...
     color=[0,0,1])
 hold on
 plot(x, Tk_reg(A, 0.002, g2_d), ...
     'DisplayName', 'alpha(visual inspection) = 0.002', ...
     color=[0,1,0])
 legend('Location','south')
 xlabel('x')
 ylabel('f')

%% Part d (a)

f1_x = f1(x);
subplot(1,2,1)
plot(x,f1_x, 'DisplayName','f1')
legend('Location','northwest')
title('Original function')
xlabel('x')
ylabel('f1')
ylim([-1.1 1.1])

g1 = A*f1_x';

% seed to make sure noise is the same every time the code is run
rng(314)
i = 1:1:length(g1);
d1 = datasample(i,length(g1)/2);
d3 = ones(length(g1),1);
d3(d1) = -1;
g1_d = g1 + (0.05 .* d3 .*g1);
if norm(g1_d-g1)/norm(g1)==0.05
    disp('5% noise in g')
end

subplot(1,2,2)
hold off
plot(x,g1, 'DisplayName', 'g_1 = A*f_1')
hold on
plot(x,g1_d, 'DisplayName','g_1^d = g_1 + noise')
legend('Location','southeast')
xlabel('x')
ylabel('g1')
ylim([-0.65 0.0])

%% Part d (b)
hold off
subplot(1,2,1)
scatter(x,f1_x, 5, ...
    'MarkerEdgeColor',[0,0,0], ...
    'MarkerEdgeAlpha', 0.2, ...
    'DisplayName', 'Original function')
hold on
plot(x, Tk_reg(A, 0, g1), ...
    'DisplayName', 'alpha=0')
hold on
plot(x, Tk_reg(A, 0.005, g1), ...
    'DisplayName', 'alpha=0.005')
hold on
plot(x, Tk_reg(A, 0.05, g1), ...
    'DisplayName', 'alpha=0.05')
title('A^+g_1 (no noise)')
legend('Location','northwest')
ylim([-1.1 1.1])
xlabel('x')

hold off
subplot(1,2,2)
scatter(x,f1_x, 5, ...
    'MarkerEdgeColor',[0,0,0], ...
    'MarkerEdgeAlpha', 0.2, ...
    'DisplayName', 'Original function')
hold on
plot(x, Tk_reg(A, 0.0001, g1_d), ...
    'DisplayName', 'alpha=0.0001', ...
    color=[1,0.75,0])
hold on
plot(x, Tk_reg(A, 0.001, g1_d), ...
    'DisplayName', 'alpha=0.001', ...
    color='blue')
hold on
plot(x, Tk_reg(A, 0.01, g1_d), ...
    'DisplayName', 'alpha=0.01', ...
    color=[1,0,0])
title('A^+g_1^d (5% noise)')
ylim([-1.5 1.5])
legend('Location','northwest')
xlabel('x')

%% Part d (c)

% slowly monotonically decreasing function : 1/ln(x)
y = 1.02:0.005:1.5;
y = (log(y).^(-1));
y = y ./2000;

%{
Let {t_k}_(k in N) be a strictly monotone decreasing sequence and tau>1
fixed. Determine k* such that
norm(Afd_(t_k*) - gd) <= tau*delta < norm(Afd_(t_i) - gd)
gamma = t_k*
%}

t = y;
tau = 1.000001;
delta = 0.4;
barrier = tau*delta*ones(length(t), 1);
% g_delta  = g2_d is defined above
% Tk_reg(A, t_k*, g_delta) is the function that gives fd_(t_k*)

F = zeros(length(f1_x) , length(t));
for i=1:length(t)
     f_temp = Tk_reg(A, t(i), g1_d);
     F(:,i) = f_temp;
end

Af_t = zeros(length(g1_d), length(t));
for i=1:length(t)
    F_temp = F(:,i);
    Af_temp = A*F_temp;
    Af(:,i) = Af_temp;
end

norms = zeros(length(t), 1);
for i=1:length(t)
    norms(i, 1) = norm(Af(:,i) - g1_d);
end

%Plot norms to find optimal alpha
hold off
plot(barrier, 'DisplayName', 'tau*delta')
hold on
plot(norms, 'DisplayName','norm(Afd_(t_k) - gd)')
xlabel('k')
legend


% From the curve, we see that the optimal alpha is t(28)~0.003
op_alpha = t(27);

hold off
plot(x, Tk_reg(A, op_alpha, g1_d), ...
    'DisplayName', 'alpha(Morozov) ~ 0.003', ...
    color=[1,0,0])
hold on 
plot(x,f1_x,  ...
    'DisplayName', 'Original function', ...
    color=[0,0,1])
hold on
plot(x, Tk_reg(A, 0.001, g1_d), ...
    'DisplayName', 'alpha(visual inspection) = 0.002', ...
    color=[0,1,0])
legend('Location','northwest')
xlabel('x')
ylabel('f')


%% Functions

function ans_f1 = f1(x)
    ans_f1 = sign(x - 0.5);
end

function ans_f2 = f2(x)
    ans_f2 = sin(pi*x);
end

function Aplusg = Tk_reg(A, alpha, g2)
    [U, S, V] = svd(A);
    N = length(A);
    imp = zeros(N,N);
    sigma = diag(S);
    
    j = 1;
    K = N;
    while j<=K
        imp(:,j) = (sigma(j)/((sigma(j)^2) +alpha)) * dot(g2, U(:, j)) * V(:,j);
        j=j+1;
    end
    Aplusg = sum(imp,2);
end