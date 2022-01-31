%% Group:
    % Kalle, Tadinada, Ravindran
%The optimal value of K was found to be 15. The minimal norm solution gives
%highly oscillating results.
%% Preparing the system

N = 100;
h = 1/N;
i = 1.0:1.0:N;
x = (i-0.5)*h;

% Create A as the sum of a low triangular matrix and a diagonal matrix
A = tril(ones(N,N), -1) + diag(0.5*ones(N, 1));

%% Part a

% define f as given
f = (exp(-2.0 .*x)) .*cos(5.0 .*x);
f = f';
% noise
rng(3141592)
noise = 0.1*normrnd(0, 0.19, N, 1);

% Computing g ; g = Af
g = A*f;

%Computing noisy g = g_d
g_d = g + noise;
% Computing norm (g - g_d) and checking if it is less than or equal to 0.19
if norm(g - g_d)<=0.19
    disp('||g - g_d|| <= 0.19');
end

%% Implementing the truncated singular value decomposition
[U, S, V] = svd(A);
imp = zeros(N,N);
sigma = diag(S);

% Without noise
j = 1;
K = N;
while j<=K
    imp(:,j) = (1.0/sigma(j)) * dot(g, U(:, j)) * V(:,j);
    j=j+1;
end

Aplusg = sum(imp,2);

% With noise
imp = zeros(N,N);
K = 15;
j = 1;
while j<=K
    imp(:,j) = (1.0/sigma(j)) * dot(g_d, U(:, j)) * V(:,j);
    j=j+1;
end

Aplusg_d15 = sum(imp,2);

%Different K
imp = zeros(N,N);
K = 1;
j = 1;
while j<=K
    imp(:,j) = (1.0/sigma(j)) * dot(g_d, U(:, j)) * V(:,j);
    j=j+1;
end

Aplusg_d1 = sum(imp,2);

% K = 3
imp = zeros(N,N);
K = 3;
j = 1;
while j<=K
    imp(:,j) = (1.0/sigma(j)) * dot(g_d, U(:, j)) * V(:,j);
    j=j+1;
end

Aplusg_d3 = sum(imp,2);

% K = 90
imp = zeros(N,N);
K = 90;
j = 1;
while j<=K
    imp(:,j) = (1.0/sigma(j)) * dot(g_d, U(:, j)) * V(:,j);
    j=j+1;
end

Aplusg_d90 = sum(imp,2);

%% Plotting results
%{
hold off
subplot(1,2,1)
plot(x,g, DisplayName='g')
hold on
scatter(x,g_d, 10, DisplayName='g^d')
legend
%}

hold off
subplot(1,2,1)
plot(x, f, DisplayName='f')
hold on
scatter(x, Aplusg, 10, DisplayName='A(+)g', MarkerEdgeColor='c')
plot(x, Aplusg_d15, '-.',DisplayName='A(+)g_d (K=15)', color=[1,0,0], LineWidth =1)

% De-comment the next 2 lines to see the result of using different K values
% plot(x, Aplusg_d3, '-.',DisplayName='A(+)g_d (K=3)', color=	'#808080')
% plot(x, Aplusg_d90, '-.',DisplayName='A(+)g_d (K=90)', color='#C0C0C0')

legend

%% Minimal norm solution

invA = inv(A);
f_ti = invA*g_d;

hold off
subplot(1,2,2)
plot(x,f, DisplayName='f')
hold on
plot(x,f_ti,  DisplayName='inv(A)*g_d')
legend