%% Discretization of the convolution and filtering
%{
Let k in L2 ([0, 1]) real and periodic continued on R. 
The periodic convolution K : L2 ([0; 1])-> L2 ([0; 1]) is defined as
Kf(x):=(f*k)(x)=int_0^1 {f(y) k(x-y) dy}
%}

%% Part a and b:
%{
Discretize the convolution operator K using the points 
        y_j = j/N ; j = 0,... ,N-1 
and a box rule for approximating the integral such that the
convolution is approximated by matrix-vector product of the form 
        Kf=g 
with convolution matrix K and sampled vector f = (f0, f1, ... , f(N-1)). 
The matrix K is a Toeplitz matrix if x is discretized at the same grid as y
%}
N=500;
j = 0:1:N-1;
y = j/N;
x = y';
%{
For k, convolution matrix K=toeplitz(k). The convolution is defined as the 
function MConvolution(k,f)
%}

%% Part c:
b = (1/10);
f = cos(4*pi*x);
k = kq(x,b);
g =  MConvolution(k, f_delta);
hold off
plot(x,g)

%% Part d: 
% Add noise and denoise f

noise = 0.05*randn(size(f));
f_delta = f + noise;
b = (1/10);
k = kq(x,b);
g_delta = MConvolution(k, f_delta);
hold off
plot(x,g_delta);

%% part (e)
%{
f(x) = window function_T(x-0.5); T= 0.25
%}
b = (1/10);
T = 0.25;
f_sq = zeros(size(x));
for i=1:length(f_sq)
    a = x(i)-0.5;
    if a<T && a>-T
        f_sq(i) = 1;
    %elseif a==T || a==-T
     %   f_sq(i) = 0.5;
    else
        f_sq(i) = 0.0;
    end
end

g_sq = MConvolution(k, f_sq');
hold off
plot(x,g_sq/max(g_sq))
hold on
plot(x,f_sq)
legend('g_sq=Kf_sq', 'f_sq=window_0.25(x-0.5)', 'Location', 'south')

% With noise
f_sq_n = f_sq + 0.05*randn(size(f_sq));
g_sq_n = MConvolution(k, f_sq_n');

hold off
plot(x, f_sq_n,'LineStyle','--')
hold on
plot(x, g_sq_n/max(g_sq_n),'LineStyle','-')
legend('f_sq (with noise)', 'g_sq_n/max(g_sq_n)', 'Location','south')

%% Functions

function g = MConvolution(k_giv, f)
    K = toeplitz(k_giv);
    g = K*f;
end

% indicator function of the interval [0,b_giv)
function k_ans = kq(x_giv,b_giv)
   k_b = zeros(size(x_giv));
   for i = 1:length(x_giv)
       if x_giv(i)>=0 && x_giv(i)<b_giv
            k_b(i) = 1;
       else
           k_b(i) = 0;
       end
   end
   k_ans = k_b;
end