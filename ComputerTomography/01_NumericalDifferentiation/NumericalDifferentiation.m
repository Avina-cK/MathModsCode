%% Exercise Sheet 1
%{
Implement differentiation method.
g = exp(-1,1)
noisy g = gdelta = g+ normal distributed noise with mean 0 and std 0.01 to
g
1. Find h+ by visual comparison with true derivative g'
2. Plot g' and D_(h+)(gdelta)

%}

%% Preliminary 
divz = 0.01;
x = -1.0:divz:1.0;
g_vec = exp(x);

%Creating noise : std*randn() + mean
rng(31415);
noise = 0.01*randn(1,length(g_vec)) + 0.0;
% disp([mean(noise), std(noise)]);

%noisy data
gdelta = g_vec+noise;

%Plot g and noisy g
hold off
plot(x,gdelta, "blue",'LineStyle', '--')
hold on
plot(x,g_vec,color="red")
hold off
legend('gdelta', 'g')

%% Part 1 : Finding h+
h_div = 1:1:50;
norms = zeros(size(h_div));
for j=h_div
    dgd_v = zeros(size(gdelta));
    h = divz*h_div(j);
    divtemp = 0:divz:1.0;
    h_i = find(divtemp==h,1);
    
    for i = 1:length(dgd_v)
        if ismember(i,1:h_i)
            dgd_v(i) = gdelta(i);
        elseif ismember(i,(length(dgd_v)-h_i):length(dgd_v))
            dgd_v(i)=gdelta(i);
        end
    end
    
    inter = (h_i+1):(length(g_vec)-h_i); 
    for i=inter
        dgd_v(i) = dg_d(gdelta,h,i);
    end
    
    norms(j) = norm(g_vec(inter)-dgd_v(inter));
end

%% Plotting with optimal h: h+
dgd_v = zeros(size(gdelta));
    h = divz*30;
    divtemp = 0:divz:1.0;
    h_i = find(divtemp==h,1);
    
    for i = 1:length(dgd_v)
        if ismember(i,1:h_i)
            dgd_v(i) = gdelta(i);
        elseif ismember(i,(length(dgd_v)-h_i):length(dgd_v))
            dgd_v(i)=gdelta(i);
        end
    end
    
    inter = (h_i+1):(length(g_vec)-h_i); 
    for i=inter
        dgd_v(i) = dg_d(gdelta,h,i);
    end
    
hold off
plot(x, g_vec, "red")
hold on
plot(x, dgd_v, "green")

%%
disp('done')
%% Functions
function dgd = dg_d(gdelta, h,z)
    divz = 0.01;
    x1 = int32(z+((h/divz)+1));
    x2 = int32(z-((h/divz)+1));
    dgd = (1/(2*h))*(gdelta(x1) - gdelta(x2));
end



