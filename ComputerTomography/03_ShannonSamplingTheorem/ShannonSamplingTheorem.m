
% d is assumed to be d = pi/L 
%% Part (a)
%{ 
Band limited function:
    f(t)=(4/pi)*((sin(2t)/(2t))^2)*cos(2t)
Fourier transform:
                0.25*(v+6)  -6<=v<=-2
    Ff(v) = {   1           |v|<2
                0.25*(v-6)  2<=v<=6
                0           |v|>-6

%}

L = 100*pi;
d = pi/L;

t = -8:d:8;
pf = zeros(1,length(t));
for i=1:length(pf)
   pf(i) = f(t(i));
end

pSn = zeros(1,length(t));
for i=1:length(pSn)
   pSn(i) = S(10000, t(i), L);
end

figure
subplot(1,4,1)
plot(t, pf)
xlabel('t')
ylabel('f(t)')
subplot(1,4,3)
plot(t, pSn)
xlabel('t')
ylabel('S(10000)(t;f)')
title('L=100*pi')

L = 1000*pi;
d = pi/L;

t = -8:d:8;
pf = zeros(1,length(t));
for i=1:length(pf)
   pf(i) = f(t(i));
end

pSn = zeros(1,length(t));
for i=1:length(pSn)
   pSn(i) = S(10000, t(i), L);
end
subplot(1,4,4)
plot(t, pSn)
xlabel('t')
ylabel('S(10000)(t;f)')
title('L=1000*pi')

L = 50*pi;
d = pi/L;

t = -8:d:8;
pf = zeros(1,length(t));
for i=1:length(pf)
   pf(i) = f(t(i));
end

pSn = zeros(1,length(t));
for i=1:length(pSn)
   pSn(i) = S(10000, t(i), L);
end
subplot(1,4,2)
plot(t, pSn)
xlabel('t')
ylabel('S(10000)(t;f)')
title('L=50*pi')

% Higher the L, higher the resolution of the approximation.


%% Part b
L = 2*pi;
d = pi/L;

t = -3:d:3;
pf = zeros(1,length(t));
for i=1:length(pf)
   pf(i) = fb(t(i));
end

pSn = zeros(1,length(t));
for i=1:length(pSn)
   pSn(i) = Sb(10000, t(i), L);
end

figure

subplot(1,2,1)
plot(t, pf)
xlabel('t')
ylabel('f(t)')
subplot(1,2,2)
plot(t, pSn)
xlabel('t')
ylabel('S(10000)(t;f)')

% With L = 2pi, nothing is plotted.


%% Part c
L = 100*pi;
d = pi/L;

t = -3:d:3;
pf = zeros(1,length(t));
for i=1:length(pf)
   pf(i) = fc(t(i));
end

pSn = zeros(1,length(t));
for i=1:length(pSn)
   pSn(i) = Sc(10000, t(i), L);
end

figure

subplot(1,4,1)
plot(t, pf)
xlabel('t')
ylabel('f(t)')
subplot(1,4,4)
plot(t, pSn)
xlabel('t')
ylabel('S(10000)(t;f)')
title('L=100*pi')


L = 50*pi;
d = pi/L;

t = -3:d:3;
pf = zeros(1,length(t));
for i=1:length(pf)
   pf(i) = fc(t(i));
end

pSn = zeros(1,length(t));
for i=1:length(pSn)
   pSn(i) = Sc(10000, t(i), L);
end

subplot(1,4,3)
plot(t, pSn)
xlabel('t')
ylabel('S(10000)(t;f)')
title('L=50*pi')


L = 10*pi;
d = pi/L;

t = -3:d:3;
pf = zeros(1,length(t));
for i=1:length(pf)
   pf(i) = fc(t(i));
end

pSn = zeros(1,length(t));
for i=1:length(pSn)
   pSn(i) = Sc(10000, t(i), L);
end

subplot(1,4,2)
plot(t, pSn)
xlabel('t')
ylabel('S(10000)(t;f)')
title('L=10*pi')

% Higher the L, higher the resolution of the approximation.


%% Functions

function ans_f = f(t)
    if t==0
        ans_f = 4/pi;
    else
        ans_f = (4/pi)*((sin(2*t)/(2*t))^2)*(cos(2*t));
    end
end

function ans_fb = fb(t)
    ans_fb = exp(-(t^2));
end

function ans_fc = fc(t)
    if t>=(-1)
        if t<=1
            ans_fc = 1;
        end
    end
    if t<(-1)
        ans_fc = 0;
    end
    if t>1
        ans_fc = 0;
    end
end

function ans_Ff = Ff(v)
    if v>=-6
       if v<=-2
          ans_Ff = 0.25*(v+6); 
       end
    end
    
    if abs(v)<2
        ans_Ff = 1;
    end
    
    if v>=2
        if v<=6
            ans_Ff = 0.25*(v-6);
        end
    end
    
    if abs(v)>6
        ans_Ff = 0;
    end
end

function ans_Sn = S(n,t, L)
    Snseq = zeros(length((2*n)));
    for i=-n:(n)
       
           Snseq(i+n+1) = f(i*pi/L)*(sin((L*t)-(i*pi))/((L*t)-(i*pi)));

    end
    ans_Sn = sum(Snseq);
end


function ans_Snb = Sb(n,t, L)
    Snseq = zeros(length((2*n)));
    for i=-n:(n)
           Snseq(i+n+1) = fb(i*pi/L)*(sin((L*t)-(i*pi))/((L*t)-(i*pi)));
    end
    ans_Snb = sum(Snseq);
end

function ans_Snc = Sc(n,t, L)
    Snseq = zeros(length((2*n)));
    for i=-n:(n)
       
           Snseq(i+n+1) = fc(i*pi/L)*(sin((L*t)-(i*pi))/((L*t)-(i*pi)));
       
    end
    ans_Snc = sum(Snseq);
end
