% v=-1:0.01:1;
% H = abs((1-exp(sqrt(-1)*2*pi*v))./(1-exp(sqrt(-1)*2*pi*v/256)));
% plot(v,H)
% clear all
% clc
% close all;
% v = 25;
% N=256;
% H = zeros(1,101);
% for m = -50:1:50
%     for n = 1:N
%         H(m+51) = H(m+51) + exp(sqrt(-1)*2*pi*(v-m)*n/N);
%     end
% end
% figure 
% plot(-50:1:50,abs(H));
clear all;
close all;
clc
%==========获取OFDM系统仿真配置===========
config = OFDMSystemConfig;

simu_times = config.simu_times;%仿真次数
N = config.N;%OFDM符号周期长度
Ng = config.Ng;%循环前缀长度
ffo_df = 0;%小数倍频偏
ifo_df = 50;%整数倍频偏
SNR = config.SNR;%信噪比
%============构建训练符号=================
%[transmit_data,TS] = generate_TS_method1(N,Ng);
TS = generate_TS(N);
transmit_data = TS;
%===========添加频偏===========
df =ffo_df + ifo_df;%总频偏
recv_sig = zeros(1,length(transmit_data));
%Shao_max_offset = zeros(1,length(SNR));%记录最大偏差率。
for k=1:length(transmit_data)
    %添加频偏
    recv_sig(k) = transmit_data(k)*exp(sqrt(-1)*2*pi*df*k/N);
end    
recv_sig = awgn(recv_sig,-5,'measured');
 v = zeros(1,N+1);
% for k=1:N/32
%     v(k) = N*angle(recv_sig(k)*TS(k))/(2*pi*k);
% end
% vv = sum(abs(v)) / (N/32)
for  k = 1:N
    Q1 = recv_sig(k)*conj(TS(k));
    if k+1 >N
        Q2 = recv_sig(1)*conj(TS(1));
    else
        Q2 = recv_sig(k+1)*conj(TS(k+1));
    end
    v(k) = N*angle(conj(Q1)*Q2)/(2*pi);
end
vv = sum(v)/(N);
if vv >= 0
    v1 = floor(vv);
    df = vv - v1;
    v2 = v1-N;
else
    v1 = ceil(vv);
    df = vv - v1;
    v2 = vv+N;
end

W = N/8;
P = zeros(1,2*W+1);
E =0 ;
r = recv_sig;
for n = 1:N
    E = E + abs(r(n))^2;%计算信号能量
end
index = 0;
for k = v1 - W : v1+W
    index = index+1;
    for n = 1:N
        P(index) = P(index) + r(n)*conj(TS(n))*exp(-sqrt(-1)*2*pi*k*n/N);
    end
end
H = abs(P).^2/E;
[maxvalue1,index] = max(H);
ifo1 = v1 - W + index - 1;
figure
plot(v1 - W : v1+W,H);
hold on


index = 0;
P = zeros(1,2*W+1);
for k = v2 - W : v2+W
    index = index+1;
    for n = 1:N
        P(index) = P(index) + r(n)*conj(TS(n))*exp(-sqrt(-1)*2*pi*k*n/N);
    end
end
H = abs(P).^2/E;
[maxvalue2,index] = max(H);
ifo2 = v2 - W + index - 1;
plot(v2 - W : v2+W,H);

if maxvalue1>maxvalue2
    ifo = ifo1;
else
    ifo = ifo2;
end
ifo
df
v1
v2
H = zeros(1,2*N+1);
for m=-N:N
    for n=1:N
        H(m+N+1) = H(m+N+1) + recv_sig(n)*conj(TS(n))*exp(-sqrt(-1)*2*pi*m*n/N);
    end
end
M= (abs(H).^2)./E;
[~,df] = max(M);
df = df -N-1;
figure
plot(-N:N,M);

function [TS] = generate_TS(N)
a = zeros(1,N/4);
mu = N/4-1;
for n=0:N/4-1
   a(n+1) =7*exp(4*sqrt(-1)*mu*pi*n*n/N);
end
A = ifft(a);
B = conj(A(1,N/4:-1:1));
C= zeros(1,N/4);
for n=1:1:N/4 
    if mod(n,2)
       C(n) = (-1)*B(n);
    else
       C(n) = B(n);
    end
end
D=C;
signal = [A B C D];
TS =signal;%时域上的训练符号
end