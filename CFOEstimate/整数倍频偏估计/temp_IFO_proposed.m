clear all;
close all;
clc
%==========获取OFDM系统仿真配置===========
config = OFDMSystemConfig;
ffo_df = config.ffo_df;%小数倍频偏
ifo_df = config.ifo_df;%整数倍频偏
simu_times = config.simu_times;%仿真次数
N = config.N;%OFDM符号周期长度
Ng = config.Ng;%循环前缀长度
SNR = config.SNR;%信噪比
%============构建训练符号=================

%[transmit_data,TS] = generate_TS_method1(N,Ng);
[TS_TimeDomain,TS_FreDomain] = generate_TS(N,Ng);
transmit_data = TS_TimeDomain;
%===========添加频偏===========
df =ffo_df + ifo_df;%总频偏
recv_sig = zeros(1,length(transmit_data));
%Shao_max_offset = zeros(1,length(SNR));%记录最大偏差率。
for k=1:length(transmit_data)
    %添加频偏
    recv_sig(k) = transmit_data(k)*exp(sqrt(-1)*2*pi*df*k/N);
end    
recv_sig = awgn(recv_sig,0,'measured');
%ifo_df = get_IFO_method1(recv_sig,N,Ng,TS)
%ifo_fre_df1 = get_IFO_FreDomain(recv_sig,N,Ng,TS_FreDomain);%频域估计算法
ifo_time_df1 = get_IFO_TimeDomain1(recv_sig,N,TS_TimeDomain);%时域估计算法
ifo_time_df2 = get_IFO_TimeDomain2(recv_sig,N,TS_TimeDomain);%时域估计算法
ifo_time_df3 = get_IFO_TimeDomain3(recv_sig,N,TS_TimeDomain);%时域估计算法
%ifo_fre_df1
ifo_time_df1
ifo_time_df2
ifo_time_df3
function [df] = get_IFO_TimeDomain1(recv_sig,N,TS_TimeDomain)
P = zeros(1,N/2+1);
R = zeros(1,N/2+1);
d=1;
for m=-N/4:N/4
    for n=0:N-1
        if m<0
            v = m-0.25;
        else
            v = m+0.25;
        end
%     v=m;    
     P(m+N/4+1) = P(m+N/4+1) + recv_sig(d+n)*conj(TS_TimeDomain(n+1))*exp(-sqrt(-1)*2*pi*v*n/N);
     R(m+N/4+1) = R(m+N/4+1)+abs(recv_sig(d+n))^2;
    end
end
H = (abs(P).^2)./(R.^2);
[~,df] = max(H);
% figure
% plot(-N/4:N/4,abs(H));
% xlabel('Method1:整数频偏'); 
% ylabel('时域估计值'); 
df = df - N/4 - 1;
end
function [df] = get_IFO_TimeDomain2(recv_sig,N,TS_TimeDomain)
P = zeros(1,N/2+1);
R = zeros(1,N/2+1);
d=1;
for m=-N/4:N/4
    for n=0:N-1
    % P(m+N/4+1) = P(m+N/4+1) + recv_sig(d+n)*conj(TS_TimeDomain(n+1))*exp(-sqrt(-1)*2*pi*m*n/N);
    
    P(m+N/4+1) = P(m+N/4+1) + abs((recv_sig(d+n)-TS_TimeDomain(n+1)*exp(sqrt(-1)*2*pi*m*n/N)))^2;
    R(m+N/4+1) = R(m+N/4+1)+abs(recv_sig(d+n))^2;
    end
end
H = (P.^2)./(R.^2);
[~,df] = min(H);
% figure
% plot(-N/4:N/4,H);
% xlabel('Method2:整数频偏'); 
% ylabel('时域估计值'); 
df = df - N/4 - 1;
end
function [df] = get_IFO_TimeDomain3(recv_sig,N,TS_TimeDomain)
P = zeros(1,N/2+1);
R = zeros(1,N/2+1);

d=1;
for m=-N/4:N/4
    for n=0:N-1
     % H(m+N/4+1) = H(m+N/4+1) + recv_sig(d+n)*conj(TS_TimeDomain(n+1))*exp(-sqrt(-1)*2*pi*m*n/N);
      P(m+N/4+1) = P(m+N/4+1) + abs((recv_sig(d+n)-TS_TimeDomain(n+1)*exp(sqrt(-1)*2*pi*m*n/N)))^2;
      R(m+N/4+1) = R(m+N/4+1)+abs(recv_sig(d+n))^2;
    end
end
M2 = (P.^2)./(R.^2);
[~,df] = min(M2);
df = df - N/4 - 1;
H = zeros(1,5);
R = zeros(1,5);
for m=df-2:df+2
        if m<0
            v = m-0.25;
        else
            v = m+0.25;
        end
    for n=0:N-1
        H(m-(df-2)+1) = H(m-(df-2)+1) + recv_sig(d+n)*conj(TS_TimeDomain(n+1))*exp(-sqrt(-1)*2*pi*v*n/N);
        R(m-(df-2)+1) = R(m-(df-2)+1)+abs(recv_sig(d+n))^2;
    end
end
M1= (abs(H).^2)./R;

%M=M1-M2;
[~,df2] = max(M1);
% figure
% plot(df-2:df+2,M1);
% hold on

%plot(-N/4:N/4,M);
% xlabel('Method3:整数频偏'); 
% ylabel('时域估计值'); 
df = df2 +(df-2)-1;
end
function [ifo_df] = get_IFO_FreDomain(recv_sig,N,Ng,TS_FreDomain)
    d=1;
    r = recv_sig(1,d:d+N-1);
    R = fft(r);
    H = zeros(1,N+1);
    D= zeros(1,N);
    E= zeros(1,N);
    for m=-N/2:N/2
        if m>=0
            Rs = [TS_FreDomain(1,N-m+1:N) TS_FreDomain(1,1:N-m)];
        end
        if m<0
            n = abs(m);
            Rs = [TS_FreDomain(1,n+1:N) TS_FreDomain(1,1:n)];
        end
        for k=1:N
            D(k) = abs(Rs(k) - R(k))^2; 
            E(k) = abs(R(k))^2;
        end
        H(m+N/2+1) = sum(D)^2/sum(E)^2;
    end
    figure
    plot(-N/2:N/2,abs(H));
    xlabel('整数频偏'); 
    ylabel('频域估计值'); 
    [~,ifo_df] = min(H);
    ifo_df = ifo_df-N/2-1;
end
function [TS_TimeDomain,TS_FreDomain] = generate_TS(N,Ng)
a = zeros(1,N/4);
mu = N/4-1;
for n=0:N/4-1
   a(n+1) =exp(4*sqrt(-1)*mu*pi*n*n/N);
end
A = ifft(a);
%A =a;
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
TS_TimeDomain =signal;%时域上的训练符号
TS_FreDomain = fft(signal);%对应频域上的训练符号
end
%TS = [train_symbol(1,N-Ng+1:N) train_symbol];

% function [ifo_df] = get_IFO_method1(recv_sig,N,Ng,TS)% Cover 韩留斌
%     d = 1;
%     r = recv_sig(1,d:d+N-1);
%     R = fft(r);
%     for m=0:10
%         Rs = [TS(1,N-m+1:N) TS(1,1:N-m)];
%         for k=1:N
%             D(k) = abs(Rs(k) - R(k))^2; 
%             E(k) = abs(R(k))^2;
%         end
%         H(m+1) = sum(D)^2/sum(E)^2;
%     end
% 
%     plot(0:10,H);
%     [~,ifo_df] = min(H);
%     ifo_df = ifo_df-1;
% end
% function [train_symbol,ts] = generate_TS_method1(N,Ng)
% a = zeros(1,N);
% for n=0:N-1
%     a(n+1) = exp(sqrt(-1)*pi*n*n/N);
% end
% ts = a;
% train_symbol = ifft(a);
% transmit_data = train_symbol;
% %TS = [train_symbol(1,N-Ng+1:N) train_symbol];
% end