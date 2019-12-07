close all; 
clear all; 
clc; 
%参数定义 
N=256;       %FFT/IFFT 变换的点数或者子载波个数（Nu=N） 
Ng=N/8;      %循环前缀的长度 (保护间隔的长度) 
Ns=Ng+N;     %包括循环前缀的符号长度 
SNR=15;
%SNR = [-10,-5,0,5,10,15,20,25];
QAMTable =[7+7i,-7+7i,-7-7i,7-7i]; 
deltaf = 0.3;
%QAMTable=[1+1i,-1+1i,-1-1i,1-1i]; 
%=============生成数据符号=============
src = QAMTable(randi([0,3],N,1)+1); 
sym = ifft(src); 
cp_sym=[sym(1,N-Ng+1:N) sym];
predata = cp_sym;
suffixdata = cp_sym;
%=============生成训练序列======================
buf=QAMTable(randi([0,3],N/2,1)+1);
x=zeros(1,N/2 ); 
index = 1; 
for n=1:2:N/2 
     x(n)=buf(index); 
     index=index+1; 
end; 
sch = ifft(x);   %[A A]的形式 
sch2=[sch (-1).*sch];%[A A -A -A]形式
cp_train = [sch2(1,N-Ng+1:N) sch2];
transmit_data = [predata cp_train suffixdata];  
%=============经过信道=========================
recv_sig = transmit_data;
for k=1:length(recv_sig)
    recv_sig(k) = recv_sig(k)*exp(i*2*pi*deltaf*(k-1)/N);
end
recv = awgn(recv_sig,SNR);
simulation_times = 5;
df = zeros(1,length(SNR));
df2 = zeros(1,length(SNR));
% for snr = 1:length(SNR) 
% for i=1:simulation_times
%     recv = awgn(recv_sig,SNR(snr));
    %============定时粗估计=====================
    [P,M] = minn(recv,N,Ng);
    [~,m] = max(M);
    
    %============小数倍频偏估计======================
    2*angle(P(m))/pi
    [mP,df] = CFO_estimate(recv,m,N,Ng) ;
%     2*angle(P(m))/pi
%     df(sn)= df(sn) + 2*angle(P(m))/pi;
%     df2(sn) = df2(sn) + power(2*angle(P(m))/pi-deltaf,2);
%end
%      df(sn)= df(sn)/simulation_times;
%      df2(sn)= df2(sn)/simulation_times;
%end
df
%mP
%df2
% d= 1:length(SNR) ;
% figure(1)
% plot(SNR,df(d))
% figure(2)
% plot(SNR,df2(d))

function [P,M] = minn(recv,N,Ng)
    Ns = N+Ng;
    P=zeros(1,2*Ns); 
    R=zeros(1,2*Ns);
    for d = Ns/2+1:1:2*Ns
    for k=1:2
        for m=0:1:N/4-1  
            P(d-Ns/2) = P(d-Ns/2) + conj(recv(d+m+(k-1)*N/2))*recv(d+N/4+(k-1)*N/2+m);  
            R(d-Ns/2) = R(d-Ns/2) + power(abs(recv(d+N/4+(k-1)*N/2+m)),2);
        end 
    end
    end 
    M=power(abs(P),2)./power(abs(R),2); 
end

function [P,df] = CFO_estimate(recv,m,N,Ng)
P=0; 
for d=m-Ng:m-1
    P = P + conj(recv(d)) * recv(d+N);
end
df = angle(P)/(2*pi);
% for d=0:N/2-1
%     P = P + conj(recv(d+m))*(-1)*recv(d+N/2+m);  
% end
% df = angle(P)/pi
end
%===========小数倍频偏补偿========================
% for k=1:length(recv)
%     recv(k) = recv(k)*exp(-i*2*pi*df*k/N);
% end
% 
% recv_preamble = recv(1,d:d+N-1);
% 
% recv_sig = fft(recv_preamble);

