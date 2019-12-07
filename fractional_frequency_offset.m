close all; 
clear all; 
clc; 
%参数定义 
N=256;       %FFT/IFFT 变换的点数或者子载波个数（Nu=N） 
Ng=N/8;      %循环前缀的长度 (保护间隔的长度) 
Ns=Ng+N;     %包括循环前缀的符号长度 
SNR=15;
QAMTable =[7+7i,-7+7i,-7-7i,7-7i]; 
%QAMTable=[1+1i,-1+1i,-1-1i,1-1i]; 
src = QAMTable(randi([0,3],N,1)+1); 
sym = ifft(src); 
cp_sym=[sym(1,N-Ng+1:N) sym];
predata = cp_sym;
suffixdata = cp_sym;
[transmit_data_cazac,] = generate_data(N,Ng,predata,suffixdata);
delta = -0.3;
freq_offset_item = exp(sqrt(-1)*pi*delta/N);
for i=0:3*Ns-1
    transmit_data_cazac(i+1) = transmit_data_cazac(i+1) *  exp(sqrt(-1)*2*pi*i*delta/N);
end
if SNR<100
    recv = awgn(transmit_data_cazac ,SNR);
else
    recv = transmit_data_cazac;
end
M = cazac( recv,N,Ng);
d = 1:length(M);
plot(d,M);
[~,loc] = max(M);
delta =estimate_freq_offset(recv,N,Ng,loc);
delta


function [delta] = estimate_freq_offset( recv_data,N,Ng,d)
Ns = N+Ng;
%P=zeros(1,N/4);
% for i = 0:1:N/4-1
%     for m=0:N/4-1
%         P(i+1) = P(i+1) + recv_data(d+m)*conj(recv_data(d+N/2-1-m));
%     end
% end
% [~,b] = max(abs(P));
P = 0;
for m=0:N/4-1
    P = P + recv_data(d+m)*conj(recv_data(d+N/2-1-m));
end
P
delta = angle(P)/pi;
end

function [trans_data] = generate_data(N,Ng,predata,suffixdata,delta)

    a = zeros(1,N/4);
    for n=0:N/4-1
       a(n+1) =exp(4*sqrt(-1)*pi*n*n/N);
    end
    A =a;
    B = A(1,N/4:-1:1);
    C= zeros(1,N/4);
    for n=1:1:N/4
        if mod(n,2)
            C(n) = (-1)*conj(B(n));
        else
            C(n) = conj(B(n));
        end
    end
   D = conj(C); 
    signal = [A B C D];
    cp_train = [signal(1,N-Ng+1:N) signal];
    trans_data= [predata cp_train suffixdata];
end
function [ M] = cazac(recv,N,Ng )
%UNTITLED3 此处显示有关此函数的摘要
%
 Ns = N + Ng;
 P1=zeros(1,2*Ns);
 R1=zeros(1,2*Ns); 
 P2=zeros(1,2*Ns);
 R2=zeros(1,2*Ns); 
 R=zeros(1,2*Ns); 
 P=zeros(1,2*Ns);
 stime = 0;

for d = Ns/2+1:1:2*Ns
    
    for m=0:N/4-1
       % tic;
        P1(d-Ns/2) = P1(d-Ns/2) + recv(d+m)*conj(recv(d+N/2-m-1));
        R1(d-Ns/2) = R1(d-Ns/2) + power(abs(recv(d+m)),2);
        P2(d-Ns/2) = P2(d-Ns/2) + recv(d+m+N/2)*recv(d+N/2+N/4+m);
        R2(d-Ns/2) = R2(d-Ns/2) + power(abs(recv(d+m+N/2)),2);
        stime=stime + toc;
    end
   % tic;
     P(d-Ns/2) = abs(P1(d-Ns/2))*abs(P2(d-Ns/2));
     R(d-Ns/2) =  0.5*abs(R1(d-Ns/2) +  R2(d-Ns/2))*abs(R2(d-Ns/2));
     stime=stime + toc;
end 
M=power(abs(P),2)./power((R),2);
%time = stime/(Ns*3/2);
end