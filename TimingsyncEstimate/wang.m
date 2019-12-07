function [ M,time ] = wang( transmit_data,N,Ng,SNR,args )
%UNTITLED3 此处显示有关此函数的摘要
%
Ns=Ng+N;     %包括循环前缀的符号长度 
if SNR<100
recv = awgn(transmit_data,SNR);
else
    recv = transmit_data;
end
 P=zeros(1,2*Ns);
 R=zeros(1,2*Ns); 
 pn = args;
 stime = 0;
% tic;
for d = Ns/2+1:1:2*Ns
    tic;
    for m=0:N/4-1
        P(d-Ns/2) = P(d-Ns/2) + conj(recv(d+m))*recv(d+3*N/4+m)*pn(m+1);
    end
    for i =0:3
        for m=0:N/4-1
            R(d-Ns/2) = R(d-Ns/2) + power(abs(recv(d+m+i*N/4)),2);
        end
    end
    stime = stime + toc;
end         
M=power(abs(P),2)./(0.25*power(abs(R),2)); 
%time=toc;
time = stime/(3*Ns/2);
end

