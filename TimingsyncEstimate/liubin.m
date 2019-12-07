function [ M,time] = liubin(transmit_data,N,Ng,SNR,args )
Ns = N + Ng;
if SNR<100
recv = awgn(transmit_data,SNR);
else
    recv = transmit_data;
end
%*****************¼ÆËã·ûºÅ¶¨Ê±***************************** 
P1=zeros(1,2*Ns);
P2=zeros(1,2*Ns);
R=zeros(1,2*Ns); 
M=zeros(1,2*Ns);
pn=args;
stime = 0;
%tic;
for d = Ns/2+1:1:2*Ns
    tic;
    for m=0:N/2-1
        P1(d-Ns/2) = P1(d-Ns/2) + conj(recv(d+m))*recv(d+N/2+m)*pn(m+1);
        P2(d-Ns/2) = P2(d-Ns/2) + recv(d+m)*recv(d+N-m-1)*pn(m+1);
    end
    for m =0:N-1
        R(d-Ns/2) = R(d-Ns/2) + power(abs(recv(d+m)),2);
    end
%     M(d-Ns/2)=(abs(P1(d-Ns/2)).*abs( P2(d-Ns/2)))./(0.5*power(abs(R(d-Ns/2),2))); 
    stime = stime + toc;
end         
M=(abs(P1).*abs(P2))./(0.5*power(abs(R),2)); 
time = stime/(3*Ns/2);
end

