function [ M,time ] = shao(transmit_data,N,Ng,SNR )
Ns = N + Ng;
if SNR<100
recv = awgn(transmit_data,SNR);
else
    recv = transmit_data;
end
%*****************¼ÆËã·ûºÅ¶¨Ê±***************************** 
P=zeros(1,2*Ns);
 R=zeros(1,2*Ns); 
 v=zeros(1,N/4);
for k=0:N/4-1
 v(k+1)= exp((-1)*sqrt(-1)*2*pi*k*k/N);
end
stime =0;
% tic;
for d = Ns/2+1:1:2*Ns
%     tic;
    for m=0:N/4-1
        tic;
        P(d-Ns/2) = P(d-Ns/2) + conj(recv(d+m))*recv(d+N/2+m)*conj(v(m+1))+conj(recv(d+N/4+m))*recv(d+3*N/4+m)*v(N/4-m);
        stime = stime + toc;
    end
    for m=0:N-1
        tic;
        R(d-Ns/2) = R(d-Ns/2) + power(abs(recv(d+m)),2);
        stime = stime + toc;
    end
%     stime = stime + toc;
end        
M=power(abs(P),2)./power((0.5*R),2); 
% time=toc;
time = stime/(3*Ns/2);
end

