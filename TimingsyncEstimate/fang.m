function [ M,time ] = fang(transmit_data,N,Ng,SNR,args )
Ns = N + Ng;
if SNR<100
recv = awgn(transmit_data,SNR);
else
    recv = transmit_data;
end
%*****************¼ÆËã·ûºÅ¶¨Ê±***************************** 
 P=zeros(1,2*Ns);
 R=zeros(1,2*Ns);
 v=args;
stime=0;
for d = Ns/2+1:1:2*Ns
%      tic;
    for m=0:N/2-1
        tic;
        P(d-Ns/2) = P(d-Ns/2) + recv(d+m)*conj(recv(d+N/2+m))*conj(v(m+1));
        stime = stime + toc;
    end
    for m =0:N-1
        tic;
        R(d-Ns/2) = R(d-Ns/2) + power(abs(recv(d+m)),2);
        stime = stime + toc;
    end
%     stime = stime + toc;
end         
M=power(abs(P),2)./(0.5*power(abs(R),2)); 
time = stime/(3*Ns/2);

end

