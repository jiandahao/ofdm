function [ M,time ] = cazacABCD( transmit_data,N,Ng,SNR )
%UNTITLED3 此处显示有关此函数的摘要
%
Ns=Ng+N;     %包括循环前缀的符号长度 
if SNR<100
recv = awgn(transmit_data,SNR);
else
    recv = transmit_data;
end
 P1=zeros(1,2*Ns);
 R1=zeros(1,2*Ns); 
 P2=zeros(1,2*Ns);
 R2=zeros(1,2*Ns); 
 R=zeros(1,2*Ns); 
 P=zeros(1,2*Ns);
 stime = 0;
% tic;

% %%===================method first===============================
% for d = Ns/2+1:1:2*Ns
%     
%     for m=0:N/4-1
%         tic;
% 
%         P1(d-Ns/2) = P1(d-Ns/2) + recv(d+m)*recv(d+N/2-m-1);
%         R1(d-Ns/2) = R1(d-Ns/2) + power(abs(recv(d+m+N/4)),2);
%         P2(d-Ns/2) = P2(d-Ns/2) + recv(d+m+N/2)*recv(d+N/2+N/4+m);
%         R2(d-Ns/2) = R2(d-Ns/2) + power(abs(recv(d+m+N/2)),2);
%         stime=stime + toc;
%     end
% end 
% M=(power(abs(P1),2)./power(abs(R1),2)).*(power(abs(P2),2)./power(abs(R2),2));
% %===============================================================



% %%===================method second===============================
% for d = Ns/2+1:1:2*Ns
%     
%     for m=0:N/4-1
%         tic;
%         P1(d-Ns/2) = P1(d-Ns/2) + recv(d+m)*recv(d+N/2-m-1);
%         P2(d-Ns/2) = P2(d-Ns/2) + recv(d+m+N/2)*recv(d+N/2+N/4+m);
%         R2(d-Ns/2) = R2(d-Ns/2) + power(abs(recv(d+m+N/2)),2);
%         stime=stime + toc;
%     end
%     for m=0:N/2-1
%         R1(d-Ns/2) = R1(d-Ns/2) + power(abs(recv(d+m)),2);
%      end  
% end 
% M=(power(abs(P1),2)./power(abs(0.5*(R1)),2)).*(power(abs(P2),2)./power(abs(R2),2));
% %===============================================================


% %%===================method third===============================
% for d = Ns/2+1:1:2*Ns
%     
%     for m=0:N/4-1
%         tic;
%         P1(d-Ns/2) = P1(d-Ns/2) + recv(d+m)*recv(d+N/2-m-1);
%         P2(d-Ns/2) = P2(d-Ns/2) + recv(d+m+N/2)*recv(d+N/2+N/4+m);
%         stime=stime + toc;
%     end
%      for m=0:N-1
%          R(d-Ns/2) = R(d-Ns/2) + power(abs(recv(d+m)),2);
%       end  
% end 
%  R1=0.25*R;
%  R2=0.25*R;
% M=(power(abs(P1),2)./power(abs(R1),2)).*(power(abs(P2),2)./power(abs(R2),2));
% %===============================================================



%%===================method forth===============================
for d = Ns/2+1:1:2*Ns
    
    for m=0:N/4-1
        tic;
        P1(d-Ns/2) = P1(d-Ns/2) + recv(d+m)*recv(d+N/2-m-1);
        R1(d-Ns/2) = R1(d-Ns/2) + power(abs(recv(d+m)),2);
        %P2(d-Ns/2) = P2(d-Ns/2) + recv(d+m+N/2)*recv(d+N/2+N/4+m);
        P2(d-Ns/2) = P2(d-Ns/2) + recv(d+m+N/2)*conj(recv(d+N/2+N/4+m));
        R2(d-Ns/2) = R2(d-Ns/2) + power(abs(recv(d+m+N/2)),2);
        stime=stime + toc;
    end
    tic;
     P(d-Ns/2) = abs(P1(d-Ns/2))*abs(P2(d-Ns/2));
     R(d-Ns/2) =  0.5*abs(R1(d-Ns/2) +  R2(d-Ns/2))*abs(R2(d-Ns/2));
     stime=stime + toc;
end 

%M=(power(abs(P1),2)./power(abs(0.5*(R1+R2)),2)).*(power(abs(P2),2)./power(abs(R2),2));
M=power(abs(P),2)./power((R),2); 



time = stime/(Ns*3/2);
end

