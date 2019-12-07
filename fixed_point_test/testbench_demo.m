close all; 
clear all; 
clc; 

N=256;      
Ng=N/8;    
Ns=Ng+N;    
SNR=15;
QAMTable=[7+7i,-7+7i,-7-7i,7-7i]; 
%%%-------------generate data symbol-------------
src = QAMTable(randi([0,3],N,1)+1); 
sym = ifft(src); 
%sym = src;
cp_sym=[sym(1,N-Ng+1:N) sym];
predata = cp_sym;
suffixdata = cp_sym;
[transmit_data_cazac,] = generate_trasmit_data(N,Ng,predata,suffixdata);
% [M_cazacABCD,time_cazac]=cazacABCD(transmit_data_cazac,N,Ng,1000);
% draw('Proposed',M_cazacABCD,1);
% if SNR<100
 transmit_data_cazac = awgn(transmit_data_cazac,SNR);
% end
M = zeros(1,2*Ns);
for d = Ns/2+1:1:2*Ns
    
    [ M(d-Ns/2) ] = cazacABCDForC(transmit_data_cazac(1,d:d+N-1));
end
figure(1);
plot(real(transmit_data_cazac));
figure(2);
plot(imag(transmit_data_cazac));
figure(3);
d=1:1:400;
plot(d,M(d));