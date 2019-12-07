%参考文献：
% @inproceedings{Wei2010Carrier,
%   title={Carrier Frequency Offset Estimation Using PN Sequence Iteration in OFDM Systems},
%   author={Wei, Jianghua and Yuan, Liu},
%   booktitle={Second International Conference on Networks Security},
%   year={2010},
% }
%函数功能：生成Jianghua Wei算法仿真数据
%输入参数 N:OFDM系统有效符号长度，即FFT电视
%         Ng:循环前缀的符号长度
function [recv_sig] = generate_RecvData_for_CFO_JianghuaWei(N,Ng,deltaf)
Ns=Ng+N;     %包括循环前缀的符号长度 
QAMTable =[7+7i,-7+7i,-7-7i,7-7i]; 
%=============生成数据符号=============
src = QAMTable(randi([0,3],N,1)+1); 
sym = ifft(src,N); 
cp_sym=[sym(1,N-Ng+1:N) sym];
%=============生成训练序列======================
buf=QAMTable(randi([0,3],N/4,1)+1); %1x(N/2)的矩阵
train1=zeros(1,N/2); 
index = 1; 
for n=1:2:N/2 
     train1(n)=buf(index); 
     index=index+1; 
end; 
sch = ifft(train1,N/2);   %[A A]的形式 
cp_train = [sch(1,N/4-Ng+1:N/4) sch sch];%[cp A A A A]格式
transmit_data = [cp_train cp_sym];
recv = transmit_data;
for k=1:length(recv)
    recv(k) = recv(k)*exp(sqrt(-1)*2*pi*deltaf*k/N);
end

recv_sig = recv;
end