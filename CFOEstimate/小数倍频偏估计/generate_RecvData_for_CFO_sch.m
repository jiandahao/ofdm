%函数功能：生成schmidl算法仿真使用的加入频偏的传输数据（recv_sig）,以及差分序列（v_sig）
%输入参数 N:OFDM系统有效符号长度，即FFT电视
%         Ng:循环前缀的符号长度
function [recv_sig,v_sig] = generate_RecvData_for_CFO_sch(N,Ng,deltaf)
Ns=Ng+N;     %包括循环前缀的符号长度 
QAMTable =[7+7i,-7+7i,-7-7i,7-7i]; 
%=============生成数据符号=============
src = QAMTable(randi([0,3],N,1)+1); 
sym = ifft(src,N); 
cp_sym=[sym(1,N-Ng+1:N) sym];
predata = cp_sym;
suffixdata = cp_sym;
%=============生成训练序列1======================
buf=QAMTable(randi([0,3],N/2,1)+1); %1x(N/2)的矩阵
train1=zeros(1,N); 
index = 1; 
for n=1:2:N 
     train1(n)=buf(index); 
     index=index+1; 
end; 
sch = ifft(train1,N);   %[A A]的形式 
cp_train1 = [sch(1,N-Ng+1:N) sch];
%=============生成训练序列2======================
buf=QAMTable(randi([0,3],N/2,1)+1); %1x(N/2)的矩阵
train2=zeros(1,N); 
index = 1; 
for n=1:2:N 
     train2(n)=buf(index); 
     index=index+1; 
end; 
sch = ifft(train2,N);   %[A A]的形式 
cp_train2 = [sch(1,N-Ng+1:N) sch];
%==============计算差分序列====================
v =sqrt(2)*train2./train1;
%=============经过信道=========================
transmit_data = [cp_train1 cp_train2 suffixdata]; 
recv = transmit_data;
for k=1:length(recv)
    recv(k) = recv(k)*exp(i*2*pi*deltaf*k/N);
end

recv_sig = recv;
v_sig = v;
end