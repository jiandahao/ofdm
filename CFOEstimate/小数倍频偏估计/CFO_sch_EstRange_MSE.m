% close all; 
clear all; 
clc; 
%============获取OFDM系统配置==========
%[N,Ng,ffo_df,ifo_df,SNR,QAMTable,simu_times] = get_OFDMSystem_config();
config = OFDMSystemConfig;
N = config.N;
Ng = config.Ng;
ffo_df = config.ffo_df;
ifo_df = config.ifo_df;
SNR = config.SNR;
simu_times = config.simu_times;
%=============生成数据符号=============
QAMTable =[7+7i,-7+7i,-7-7i,7-7i]; 
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
v_sig =sqrt(2)*train2./train1;
%=============构建传输信号====================
transmit_data = [cp_train1 cp_train2 suffixdata]; 
%==============计算无噪声情况下频偏估计的情况================
deltaf = -3:0.05:3;
CFO_sch_ffo = zeros(1,length(deltaf));
CFO_sch_ifo = zeros(1,length(deltaf));
CFO_sch_cfo = zeros(1,length(deltaf));
recv_sig = zeros(1,length(transmit_data));
for n=1:length(deltaf)  
    for k=1:length(transmit_data)
        %添加频偏
        recv_sig(k) = transmit_data(k)*exp(sqrt(-1)*2*pi*deltaf(n)*k/N);
    end
    [CFO_sch_ffo(n),CFO_sch_ifo(n),CFO_sch_cfo(n)] =  CFO_sch(N,Ng,recv_sig,v_sig);
end
%==============绘制不同频偏值下算法频偏估计值曲线（无噪声）====
figure
plot(deltaf,CFO_sch_ffo,'-','LineWidth',1);
hold on
plot(deltaf,CFO_sch_ifo,'-.','LineWidth',1);
hold on
plot(deltaf,CFO_sch_cfo,'--','LineWidth',1);
grid on
set(gca,'XTick',-3:0.5:3);
xlabel('归一化频偏'); 
ylabel('Schmidl&Cox算法频偏估计值'); 
legend('小数倍频偏估计','整数倍频偏估计','总频偏估计');
% figure(2)
% plot(deltaf,CFO_sch_ifo,'LineWidth',1);
% grid on
% set(gca,'XTick',-3:0.5:3);
% xlabel('归一化频偏'); 
% ylabel('整数倍频偏估计均值'); 
% 
% figure(3)
% plot(deltaf,CFO_sch_cfo,'LineWidth',1);
% grid on
% set(gca,'XTick',-3:0.5:3);
% xlabel('归一化频偏'); 
% ylabel('频偏估计均值'); 
%================计算不同信噪比环境下，频偏估计均方误差（deltaf=0.3）=====================
% FFO_sch_df_mse = zeros(1,length(SNR));%Mean square error
% FFO_sch_df_mfo = zeros(1,length(SNR));%Mean frequency offset
CFO_sch_df_mse = zeros(1,length(SNR));%Mean square error
CFO_sch_df_mfo = zeros(1,length(SNR));%Mean frequency offset
%simu_times = 1000;
% ffo_df = 0.5;%小数倍频偏,范围[-1,1]
% ifo_df = 3;%整数倍频偏
df =ffo_df + ifo_df;%总频偏
recv_sig = zeros(1,length(transmit_data));
sch_max_offset = zeros(1,length(SNR));%记录最大偏差率。
for k=1:length(transmit_data)
    %添加频偏
    recv_sig(k) = transmit_data(k)*exp(sqrt(-1)*2*pi*df*k/N);
end    
for n=1:length(SNR)  
    for s = 1:simu_times
        recv_sig_add_noise = awgn(recv_sig,SNR(n));
        %[FFO_df,IFO_df,CFO_df]
        %[ffo_temp,~,cfo_temp] = CFO_sch(N,Ng,recv_sig_add_noise,v_sig);
        [cfo_temp,~,~] = CFO_sch(N,Ng,recv_sig_add_noise,v_sig);
%         FFO_sch_df_mse(n) = FFO_sch_df_mse(n)+(ffo_temp - ffo_df)^2;
%         FFO_sch_df_mfo(n) = FFO_sch_df_mfo(n)+ ffo_temp;
        CFO_sch_df_mse(n) = CFO_sch_df_mse(n)+(cfo_temp - df)^2;
        CFO_sch_df_mfo(n) = CFO_sch_df_mfo(n)+ cfo_temp;
        sch_max_offset(n) = max(sch_max_offset(n),abs(cfo_temp-df));
    end
%     FFO_sch_df_mse(n) = FFO_sch_df_mse(n) / simu_times;
%     FFO_sch_df_mfo(n) = FFO_sch_df_mfo(n) / simu_times;
    CFO_sch_df_mse(n) = CFO_sch_df_mse(n) / simu_times;
    CFO_sch_df_mfo(n) = CFO_sch_df_mfo(n) / simu_times;
end

%==========================绘制MSE曲线=========================
% figure
% d = 1:length(SNR);
% semilogy(SNR,FFO_sch_df_mse(d));
% grid on
% set(gca,'XTick',0:5:25);
% xlabel('SNR(dB)'); 
% ylabel('Schmidl&Cox算法频偏估计均方误差（MSE）'); 
%====================将MSE结果记录在文件中================
save('CFO_sch_df_mse.mat','CFO_sch_df_mse','CFO_sch_df_mfo','sch_max_offset');