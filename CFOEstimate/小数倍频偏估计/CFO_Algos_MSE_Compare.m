close all
clear all
clc
%使用之前可以先通过类OFDMSystemConfig.m或函数get_OFDMSystem_config.m进行OFDM系统设置
CFO_proposed_EstRange_MSE;
CFO_sch_EstRange_MSE;
CFO_moose_EstRange_MSE;
CFO_CP_EstRange_MSE;
CFO_Jiang_EstRange_MSE;
clear all
clc

load CFO_Jiang_df_mse.mat
load CFO_proposed_df_mse.mat;
load CFO_sch_df_mse.mat;
load CFO_moose_df_mse.mat;
load CFO_CP_df_mse.mat;
label = ['-oc';'-db';'-^k';'-vr'];
d = 1:length(SNR);
%=======计算误比特率=============
% D_proposed = zeros(1,length(SNR));
% D_sch = zeros(1,length(SNR));
% D_moose = zeros(1,length(SNR));
% D_CP = zeros(1,length(SNR));
% C = 10*pi/(3*256*log(10));
% for k=1:length(SNR)
%    D_proposed(k) =  C*SNR(k)*(0.3-CFO_proposed_df_mfo(k))^2;
%    D_sch(k) =  C*SNR(k)*(0.3-CFO_sch_df_mfo(k))^2;
%    D_CP(k) =  C*SNR(k)*(0.3-CFO_CP_df_mfo(k))^2;
%    D_moose(k) =  C*SNR(k)*(0.3-CFO_moose_df_mfo(k))^2;
%     
% end

%==========绘制误比特率曲线===========

% figure 
% semilogy(SNR,D_CP(d),'-oc','LineWidth',1);
% hold on
% semilogy(SNR,D_sch(d),'-db','LineWidth',1);
% hold on
% semilogy(SNR,D_moose(d),'-^k','LineWidth',1);
% hold on
% semilogy(SNR,D_proposed(d),'-vr','LineWidth',1);
% hold on
% grid on
% set(gca,'XTick',0:5:25);
% xlabel('SNR(dB)'); 
% ylabel('误比特率D(dB)'); 
% legend('Van de Beek算法(Ng=32,N=256)','Schmidl&Cox算法','Moose算法','本文算法');
%==========绘制频偏估计均方误差曲线（MSE）========
figure
semilogy(SNR,CFO_CP_df_mse(3,d),'-oc','LineWidth',1);
hold on
semilogy(SNR,CFO_sch_df_mse(d),'-db','LineWidth',1);
hold on
semilogy(SNR,CFO_moose_df_mse(d),'-^k','LineWidth',1);
hold on
semilogy(SNR,CFO_proposed_df_mse(d),'-vr','LineWidth',1);
hold on
grid on
set(gca,'XTick',0:5:25);
xlabel('SNR(dB)'); 
ylabel('频偏估计均方误差（MSE）'); 
legend('Van de Beek算法(Ng=32,N=256)','Schmidl&Cox算法','Moose算法','本文算法');
%===========绘制频偏估计均值曲线==============
figure
plot(SNR,CFO_CP_df_mfo(3,d),'-oc','LineWidth',1);
hold on
plot(SNR,CFO_sch_df_mfo(d),'-db','LineWidth',1);
hold on
plot(SNR,CFO_moose_df_mfo(d),'-^k','LineWidth',1);
hold on
plot(SNR,CFO_proposed_df_mfo(d),'-vr','LineWidth',1);
hold on
grid on
set(gca,'XTick',0:5:25);
xlabel('SNR(dB)'); 
ylabel('小数倍频偏估计均值'); 
legend('Van de Beek算法(Ng=32,N=256)','Schmidl&Cox算法','Moose算法','本文算法');
%===========绘制总频偏估计均值曲线==============
figure
plot(SNR,CFO_sch_df_mfo(d),'-^k','LineWidth',1);
hold on
plot(SNR,CFO_proposed_df_mfo(d),'-vr','LineWidth',1);
hold on
grid on
set(gca,'XTick',0:5:25);
xlabel('SNR(dB)'); 
ylabel('频偏估计均值'); 
legend('Schmidl&Cox算法','本文算法');
 %==========绘制总频偏估计均方误差曲线（MSE）========
 %Schmidl
figure
semilogy(SNR,CFO_sch_df_mse(d),'-^k','LineWidth',1);
hold on
semilogy(SNR,CFO_proposed_df_mse(d),'-vr','LineWidth',1);
grid on
set(gca,'XTick',0:5:25);
xlabel('SNR(dB)'); 
ylabel('频偏估计均方误差（MSE）'); 
legend('Schmidl&Cox算法','本文算法');