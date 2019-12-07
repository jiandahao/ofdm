close all; 
clear all; 
clc; 
deltaf = 0.3;
simu_times = 5000;%仿真次数
SNR = [0,5,10,15,20,25];
%SNR = 25;
N = 256;
Ng=N/8;
moose_data = generate_RecvData_for_CFO_moose(N,Ng,deltaf);
CP_data = generate_RecvData_for_CFO_CP(N,Ng,deltaf);
[Sch_data,sch_v_sig] = generate_RecvData_for_CFO_sch(N,Ng,deltaf);
[proposed_data,Ts_sym] =  generate_RecvData_for_CFO_proposed(N,Ng,deltaf);

CFO_CP_df = zeros(1,length(SNR));
CFO_moose_df = zeros(1,length(SNR));
CFO_sch_df = zeros(1,length(SNR));
CFO_proposed_df = zeros(1,length(SNR));
for n = 1:length(SNR)
    for k = 1:simu_times
       CP_data_add_noise = awgn(CP_data,SNR(n));
       CFO_CP_df(n) = CFO_CP_df(n) + CFO_CP(N,Ng,CP_data_add_noise);
 
       moose_data_add_noise = awgn(moose_data,SNR(n));
       CFO_moose_df(n) = CFO_moose_df(n) + CFO_moose(N,Ng,moose_data_add_noise);
       
       Sch_data_add_noise = awgn(Sch_data,SNR(n));
       CFO_sch_df(n) =  CFO_sch_df(n) + CFO_sch(N,Ng,Sch_data_add_noise,sch_v_sig);
       
        proposed_data_add_noise = awgn(proposed_data,SNR(n));
        CFO_proposed_df(n) = CFO_proposed_df(n) + CFO_proposed(N,Ng,proposed_data_add_noise,Ts_sym);
    end
    CFO_CP_df(n) = CFO_CP_df(n)/simu_times;
    CFO_moose_df(n) =  CFO_moose_df(n)/simu_times;
    CFO_sch_df(n) =  CFO_sch_df(n)/simu_times;
    CFO_proposed_df(n) = CFO_proposed_df(n)/simu_times;
end
% moose_df = CFO_moose(SNR,deltaf,simu_times);
% [sch_FFO_df,sch_IFO_df,sch_CFO_df] = CFO_sch(SNR,deltaf,simu_times);
% [CP_FFO_df,] = CFO_CP(SNR,deltaf,simu_times,256,32);
% [p_FFO_df,p_IFO_df,p_CFO_df] = CFO_proposed(SNR,deltaf,simu_times);
d = 1:1:length(SNR);
figure(1);
plot(SNR, CFO_moose_df(d),'-*m','LineWidth',1);
hold on;
plot(SNR, CFO_sch_df(d),'-dg','LineWidth',1);
hold on;
plot(SNR,CFO_CP_df(d),'-o','LineWidth',1);
hold on;
plot(SNR,CFO_proposed_df(d),'-xb','LineWidth',1);
hold on;
legend('moose算法','Schmidl&Cox算法','基于循环前缀','改进算法');
xlabel('SNR(dB)'); 
ylabel('小数倍频偏估计值'); 