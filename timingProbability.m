close all; 
clear all; 
clc; 
%参数定义 
N=256;       %FFT/IFFT 变换的点数或者子载波个数（Nu=N） 
Ng=N/8;      %循环前缀的长度 (保护间隔的长度) 
Ns=Ng+N;     %包括循环前缀的符号长度 
% step = 5;
% max_snr = 25;
%SNR= (-1)*max_snr:step:max_snr;
%SNR= 0:step:max_snr;

SNR = [-5,-4,-3,-2,-1,0,5,10,15,20,25];
%SNR = [-6,-4,-2,0,2,4,6,8,10,12,14,16,18,20,22,24,26];
loop = 1000;
%QAMTable=[1+1i,-1+1i,-1-1i,1-1i];
QAMTable=[7+7i,-7+7i,-7-7i,7-7i]; 
%%%-------------generate data symbol-------------
src = QAMTable(randi([0,3],N,1)+1); 
sym = ifft(src); 
cp_sym=[sym(1,N-Ng+1:N) sym];
predata = cp_sym;
suffixdata = cp_sym;
%%------------------------------------------------------
ren_mean_square = zeros(1,length(SNR));
fang_mean_square = zeros(1,length(SNR));
shao_mean_square = zeros(1,length(SNR));
minn_mean_square = zeros(1,length(SNR));
park_mean_square = zeros(1,length(SNR));
cazac_mean_square = zeros(1,length(SNR));
sch_mean_square = zeros(1,length(SNR));
sch_mean  = zeros(1,length(SNR));
ren_mean  = zeros(1,length(SNR));
minn_mean  = zeros(1,length(SNR));
park_mean  = zeros(1,length(SNR));
fang_mean  = zeros(1,length(SNR));
shao_mean  = zeros(1,length(SNR));
cazac_mean  = zeros(1,length(SNR));
%--------generate transmit symbol----------------
[transmit_data_schmidl,] = generate_trasmit_data(N,Ng,predata,suffixdata,'schmidl');
% transmit_data_minn = generate_trasmit_data(N,Ng,predata,suffixdata,'minn');
% transmit_data_park = generate_trasmit_data(N,Ng,predata,suffixdata,'park');
% transmit_data_cazac = generate_trasmit_data(N,Ng,predata,suffixdata,'cazac');
[transmit_data_minn,] = generate_trasmit_data(N,Ng,predata,suffixdata,'minn');
[transmit_data_park,]= generate_trasmit_data(N,Ng,predata,suffixdata,'park');
[transmit_data_ren,args_ren]= generate_trasmit_data(N,Ng,predata,suffixdata,'ren');
[transmit_data_fang,args_fang] = generate_trasmit_data(N,Ng,predata,suffixdata,'fang');
[transmit_data_shao,] = generate_trasmit_data(N,Ng,predata,suffixdata,'shao');
[transmit_data_liu,args_liu] = generate_trasmit_data(N,Ng,predata,suffixdata,'liubin');
[transmit_data_wang,args_wang] = generate_trasmit_data(N,Ng,predata,suffixdata,'wang');
[transmit_data_cazac,] = generate_trasmit_data(N,Ng,predata,suffixdata,'cazac');
%---------------------------------------------
for i=1:length(SNR)
    for j=1:loop
     [M_sch,] = schmidl(transmit_data_schmidl,N,Ng,SNR(i));
%     M_minn = minn(transmit_data_minn,N,Ng,SNR(i));
%     M_park = park(transmit_data_park,N,Ng,SNR(i));
%     M_cazacABCD=cazacABCD(transmit_data_cazac,N,Ng,SNR(i));
    [M_minn,time_minn] = minn(transmit_data_minn,N,Ng,SNR(i));
    [M_park, time_park] = park(transmit_data_park,N,Ng,SNR(i));
    [M_ren, time_ren] = ren(transmit_data_ren,N,Ng,SNR(i),args_ren);
    [M_fang, time_fang] = fang(transmit_data_fang,N,Ng,SNR(i),args_fang);
    [M_shao,time_shao] = shao(transmit_data_shao,N,Ng,SNR(i));
    [M_cazacABCD,time_cazac]=cazacABCD(transmit_data_cazac,N,Ng,SNR(i));
    
    [~, b] = max(M_sch) ;
    if b-177 == 0
        sch_mean_square(i) = sch_mean_square(i) + 1;
    end
    
    [~, b] = max(M_minn) ;
    if b-177 == 0
        minn_mean_square(i) = minn_mean_square(i) + 1;
    end
    
    [~, b] = max(M_park); 
    if b-177-N/2 == 0
        park_mean_square(i) = park_mean_square(i) + 1;
    end
    
    [~, b] = max(M_ren) ;
    if b-177 == 0
       ren_mean_square(i) = ren_mean_square(i) + 1;
    end
    
    [~, b] = max(M_fang) ;
    if b-177 == 0
        fang_mean_square(i) = fang_mean_square(i) + 1;
    end
    
    [~, b] = max(M_shao) ;
    if b-177 == 0
        shao_mean_square(i) = shao_mean_square(i) + 1;
    end
  
    [~, b] =max(M_cazacABCD); 
    if b-177 == 0
        cazac_mean_square(i) = cazac_mean_square(i) + 1;
    end
    end 
    ren_mean_square(i) = ren_mean_square(i)/loop;
    fang_mean_square(i) = fang_mean_square(i)/loop;
    shao_mean_square(i) = shao_mean_square(i)/loop;
    minn_mean_square(i) = minn_mean_square(i)/loop;
    park_mean_square(i) = park_mean_square(i)/loop;
    cazac_mean_square(i) = cazac_mean_square(i)/loop;
    sch_mean_square(i) = sch_mean_square(i)/loop;
end
figure(1);
d = 1:length(SNR);
plot(SNR, sch_mean_square(d) ,'LineWidth',1);
hold on
plot(SNR, minn_mean_square(d) ,'-dg','LineWidth',1);
hold on;
plot(SNR,park_mean_square(d) ,'-o','LineWidth',1);
hold on;
plot(SNR,ren_mean_square(d) ,'-xb','LineWidth',1);
hold on;
plot(SNR,fang_mean_square(d) ,'-s','LineWidth',1);
hold on;
plot(SNR,shao_mean_square(d) ,'-^k','LineWidth',1);
hold on;
plot(SNR, cazac_mean_square(d) ,'-vr','LineWidth',1);
hold on;
legend('schmidl''s method','minn''s method','park''s method','ren''s method','fang''s method','shao''s method','proposed''s method');
xlabel('SNR(dB)'); 
ylabel('Timing extimation hit Rate '); 