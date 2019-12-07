close all; 
clear all; 
clc; 
%参数定义 
N=256;       %FFT/IFFT 变换的点数或者子载波个数（Nu=N） 
Ng=N/8;      %循环前缀的长度 (保护间隔的长度) 
Ns=Ng+N;     %包括循环前缀的符号长度 
SNR=15;
QAMTable =[7+7i,-7+7i,-7-7i,7-7i]; 
%QAMTable=[1+1i,-1+1i,-1-1i,1-1i]; 
src = QAMTable(randi([0,3],N,1)+1); 
sym = ifft(src); 
cp_sym=[sym(1,N-Ng+1:N) sym];
predata = cp_sym;
suffixdata = cp_sym;
transmit_data_schmidl = generate_trasmit_data(N,Ng,predata,suffixdata,'schmidl');
[transmit_data_minn,] = generate_trasmit_data(N,Ng,predata,suffixdata,'minn');
[transmit_data_park,]= generate_trasmit_data(N,Ng,predata,suffixdata,'park');
[transmit_data_ren,args_ren]= generate_trasmit_data(N,Ng,predata,suffixdata,'ren');
[transmit_data_fang,args_fang] = generate_trasmit_data(N,Ng,predata,suffixdata,'fang');
[transmit_data_shao,] = generate_trasmit_data(N,Ng,predata,suffixdata,'shao');
[transmit_data_liu,args_liu] = generate_trasmit_data(N,Ng,predata,suffixdata,'liubin');
[transmit_data_wang,args_wang] = generate_trasmit_data(N,Ng,predata,suffixdata,'wang');
[transmit_data_cazac,] = generate_trasmit_data(N,Ng,predata,suffixdata,'cazac');
Mtime=zeros(6,1);
%=============偏差概率分布相关==============
record_sch = zeros(1,2*N);
record_minn = zeros(1,2*N);
record_park = zeros(1,2*N);
record_ren = zeros(1,2*N);
record_fang = zeros(1,2*N);
record_shao = zeros(1,2*N);
record_cazac = zeros(1,2*N);
%========================================
simulation_times = 1;
for i = 1:simulation_times
M_sch = schmidl(transmit_data_schmidl,N,Ng,SNR);

%counting_time = 1;
 %for i=1:counting_time
[M_minn,time_minn] = minn(transmit_data_minn,N,Ng,SNR);
[M_park, time_park] = park(transmit_data_park,N,Ng,SNR);
[M_ren, time_ren] = ren(transmit_data_ren,N,Ng,SNR,args_ren);
[M_fang, time_fang] = fang(transmit_data_fang,N,Ng,SNR,args_fang);
[M_shao,time_shao] = shao(transmit_data_shao,N,Ng,SNR);
%$[M_liu,time_liu] = liubin(transmit_data_liu,N,Ng,SNR,args_liu);
%[M_wang,time_wang] = wang(transmit_data_wang,N,Ng,SNR,args_wang);
[M_cazacABCD,time_cazac]=cazacABCD(transmit_data_cazac,N,Ng,SNR);
% [~,b]=max(M_sch);
% record_sch(b) = record_sch(b) + 1;
% [~,b]=max(M_minn);
% record_minn(b) = record_minn(b) + 1;
% [~,b]=max(M_park);
% record_park(b) = record_park(b) + 1;
% [~,b]=max(M_ren);
% record_ren(b) = record_ren(b) + 1;
% [~,b]=max(M_fang);
% record_fang(b) = record_fang(b) + 1;
% [~,b]=max(M_shao);
% record_shao(b) = record_shao(b) + 1;
% [~,b]=max(M_cazacABCD);
% record_cazac(b) = record_cazac(b) + 1;
% Mtime(1)=Mtime(1)+time_minn;
% Mtime(2)=Mtime(2)+time_park;
% Mtime(3)=Mtime(3)+time_ren;
% Mtime(4)=Mtime(4)+time_fang;
% Mtime(5)=Mtime(5)+time_shao;
% Mtime(6)=Mtime(6)+time_cazac;
end
 % vpa(Mtime/counting_time,8)
%  draw('schmidl',M_sch,1);
%  draw('minn',M_minn,2);
%  draw('cazacABCD',M_cazacABCD,7);
%  draw('park',M_park,3);
% figure(1);
% d=1:1:350;
% % plot(d,M_sch(d),':');
% % 
% % hold on;
% % plot(d,M_minn(d),'--');
% % hold on;
% plot(d,M_park(d+N/2),'-^');
% hold on;
% plot(d,M_cazacABCD(d),'-x');
% hold on;
% legend('schmidl''s method','minn''s method','park''s method','proposed''s method');
% xlabel('Time (sample)'); 
% ylabel('Timing Metric'); 
% axis([0,350,0,1.1]); 
%================绘制偏差概率分布图============
% d=-176:2*N-177;
% subplot(7,1,1);
% plot(d,record_sch/simulation_times);
% subplot(7,1,2);
% plot(d,record_minn/simulation_times);
% subplot(7,1,3);
% plot(d,record_park/simulation_times);
% subplot(7,1,4);
% plot(d,record_ren/simulation_times);
% subplot(7,1,5);
% plot(d,record_fang/simulation_times);
% subplot(7,1,6);
% plot(d,record_shao/simulation_times);
% subplot(7,1,7);
% plot(d,record_cazac/simulation_times);
%=============================================
%================Draw Metric function curves=============
  draw('schmidl',M_sch,5);
  draw('Minn',M_minn,8);
  draw('Park',M_park,1);
  draw('Ren',M_ren,2);
  draw('Fang',M_fang,3);
  draw('Shao',M_shao,4);
% %  %draw('Liubin',M_liu,5);
%  %draw('wang',M_wang,6);
  draw('Proposed',M_cazacABCD,7);
%=====================================
time_minn
time_park
time_ren
time_fang
time_shao
time_cazac
% grid on; 
