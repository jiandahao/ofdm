%@INPROCEEDINGS{8049990, 
% author={S. Bhowmick and K. Vasudevan}, 
% booktitle={2017 4th International Conference on Signal Processing and Integrated Networks (SPIN)}, 
% title={Matched filter based integer frequency offset estimation in OFDM systems}, 
% year={2017}, 
% volume={}, 
% number={}, 
% pages={442-447}, 
% keywords={AWGN;matched filters;OFDM modulation;probability;Rayleigh channels;matched filter based integer frequency offset estimation;OFDM systems;carrier frequency offset;integer frequency offset;IFO;CFO;multipath Rayleigh fading channel;additive white Gaussian noise;AWGN;fractional frequency offset;FFO;improved probability;preamble aided integer frequency synchronization scheme;OFDM;Frequency-domain analysis;Time-domain analysis;Estimation;Frequency estimation;Fading channels;Convolution;Frequency synchronization;OFDM;Preamble;Matched filter}, 
% doi={10.1109/SPIN.2017.8049990}, 
% ISSN={}, 
% month={Feb},}
function [failure_rate] =IFO_Sch(config)
    %==========获取OFDM系统仿真配置===========
    %config = OFDMSystemConfig;
    %ffo_df = config.ffo_df;%小数倍频偏
    %ifo_df = config.ifo_df;%整数倍频偏
    simu_times = config.simu_times;%仿真次数
    N = config.N;%OFDM符号周期长度
    Ng = config.Ng;%循环前缀长度
    SNR = config.SNR;%信噪比
    
    [transmit_data,v] = generate_data(N,Ng);
    
    success_rate = zeros(1,length(SNR));
    failure_rate = zeros(1,length(SNR));
    %===========进行蒙特卡洛仿真============
    for n=1:length(SNR)  
        cnt = 0;
        for s = 1:simu_times
            ffo_df = -0.5;
            ifo_df = randi([0,N/4])*2;
            df = ifo_df + ffo_df;
            %[transmit_data,v] = generate_data(N,Ng);
            recv_sig = add_CFO(transmit_data,df,N);
            recv_sig_add_noise = awgn(recv_sig,SNR(n),'measured');
            ifo_temp =  get_IFO(recv_sig_add_noise,N,Ng,v);%频域估计算法
            if (ifo_temp == ifo_df ) 
                cnt = cnt + 1;
            end
        end
        success_rate(n) = cnt/simu_times;
        failure_rate(n) = 1 - success_rate(n);
    end
    figure
    d = 1:length(SNR);
    plot(SNR,success_rate(d));
    %semilogy(SNR,1-success_rate(d));
    xlabel('SNR(dB)'); 
    ylabel('Probabilty of Success');  
    legend('schmidl算法');
    grid on
end

function[signal,v] = generate_data(N,Ng)
%============构建数据=================
    QAMTable =[7+7i,-7+7i,-7-7i,7-7i]; 
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
    transmit_data = [cp_train1 cp_train2];
    signal = transmit_data;
end

function[recv_sig] = add_CFO(trans_sig,df,N)
%===========添加频偏===========
recv_sig = zeros(1,length(trans_sig));
for k=1:length(trans_sig)
    %添加频偏
    recv_sig(k) = trans_sig(k)*exp(sqrt(-1)*2*pi*df*k/N);
end  
end

function [IFO_df] =  get_IFO(recv_sig,N,Ng,v)
recv = recv_sig;
Ns=N+Ng;
r1 = recv(1,Ng+1:Ns);%获取第一个OFDM符号
r2 = recv(1,N+2*Ng+1:2*Ns);%获取第二个OFDM符号
R1 = fft(r1);
R2 = fft(r2);

for g = -N/4:N/4
    for k = 1:2:N
        if k+2*g <=0
             B1(k) = conj(R1(k+2*g+N))*conj(v(k))*R2(k+2*g+N);
        end
        if (k+2*g<=N) &&  (k+2*g>0)
            B1(k) = conj(R1(k+2*g))*conj(v(k))*R2(k+2*g);
        end
        if k+2*g>N
            B1(k) =conj(R1(k+2*g-N))*conj(v(k))*R2(k+2*g-N);
        end
        B2(k) = abs(R2(k))^2;
    end
    B(g+N/4+1) = (abs(sum(B1)))^2/(2*sum(B2)^2);
end
[~,IFO_df] = max(B);
% figure
% plot(-N/4:N/4,B);
IFO_df = 2*(IFO_df -(N/4+1));
%CFO_df = 2*(IFO_df-11) + FFO_df
end