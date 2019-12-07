function [failure_rate] =IFO_Shao(config)
    %==========获取OFDM系统仿真配置===========
    %config = OFDMSystemConfig;
    %ffo_df = config.ffo_df;%小数倍频偏
    %ifo_df = config.ifo_df;%整数倍频偏
    simu_times = config.simu_times;%仿真次数
    N = config.N;%OFDM符号周期长度
    Ng = config.Ng;%循环前缀长度
    SNR = config.SNR;%信噪比
    %============构建训练符号=================
    a = zeros(1,N/4);
    for n=0:N/4-1
        a(n+1) = exp(sqrt(-1)*pi*n*n/N);
    end
    A=a;
    B = conj(a(1,N/4:-1:1));
    train_symbol = [A B conj(A) conj(B)];
    cp_train = [train_symbol(1,N-Ng+1:N) train_symbol];
    %============添加传输数据===================
    QAMTable=[7+7i,-7+7i,-7-7i,7-7i]; 
    src = QAMTable(randi([0,3],N,1)+1); 
    sym = ifft(src); 
    cp_sym=[sym(1,N-Ng+1:N) sym];
    predata = cp_sym;
    suffixdata = cp_sym;
    transmit_data = [predata cp_train suffixdata];

    success_rate = zeros(1,length(SNR));
    failure_rate = zeros(1,length(SNR));
    for n=1:length(SNR) 
        cnt = 0;
        for s = 1:simu_times
            %============添加频偏====================
            ffo_df = -0.5 + rand(1,simu_times);
            ifo_df = 16;
            df =ffo_df + ifo_df;%总频偏
            recv_sig = add_CFO(transmit_data,df(s),N);
            recv_sig_add_noise = awgn(recv_sig,SNR(n),'measured');
            ifo_temp =  Shao_GetIFO(recv_sig_add_noise,train_symbol,N,Ng);
            if ifo_df == ifo_temp
                cnt = cnt + 1;
            end
        end
        success_rate(n) = cnt/simu_times;
        failure_rate(n) = 1 - success_rate(n);
    end
%     figure
%     d = 1:length(SNR);
%     plot(SNR,success_rate(d));
end
function[recv_sig] = add_CFO(trans_sig,df,N)
%===========添加频偏===========
    recv_sig = zeros(1,length(trans_sig));
    for k=1:length(trans_sig)
        %添加频偏
        recv_sig(k) = trans_sig(k)*exp(sqrt(-1)*2*pi*df*k/N);
    end  
end

function [ifo_df] = Shao_GetIFO(recv_sig,train_symbol,N,Ng)
Ns = N+Ng;
R1 = zeros(1,2*Ns);
R2 = zeros(1,2*Ns);
R3 = zeros(1,2*Ns);
R4 = zeros(1,2*Ns);
A = train_symbol(1,1:N/4);
B = train_symbol(1,N/4+1:N/2);
for d = Ns+1:1:2*Ns
    for k = 0:N/4-1
        R1(d-Ns/2) = R1(d-Ns/2) + recv_sig(d+k)*conj(A(k+1));
    end
end
for d = Ns+1:1:2*Ns
    for k = 0:N/4-1
        R2(d-Ns/2) = R2(d-Ns/2) + recv_sig(d+k)*A(k+1);
    end
end
for d = Ns+1:1:2*Ns
    for k = 0:N/4-1
        R3(d-Ns/2) = R3(d-Ns/2) + recv_sig(d+k)*conj(B(k+1));
    end
end
for d = Ns+1:1:2*Ns
    for k = 0:N/4-1
        R4(d-Ns/2) = R4(d-Ns/2) + recv_sig(d+k)*B(k+1);
    end
end
[~,d1] = max(abs(R1));
[~,d2] = max(abs(R2));
[~,d3] = max(abs(R3));
[~,d4] = max(abs(R4));
% figure
% d = Ns+1:1:2*Ns;
% subplot(4,1,1)
% plot(d,abs(R1(d-Ns/2)));
% subplot(4,1,2)
% plot(d,abs(R2(d-Ns/2)));
% subplot(4,1,3)
% plot(d,abs(R3(d-Ns/2)));
% subplot(4,1,4)
% plot(d,abs(R4(d-Ns/2)));
%ifo_df = (N/2+d1-d2+d3-d4)/4
ifo_df = ((d2-d1)+(d3-d4))/4;
% d1 = d1 + Ns/2
% Ns+Ng
% d2 = d2 + Ns/2
% Ns+Ng+N/2
% d3 = d3 + Ns/2
% Ns+Ng+N/4
% d4 = d4 + Ns/2
% Ns+Ng+3*N/4
end
