function SelfTest(obj,config)
    N = config.N;
    Ng = config.Ng;
    Ns = N+Ng;
    SNR = 25;%config.SNR;
    simu_times = 1;%config.simu_times;
    [Data,TS] = obj.generate_data(N,Ng);
    SR = zeros(1,length(SNR));
    for n = 1:length(SNR)
        for s=1:simu_times
            ffo_df = 0;%-0.05 + 0.1*rand();
            ifo_df = 4;
            df = ffo_df + ifo_df;
            sData = add_CFO(Data,df,N);
            sData = awgn(sData,SNR(n),'measured');
            IFO = obj.IFOEstimator(sData(1,Ng+1:Ns),N,TS);
            if IFO == ifo_df
                SR(n) = SR(n) + 1;
            end
        end
        SR(n) = SR(n)/simu_times;
    end
    figure 
    d = 1:length(SNR);
    semilogy(SNR,1-SR(d),'-oc','LineWidth',1); 
    xlabel('SNR(dB)'); 
    ylabel('Probabilty of Failure (POF)'); 
    %legend('ren??');
end

function[recv_sig] = add_CFO(trans_sig,df,N)
    %===========Ìí¼ÓÆµÆ«===========
    recv_sig = zeros(1,length(trans_sig));
    parfor k=1:length(trans_sig)
        %Ìí¼ÓÆµÆ«
        recv_sig(k) = trans_sig(k)*exp(sqrt(-1)*2*pi*df*k/N);
    end  
end