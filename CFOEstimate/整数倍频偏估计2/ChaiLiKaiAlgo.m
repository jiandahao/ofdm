classdef ChaiLiKaiAlgo
    methods
        function [df] = get_IFO(obj,recv_sig,N,TS)
            P = zeros(1,N+1);
            recv = recv_sig;
            for g =-N/2:N/2
                for n=1:N
                    P(g+N/2+1) = P(g+N/2+1) + recv(n)*conj(TS(n))*exp(-sqrt(-1)*2*pi*g*n/N);
                end
            end

            H = abs(P).^2;
            [~,df]= max(H);
            df = df -(N/2+1);
        end

        function [signal,TS] = generate_data(obj,N,Ng)
            a = zeros(1,N/4);
            mu = 1;
            for n=0:N/4-1
               a(n+1) =exp(4*sqrt(-1)*mu*pi*n*n/N);
            end
            A = a;
            B = A.*A;
            C= -A;
            D= -B;
            signal = [A B C D];
            TS =signal;%Ê±ÓòÉÏµÄÑµÁ··ûºÅ
            signal = [signal(1,N-Ng+1:N) signal];
        end
    end
end