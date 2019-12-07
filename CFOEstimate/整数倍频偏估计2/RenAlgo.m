classdef RenAlgo
    methods
        function [data,TS] = generate_data(obj,N,Ng)
            a = zeros(1,N/2);
            for n=1:N/2
               % a(n) = 7*sqrt(2)*exp(2*sqrt(-1)*pi*n*n/N);%original
               a(n) =  exp(2*sqrt(-1)*pi*n*n/N);
            end
            A = ifft(a);
            %A=a;
            PNFactor= 2*(rand(1,N)>0.5)-1;
            signal = PNFactor.*[A A];
            TS = signal;
            data = [signal(1,N-Ng+1:N) signal];
        end
%         function [df] = FFOEstimator(obj,recv_sig,PNFactor,N)
%             P = 0;
%             pn = PNFactor;
%             r = recv_sig;
%             for n=1:N/2
%                 P = P + pn(n)*pn(n+N/2)*conj(r(n))*r(n+N/2);
%             end
%             df = angle(P)/pi;
%         end
        function [df] = IFOEstimator(obj,recv_sig,N,TS)
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
    end
end

function[recv_sig] = add_CFO(trans_sig,df,N)
    %===========Ìí¼ÓÆµÆ«===========
    recv_sig = zeros(1,length(trans_sig));
    for k=1:length(trans_sig)
        %Ìí¼ÓÆµÆ«
        recv_sig(k) = trans_sig(k)*exp(sqrt(-1)*2*pi*df*k/N);
    end  
end