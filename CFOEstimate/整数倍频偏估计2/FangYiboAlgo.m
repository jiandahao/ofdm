% @InProceedings{fang2012novel,
%   author       = {Fang, Yibo and Zhang, Zuotao and Liu, Guanghui},
%   title        = {A novel synchronization algorithm based on CAZAC sequence for OFDM systems},
%   booktitle    = {2012 8th International Conference on Wireless Communications, Networking and Mobile Computing (WiCOM)},
%   year         = {2012},
%   pages        = {1--4},
%   organization = {IEEE},
% }
classdef FangYiboAlgo
    methods
        function [signal,v,mu] = generate_data(obj,N,Ng)
            a = zeros(1,N/2);
            for n=1:N/2
                a(n) = exp(2*sqrt(-1)*pi*n*n/N);
            end
            A = a;
            m = 1.2*(rand(1,N/2))-0.2;
            v = exp((sqrt(-1)*pi).*m);
            signal = [v.*A A];
            signal = [signal(1,N-Ng+1:N) signal];
            mu = A;
        end
        function [df] = IFOEstimator(obj,recv_sig,N,v,mu)
            recv = recv_sig;
            B = zeros(1,N/2);
            F=  zeros(1,N/2+1);
            for g=-N/4:N/4
                for k = 0:N/2-1
                    B(k+1)  = (recv(k+1)*conj(v(k+1))+recv(1+k+N/2))*conj(mu(mod(k+g,N/2)+1));
                end

                F(g+N/4+1) = abs(sum(B));
            end
            [~,df] = max(F);
            df = 2*(df - N/4 -1);
        end
    end
end