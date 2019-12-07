classdef ProposedAlgo
    properties
        
    end
    methods
        function df = FFOEstimator(obj,N,recv)
            df = FFOAlgoProposed(N,recv);
        end
        function df = IFOEstimator(obj,recv_sig,N,TS)
            Q=0;
            recv = recv_sig;
            %==============´Ö¹À¼Æ=============
            for k=1:N-1
               Q1 = recv(k)*conj(TS(k));
               Q2 = recv(k+1)*conj(TS(k+1));
               Q = Q + conj(Q1)*Q2;
            end
            idf_coarse = N*angle(Q)/(2*pi);
            if idf_coarse>=0
                idf = floor(idf_coarse);
            else
                idf = ceil(idf_coarse);
            end
            %==============Ï¸¹À¼Æ==============
            W= N/4;
            P = zeros(1,2*W+1);
            E = 0;
            for g = idf-W:idf+W
                for n=1:N
                    P(g-idf+W+1) = P(g-idf+W+1) + recv(n)*conj(TS(n))*exp(-sqrt(-1)*2*pi*g*n/N);
                end
            end

            for n=1:N
                E = E+abs(recv(n))^2;
            end

            H = abs(P).^2/E^2;
            [~,df]= max(H);
            df = df +idf-W-1;
        end

        function [data,TS] = generate_data(obj,N,Ng)
            a = zeros(1,N/4);
            mu = N/4-1;
            for n=0:N/4-1
               a(n+1) =exp(4*sqrt(-1)*mu*pi*n*n/N);
            end
            A = ifft(a);
            B = conj(A(1,N/4:-1:1));
            C= zeros(1,N/4);
            for n=1:1:N/4 
                if mod(n,2)
                   C(n) = (-1)*B(n);
                else
                   C(n) = B(n);
                end
            end
            D=C;
            TS = [A B C D];
            data =[TS(1,N-Ng+1:N) TS];%Ê±ÓòÉÏµÄÑµÁ··ûºÅ
        end
    end
end