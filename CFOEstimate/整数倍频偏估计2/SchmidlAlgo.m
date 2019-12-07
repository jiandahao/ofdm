classdef SchmidlAlgo
    properties 
    end    
    
    methods
        function[signal,v] = generate_data(obj,N,Ng)
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
            %cp_train1 = sch;
            %=============生成训练序列2======================
%            buf=QAMTable(randi([0,3],N/2,1)+1); %1x(N/2)的矩阵
            %train2=zeros(1,N);
            train2=QAMTable(randi([0,3],N,1)+1);
%             index = 1; 
%             for n=1:2:N 
%                  train2(n)=buf(index); 
%                  index=index+1; 
%             end; 
            sch = ifft(train2,N);   %[A A]的形式 
            cp_train2 = [sch(1,N-Ng+1:N) sch];
            %cp_train2 = sch;
            %==============计算差分序列====================
            v =sqrt(2)*train2./train1;
            transmit_data = [cp_train1 cp_train2];
            signal = transmit_data;
        end
        function [IFO_df] =  IFOEstimator(obj,r1,r2,N,v)
%             recv = recv_sig;
%             Ns = N + Ng;
%             r1 = recv(1,Ng+1:Ns);%获取第一个OFDM符号
%             r2 = recv(1,N+Ng+1:2*Ns);%获取第二个OFDM符号
%             R1 = fft(r1);
%             R2 = fft(r2);
%             
%             for g = -N/4:N/4
%                  B1 = zeros(1,N);
%                  B2 = zeros(1,N);
%                 for k = 1:2:N
%                     if k+2*g <=0
%                          B1(k) = conj(R1(k+2*g+N))*conj(v(k))*R2(k+2*g+N);
%                     end
%                     if (k+2*g<=N) &&  (k+2*g>0)
%                         B1(k+N/2+1) = conj(R1(k+2*g))*conj(v(k))*R2(k+2*g);
%                     end
%                     if k+2*g>N
%                         B1(k+N/2+1) =conj(R1(k+2*g-N))*conj(v(k))*R2(k+2*g-N);
%                     end
%                     B2(k) = abs(R2(k))^2;
%                 end
%                 B(g+N/4+1) = (abs(sum(B1)))^2/(2*sum(B2)^2);
%             end
%             [~,IFO_df] = max(B);
%             IFO_df = 2*(IFO_df -(N/4+1));
%             end
            %r1 = recv(1,Ng+1:Ns);%获取第一个OFDM符号
            %r2 = recv(1,N+2*Ng+1:2*Ns);%获取第二个OFDM符号
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
    end
end