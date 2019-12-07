function [ M ] = cazacABCDForC( data_N_length)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
  P1=complex(0);
  R1=0; 
  P2=complex(0);
  R2=0; 
  R=0; 
  P=0; 
    recv = data_N_length;
    for m=1:64
    
    P1 = P1 + recv(m)*recv(128-m-1);
    R1 = R1 + real(recv(m))*real(recv(m))+imag(recv(m))*imag(recv(m));%power(abs(recv(d+m)),2);%
    P2 = P2 + recv(m+128)*recv(192+m);
    R2 = R2 + real(recv(m+128))*real(recv(m+128))+imag(recv(m+128))*imag(recv(m+128));%power(abs(recv(d+m+N/2)),2);


    end
    P1_realpart2 = real(P1)*real(P1);
    P1_imagpart2 = imag(P1)*imag(P1);
    absOf_P1 = sqrt(P1_realpart2+P1_imagpart2);
    P2_realpart2 =real(P2)*real(P2);
    P2_imagpart2 =imag(P2)*imag(P2);
    absOf_P2 = sqrt(P2_realpart2 + P2_imagpart2);
    P = absOf_P1 * absOf_P2;
    R = 0.5*(R1+R2)*R2;

    M=(P*P)/(R*R);

end