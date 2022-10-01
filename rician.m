close all
clear all;
clc;

%BPSK Generation

N=10^6;                      
data=rand(1,N)>0.5;         
sym=2*data-1;  
colorstflex='kbgrm';
snr_dB=-5:2:20; %for multiple snr
K=[1 2 5 10 20 ]; 
tot_Power=1; %LOS and scatter path power

sim_ricean=zeros(1,length(snr_dB));
plotStyle={'-o','-*','-r','-s','-^'};
for index=1:length(K)
    k=K(index);
    s=sqrt(k/(k+1)*tot_Power);      % non central parameter
    sigma=tot_Power/sqrt(2*(k+1));
    
    %sim_BER_r=zeros(1,length(snr_dB));

    for i=1:length(snr_dB)
        n=(10^(-snr_dB(i)/20));  
        noise=1/sqrt(2)*(randn(1,N)+1i*randn(1,N));           % AwGN noise with mean =0 and var =1
        h=((sigma*randn(1,N)+s) + 1i*(randn(1,N)*sigma +0)); %Rician fading
     
        y_ricean=h.*sym + noise*n;  % received signal through rician channel
     
     
     
     
       %coherent reciver 
       y_rcv=y_ricean./h; %equalisation of received signalby channel information at receiver.
       yHat=real(y_rcv)>0;
     
       %receiver for AWGN channel
       sim_ricean(i)=sum(xor(data,yHat));
   
    end
    sim_ricean=sim_ricean/N;
    snr=10.^(snr_dB/10); % for linear scale
    semilogy(snr_dB,sim_ricean,plotStyle{index},'color',colorstflex(index));
    hold on;
    legendInfo{index} =['K= ',num2str(K(index))];
   
    
     
end


grid on;
axis([-5 20 10^-5 1]);
legend(legendInfo);
title('SNR Vs BER over Rician Channels with AWGN noise');
xlabel('SNR(dB)');
ylabel('Bit Error Rate');

     
     
    
    
    