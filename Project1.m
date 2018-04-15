close all
clear
clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question A

% P = 100;
% N = 10^5;
% A = 1;
% T = 0.1;
% beta = 0.02;
% SNR_dB = 5;
% SNR = 10.^(SNR_dB/10);
% BER = 0;
% 
% for p=1:P
%     x_est = zeros(1,N);
%     x = 2*round(rand(1,N))-1;
%     
%     x_i = zeros(1,N/2);
%     x_q = zeros(1,N/2);
%     x_i(1) = -(-1)*1;        %x_i(0)=-x_q(n-1)*x(2n-1)=-x_q(-1)*x(-1)=-(-1)*1
%     x_q(1) = -x_i(1)*x(1);   %x_q(0)=-x_i(0)*x(0)=-1*x(0)
%     
%     for n=2:N/2
%         x_i(n) = -x_q(n-1)*x(2*(n-1));
%         x_q(n) = -x_i(n)*x(2*n-1);
%     end
%     
%     z = x_i + 1j*x_q;
%     
%     n = randn(1,N/2) + 1j*randn(1,N/2);
%     y = A*T*z + sqrt(T^2*A^2/SNR)*n;
%     
%     y_est_r = sign(real(y));
%     y_est_i = sign(imag(y));
%     
%     x_est(1) = -y_est_i(1)*y_est_r(1);
%     for k=2:N/2
%         x_est(2*(k-1)) = -y_est_i(k-1)*y_est_r(k);
%         x_est(2*k-1) = -y_est_i(k)*y_est_r(k);
%     end
%     BER = BER + sum(x_est ~= x);
%     
% end
% 
% BER = BER/(N*P);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question B

% P = 100;
% N = 10^5;
% A = 1;
% T = 0.1;
% SNR_dB = 5:12;
% SNR = 10.^(SNR_dB/10);
% BER = zeros(1,length(SNR));
% 
% for j=1:length(SNR)
%     for p=1:P
%         x_est = zeros(1,N);
%         x = 2*round(rand(1,N))-1;
%         
%         x_i = zeros(1,N/2);
%         x_q = zeros(1,N/2);
%         x_i(1) = -(-1)*1;        %x_i(0)=-x_q(n-1)*x(2n-1)=-x_q(-1)*x(-1)=-(-1)*1
%         x_q(1) = -x_i(1)*x(1);   %x_q(0)=-x_i(0)*x(0)=-1*x(0)
%         
%         for n=2:N/2
%             x_i(n) = -x_q(n-1)*x(2*(n-1));
%             x_q(n) = -x_i(n)*x(2*n-1);
%         end
%         
%         z = x_i + 1j*x_q;
%         
%         n = randn(1,N/2) + 1j*randn(1,N/2);
%         y = A*T*z + sqrt(T^2*A^2/SNR(j))*n;
%         
%         y_est_r = sign(real(y));
%         y_est_i = sign(imag(y));
%         
%         x_est(1) = -y_est_i(1)*y_est_r(1);
%         for k=2:N/2
%             x_est(2*(k-1)) = -y_est_i(k-1)*y_est_r(k);
%             x_est(2*k-1) = -y_est_i(k)*y_est_r(k);
%         end
%         BER(j) = BER(j) + sum(x_est ~= x);
%     end
%     BER(j) = BER(j)/(N*P);
% end
% 
% figure, semilogy(SNR_dB,BER), xlabel('SNR in dB'), ylabel('Bit Error Rate'), hold on, semilogy(SNR_dB,qfunc(sqrt(SNR)),'r')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question D

% rng('default')
% SNR_dB = 12;
% SNR = 10.^(SNR_dB/10);
% N = 10^5;
% P = 1000;
% A = 1;
% T = 0.1;
% beta = 0.02;
% s1 = [A*sqrt(T); 0];
% s2 = [-2*A*sqrt(T)*1j/pi; A*sqrt(T)*sqrt(pi^2 -4)/pi];
% x = 2*round(rand(1,N))-1;
% n1 = sqrt(2*A^2*T/SNR)*(randn(1,N) + 1j*randn(1,N));
% n2 = sqrt(2*A^2*T/SNR)*(randn(1,N) + 1j*randn(1,N));
% phi = zeros(1,N);
% r = zeros(2,N);
% BER_VA = 0;
%
% for j=1:P
%     phi = zeros(1,N);
%     n1 = sqrt(2*A^2*T/SNR)*(randn(1,N) + 1j*randn(1,N));
%     n2 = sqrt(2*A^2*T/SNR)*(randn(1,N) + 1j*randn(1,N));
%
%     for n=1:N
%         phi(n+1) = phi(n) + x(n)*pi/2;
%
%         if(x(n)==1)
%             r(:,n) =  s1.*exp(1j*phi(n)) + [n1(n); n2(n)];
%         else
%             r(:,n) =  s2.*exp(1j*phi(n)) + [n1(n); n2(n)];
%         end
%     end
%
%     x_est = ViterbiAlgorithm(N,s1,s2,r);
%     BER_VA = BER_VA + sum(x~=x_est);
% end
%
% BER_VA = BER_VA/(N*P);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question E

% N = 10^5;
% P = 5;
% A = 1;
% T = 0.1;
% SNR_dB = 5:12;
% SNR = 10.^(SNR_dB/10);
% x = 2*round(rand(1,N))-1;
% BER_VA = zeros(1,length(SNR));
% r = zeros(2,N);
%
% s1 = [A*sqrt(T); 0];
% s2 = [-2*A*sqrt(T)*1j/pi; A*sqrt(T)*sqrt(pi^2 -4)/pi];
%
% for i=1:length(SNR)
%     for j=1:P
%         phi = zeros(1,N);
%         n1 = sqrt(2*A^2*T/SNR(i))*(randn(1,N) + 1j*randn(1,N));
%         n2 = sqrt(2*A^2*T/SNR(i))*(randn(1,N) + 1j*randn(1,N));
%
%         for n=1:N
%             phi(n+1) = phi(n) + x(n)*pi/2;
%
%             if(x(n)==1)
%                 r(:,n) =  s1.*exp(1j*phi(n)) + [n1(n); n2(n)];
%             else
%                 r(:,n) =  s2.*exp(1j*phi(n)) + [n1(n); n2(n)];
%             end
%         end
%
%         x_est = ViterbiAlgorithm(N,s1,s2,r);
%         BER_VA(i) = BER_VA(i) + sum(x~=x_est);
%     end
%     BER_VA(i) = BER_VA(i)/(N*P);
% end
%
% figure, semilogy(SNR_dB,BER), xlabel('SNR in dB'), ylabel('Bit Error Rate'), hold on, semilogy(SNR_dB,BER_VA), hold on, semilogy(SNR_dB,qfunc(sqrt(SNR)))
% leg = legend('ML Simulation','Viterbi Simulation','$Q(\sqrt{SNR})$');
% set(leg,'Interpreter','latex'), set(leg,'FontSize',12);

