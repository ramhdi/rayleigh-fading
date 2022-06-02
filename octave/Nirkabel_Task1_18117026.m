clear all; close all; clc; pkg load communications;
% Task 1 Jaringan Akses Nirkabel
%% Rama Rahardi / 18117026

tic
%% initialization
sr=128000.0; % Symbol rate
ml=1;        % ml:Number of modulation levels (BPSK:ml=1, QPSK:ml=2, 16QAM:ml=4)
br=sr .* ml; % Bit rate
nd = 10^5;   % Number of symbols that simulates in each loop
fd = 30; %Hertz %doppler frequency
EbN0_awgn = [0:2:8];
EbN0_fading = [0:5:25];

%% data generation
data=rand(1,nd*ml)>0.5;

%% BPSK modulation
x1 = data*2-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AWGN Channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EbN0 = EbN0_awgn;  # Eb/N0 in dB
for i = 1:length(EbN0)
  n = 1/sqrt(2)*[randn(1,nd) + j*randn(1,nd)]; % white gaussian noise, 0dB variance
  y1 = x1 + 10^(-EbN0(i)/20)*n;
  
  %% receiver decision (demodulation)
  y2 = real(y1)>0;
  
  %% error calculation
  error(i) = size(find([data- y2]),2);
end
%% BER Calculation AWGN
ber_awgn_sim = error/nd;
ber_awgn_theory = qfunc(sqrt(2*10.^(EbN0/10)));

%% plot konstelasi
figure;
scatter(real(y1(1:1e2:nd)), imag(y1(1:1e2:nd)));
grid on; axis([-2 2 -2 2]);
set(gca, 'fontsize', 14);
xlabel('Real');
ylabel('Imaginary');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fading Channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EbN0 = EbN0_fading;  # Eb/N0 in dB
for i = 1:length(EbN0)
  %fad = cxn(nd, 1);
  %fad = fading(nd, fd, 1/sr)';
  fad = fading2(nd, fd, 1/sr);
  n = 1/sqrt(2)*[randn(1,nd) + j*randn(1,nd)]; % white gaussian noise, 0dB variance
  y1 = x1.*fad + 10^(-EbN0(i)/20)*n;
  
  %% receiver decision (demodulation)
  %%% channel equalization
  y2 = y1./fad;
  y3 = real(y2)>0;
  
  %% error calculation
  error(i) = size(find([data- y3]),2);
end

%% BER Calculation AWGN
ber_fading_sim = error/nd;
ber_fading_theory = (1/2).*(1-sqrt(10.^(EbN0/10)./(10.^(EbN0/10)+1)));

%% plot konstelasi
figure;
scatter(real(y2(1:1e2:nd)), imag(y2(1:1e2:nd)));
grid on; axis([-2 2 -2 2]);
set(gca, 'fontsize', 14);
xlabel('Real');
ylabel('Imaginary');

lw = 'linewidth';
figure; 
hold on;
semilogy(EbN0_awgn, ber_awgn_theory, '--r',"marker", "x", lw, 1.75);
semilogy(EbN0_awgn, ber_awgn_sim, '-r', "marker", "*", lw, 1.75);
semilogy(EbN0_fading, ber_fading_theory, '--b',"marker", "x", lw, 1.75);
semilogy(EbN0_fading, ber_fading_sim, '-b', "marker", "*", lw, 1.75);
hold off;
ylim([-10^-4 10^0]);
h = legend('AWGN Theory','AWGN Simulated', 'AWGN Fading Theory','AWGN Fading Simulated', "location", 'NorthEast');
xlabel('Eb/N_0 (dB)');
ylabel('Bit Error Rate (P_e)');
legend boxoff;
set(gca, 'fontsize', 14);
set(h, 'fontsize', 14);
grid on;

toc