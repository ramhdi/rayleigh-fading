clear all; close all; clc; pkg load communications;
% Task 2 Jaringan Akses Nirkabel
%% Rama Rahardi / 18117026

%% initialization
tic
sr=128000.0; % Symbol rate
nd = 1e5;   % Number of symbols that simulates in each loop
fd = 30; %Hertz %doppler frequency
EbN0_fading = [0:5:25];%[0:5:25];

%% data generation
data=rand(1,nd)>0.5;

%% BPSK symbol mapping
x1 = data*2-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fading Channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EbN0 = EbN0_fading;  # Eb/N0 in dB
error = zeros(1,length(EbN0));
for i = 1:length(EbN0)
  disp(EbN0(i)); disp(i);
  %fad = cxn(nd, 1);
  %fad = fading(nd, fd, 1/sr)';
  fad1 = fading2(nd, fd, 1/sr);
  fad2 = fading2(nd, fd, 1/sr);
  %fad1 = ones(1,nd);
  %fad2 = ones(1,nd);
  n1 = 1/sqrt(2)*[randn(1,nd) + j*randn(1,nd)]; % white gaussian noise, 0dB variance
  n2 = 1/sqrt(2)*[randn(1,nd) + j*randn(1,nd)]; % white gaussian noise, 0dB variance
  y11 = x1.*fad1 + 10^(-EbN0(i)/20)*n1;
  y12 = x1.*fad2 + 10^(-EbN0(i)/20)*n2;
  
  y2 = zeros(1,nd); y3 = zeros(1,nd);
  %% selection diversity combining
  for n = 1:nd
    if (abs(y11(n)) >= abs(y12(n)))
      y2(n) = y11(n);
      y3(n) = y11(n)/fad1(n);
    else
      y2(n) = y12(n);
      y3(n) = y12(n)/fad2(n);
    end
  end
  
  y4 = real(y3)>0;
  
  error(i) = size(find([data - y4]),2);
end

%% plot daya
figure; hold on;
t = [1:5e1:nd/2]./sr;
lw = 'linewidth';
plot(t, 20*log10(abs(y11(1:5e1:nd/2))), '--', lw, 1.5);
plot(t, 20*log10(abs(y12(1:5e1:nd/2))), '--', lw, 1.5);
plot(t, 20*log10(abs(y2(1:5e1:nd/2))), 'k', lw, 1.75); hold off;
h = legend('Antena 1', 'Antena 2', 'SDC', 'location', 'SouthEast');
legend boxoff;
set(gca, 'fontsize', 14);
set(h, 'fontsize', 14);
xlabel('Time (s)');
ylabel('Received Signal Level (dB)');
grid on;

%% BER Calculation AWGN
ber_fading_sim = error/nd;
ber_fading_theory = (1/2).*(1-sqrt(10.^(EbN0/10)./(10.^(EbN0/10)+1)));

%% plot konstelasi
figure;
scatter(real(y3(1:1e2:nd)), imag(y3(1:1e2:nd)));
grid on; axis([-2 2 -2 2]);
set(gca, 'fontsize', 14);
xlabel('Real');
ylabel('Imaginary');

%% Plot BER vs EbN0
lw = 'linewidth';
figure; 
hold on;
semilogy(EbN0_fading, ber_fading_theory, '--b',"marker", "x", lw, 1.75);
semilogy(EbN0_fading, ber_fading_sim, '-b', "marker", "*", lw, 1.75);
hold off;
ylim([-1e-5 1e0]);
h = legend('Theoretical without diversity','Simulated with SDC', "location", 'NorthEast');
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (P_e)');
legend boxoff;
set(gca, 'fontsize', 14);
set(h, 'fontsize', 14);
grid on;

toc