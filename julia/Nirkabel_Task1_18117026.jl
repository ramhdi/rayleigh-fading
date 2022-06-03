#clear all; close all; clc; pkg load communications;
# Task 1 Jaringan Akses Nirkabel
# Rama Rahardi / 18117026

#tic
using Plots
include("qfunc.jl")
include("fading2.jl")

# initialization
sr = 128000.0; # Symbol rate
ml=1;        # ml:Number of modulation levels (BPSK:ml=1, QPSK:ml=2, 16QAM:ml=4)
br=sr .* ml; # Bit rate
nd = 10^5;   # Number of symbols that simulates in each loop
fd = 30; # Hertz %doppler frequency
EbN0_awgn = collect(0:1:10)';
EbN0_fading = collect(0:1:25)';
error_awgn = zeros(1,length(EbN0_awgn));
error_fading = zeros(1,length(EbN0_fading));
ber_awgn_sim = zeros(1,length(EbN0_awgn));
ber_fading_sim = zeros(1,length(EbN0_fading));
N_tries = 5; # number of simulation trials

# data generation
data = rand(1,nd*ml) .> 0.5;

# BPSK modulation
x1 = (data*2).-1;

# AWGN Channel
EbN0 = EbN0_awgn;  # Eb/N0 in dB

for tr = 1:N_tries
  for i = 1:length(EbN0)
    n = 1/sqrt(2) * (randn(1,nd) + im*randn(1,nd)); # white gaussian noise, 0dB variance
    local y1 = x1 + 10^(-EbN0[i]/20) * n;
    
    # receiver decision (demodulation)
    local y2 = real(y1) .> 0;
    
    # error calculation
    #error(i) = size(find([data - y2]),2);
    error_awgn[i] = sum(abs.(data.-y2))
  end
  global ber_awgn_sim += error_awgn/nd;
end

# BER Calculation AWGN
ber_awgn_sim = ber_awgn_sim/N_tries;
ber_awgn_theory = qfunc.(sqrt.(2 * 10 .^ (EbN0/10)));

# plot BER
plot(
  EbN0_awgn',
  [ber_awgn_theory' ber_awgn_sim'],
  xlims=(0,8),
  ylims=(1e-4, 1e-1),
  yaxis=:log,
  xlabel="Eb/N0 (dB)",
  ylabel="Bit Error Rate",
  label=["Theoretical" "Simulated"],
  title="BPSK AWGN"
)

# Fading Channel
EbN0 = EbN0_fading;  # Eb/N0 in dB
for tr = 1:N_tries
  for i = 1:length(EbN0)
    #fad = cxn(nd, 1);
    #fad = fading(nd, fd, 1/sr)';
    fad = fading2(nd, fd, 1/sr);
    n = 1/sqrt(2)*(randn(1,nd) + im*randn(1,nd)); # white gaussian noise, 0dB variance
    local y1 = x1.*fad + 10^(-EbN0[i]/20)*n;
    
    # receiver decision (demodulation)
    # channel equalization
    local y2 = y1./fad;
    local y3 = real(y2) .> 0;
    
    # error calculation
    error_fading[i] = sum(abs.(data.-y3));
  end
  global ber_fading_sim += error_fading/nd;
end

# BER Calculation AWGN
ber_fading_sim = ber_fading_sim/N_tries;
ber_fading_theory = (1/2) .* (1 .- sqrt.(10 .^ (EbN0/10) ./ (10 .^ (EbN0/10) .+ 1)));

# plot konstelasi
plot(
  EbN0_fading',
  [ber_fading_theory' ber_fading_sim'],
  xlims=(0,25),
  ylims=(1e-4, 1e0),
  yaxis=:log,
  xlabel="Eb/N0 (dB)",
  ylabel="Bit Error Rate",
  label=["Theoretical" "Simulated"],
  title="BPSK Fading"
)

#=
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
=#
#toc