EbN0 = -1:8; % power efficiency range in dB

% ======================================================================= %
% Task 1
% ======================================================================= %
figure(1);

% Load necessary plot data from .mat files.
data = load('./Task 1/uncoded_QPSK.mat');
BER_uncoded = data.BER;

data = load('./Task 1/eps2_QPSK_hard_rx.mat');
eps2_hard = data.BER;

data = load('./Task 1/eps2_QPSK_soft_rx.mat');
eps2_soft = data.BER;

% Calculate the upper bound on the soft decoder.
data = load('./Task 1/trellis.mat');
N = 20;
trellis = data.trellis;
code_rate = 1/2;
spect = distspec(trellis, N); 
dfree = spect.dfree; 
weights = spect.weight; 
Pb_soft_bound = zeros(size(EbN0));
EbN0_linear = 10.^(EbN0 ./ 10);
for i = 1:length(EbN0)
    d_values = (dfree:dfree+N-1);
    Q_values = qfunc(sqrt(2 .* d_values .* code_rate .* EbN0_linear(i)));
    Pb_soft_bound(i) = sum(weights .* Q_values);
end

% Generate the plot
semilogy(EbN0, BER_uncoded, 'o', 'DisplayName','Uncoded (simulation)');
hold on;
semilogy(EbN0, eps2_soft, '-o', 'DisplayName', 'Coded, soft Rx');
semilogy(EbN0, eps2_hard, '-x', 'DisplayName', 'Coded, hard Rx');
% Soft decoding BER upper bound
semilogy(EbN0, Pb_soft_bound, '-*', 'DisplayName', 'Soft Decoder Upper Bound');

% Plot formatting
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
legend show;
ylim([1e-4 1]);
title('BER vs E_b/N_0 using QPSK');
hold off;

% ======================================================================= %
% Task 2
% ======================================================================= %
figure(2);

% Load necessary plot data from .mat files.
data = load('./Task 2/BER_uncoded.mat');
BER_uncoded = data.BER;

data = load('./Task 2/encoder1simu.mat');
eps1 = data.encoder1simu;
data = load('./Task 2/eps1ub.mat');
eps1_UB = data.Pb_soft_bound;

data = load('./Task 2/encoder2simu.mat');
eps2 = data.encoder2simu;
data = load('./Task 2/eps2ub.mat');
eps2_UB = data.Pb_soft_bound;

data = load('./Task 2/encoder3simu.mat');
eps3 = data.encoder3simu;
data = load('./Task 2/eps3ub.mat');
eps3_UB = data.Pb_soft_bound;

% Generate the plot
semilogy(EbN0, eps1, 'b-', 'DisplayName', '$\varepsilon_1$ (simulation)');
hold on;
semilogy(EbN0, eps2, 'r-', 'DisplayName', '$\varepsilon_2$ (simulation)');
semilogy(EbN0, eps3, 'k-', 'DisplayName', '$\varepsilon_3$ (simulation)');
semilogy(EbN0, eps1_UB, 'b--','DisplayName', '$\varepsilon_1$ (upper bound)');
semilogy(EbN0, eps2_UB, 'r--','DisplayName', '$\varepsilon_2$ (upper bound)');
semilogy(EbN0, eps3_UB, 'k--','DisplayName', '$\varepsilon_3$ (upper bound)');
% Uncoded case
semilogy(EbN0, BER_uncoded, 'o', 'DisplayName','Uncoded (simulation)','Color', [0, 0.4470, 0.7410]);

% Plot formatting
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
set(legend, 'interpreter', 'latex', fontsize = 14);
legend show;
ylim([1e-4 1]);
title('BER vs E_b/N_0 using QPSK');
hold off;

% ======================================================================= %
% Task 3
% ======================================================================= %
figure(3);

% Load necessary plot data from .mat files.
data = load('./Task 3/BPSK_uncoded.mat');
BPSK_uncoded = data.BER;
data = load('./Task 3/QPSK_uncoded.mat');
QPSK_uncoded = data.BER;
data = load('./Task 3/AMPM_uncoded.mat');
AMPM_uncoded = data.BER;

data = load('./Task 3/eps3_BPSK.mat');
eps3_BPSK = data.BER;
data = load('./Task 3/eps3_QPSK.mat');
eps3_QPSK = data.BER;
data = load('./Task 3/eps4_AMPM.mat');
eps4_AMPM = data.BER;

semilogy(EbN0, eps3_BPSK, 'b-', 'DisplayName', '$\varepsilon_3$ with BPSK');
hold on;
semilogy(EbN0, eps3_QPSK, 'r-', 'DisplayName', '$\varepsilon_3$ with QPSK');
semilogy(EbN0, eps4_AMPM, 'k-', 'DisplayName', '$\varepsilon_4$ with AMPM');
semilogy(EbN0, BPSK_uncoded, 'b--','DisplayName', 'Uncoded BPSK');
semilogy(EbN0, QPSK_uncoded, 'r--','DisplayName', 'Uncoded QPSK');
semilogy(EbN0, AMPM_uncoded, 'k--','DisplayName', 'Uncoded AMPM');

% Plot formatting
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
set(legend, 'interpreter', 'latex', fontsize = 14);
legend show;
ylim([1e-4 1]);
title('BER vs E_b/N_0 using various constellations');
hold off;
