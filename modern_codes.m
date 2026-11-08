% ======================================================================= %
% SSY125 Project Bonus task
% ======================================================================= %
clc
clear

% ======================================================================= %
% Simulation Options
% ======================================================================= %
N = 32400;  % simulate N bits each transmission (one block) adjusted to fit LDPC codewords.
maxNumErrs = 100; % get at least 100 bit errors (more is better)
maxNum = 1e6; % OR stop if maxNum bits have been simulated
EbN0 = -1:8; % power efficiency range

% ======================================================================= %
% Other Options
% ======================================================================= %
% n and k are fetched from the DVB-S2 standard.
n = 64800; 
k = 32400;  


code_rate = k / n; % Code rate (R_c = 1/2)
modOrder = 4; % QPSK (log2(M) = 2)
bpsymb = log2(modOrder);

% Generate LDPC matrix using DVB-S2 standard
H = dvbs2ldpc(code_rate);
encodercfg = ldpcEncoderConfig(H);
decodercfg = ldpcDecoderConfig(encodercfg);

% ======================================================================= %
% Simulation Chain
% ======================================================================= %
BER = zeros(1, length(EbN0)); % pre-allocate a vector for BER results

for i = 1:length(EbN0)
    totErr = 0;  % Number of errors observed
    num = 0; % Number of bits processed

    while((totErr < maxNumErrs) && (num < maxNum))
        % ===================================================================== %
        % Begin processing one block of information
        % ===================================================================== %
        % [SRC] generate N information bits
        numCodewords = N / k; % Calculate number of codewords
        info_bits = randi([0 1], k, numCodewords); % Generate information bits in the shape k x N/k matrix

        % Encode using the parity check matrix H from earlier
        encodedData = ldpcEncode(info_bits,encodercfg);

        % Flatten for ”transmission"
        tx_bits = encodedData(:);

        % [MOD] symbol mapper
        const = [(1 + 1i), (1 - 1i), (-1 + 1i), (-1 - 1i)]*(1/sqrt(2));

        % Reshape and "send" the bits in packages of bpsymb
        m = reshape(tx_bits, bpsymb, [])';
        m_idx = bi2de(m, 'left-msb')'+1;
        x = const(m_idx);

        % [CHA] add Gaussian noise
        E_s = mean(abs(const).^2);
        EbN0_linear = 10^(EbN0(i) / 10);
        noise_power = E_s / (code_rate * bpsymb * EbN0_linear);
        sigma_sq = noise_power/2;

        N_noise = length(x);
        noise = sqrt(noise_power / 2) * (randn(1,N_noise) + 1i*randn(1,N_noise));

        y = x + noise;

        % [SR] Soft Receiver

        % Maximum number of decoder iterations:
        maxIterations = 50;

        % Compute the log-likelihood ratios
        llr_real = 2 * real(y) / sigma_sq;
        llr_imag = 2 * imag(y) / sigma_sq;

        % Combine LLRs into a single vector
        llr = reshape([llr_real; llr_imag], [], 1);

        % Reshape LLR vector to fit ldpcDecode() standard
        llr = reshape(llr, k/code_rate, []);

        % [decoded_bits, ~, ~] = ldpcDecode(llr, decodercfg, maxIterations);
        [decoded_bits, ~, ~] = customLDPCDecode(llr, H, maxIterations);



        % ===================================================================== %
        % End processing one block of information
        % ===================================================================== %
        BitErrs = sum(xor(decoded_bits(1:N), info_bits(:))); % count the bit errors and evaluate the bit error rate
        totErr = totErr + BitErrs;
        num = num + N;

        disp(['+++ ' num2str(totErr) '/' num2str(maxNumErrs) ' errors. '...
            num2str(num) '/' num2str(maxNum) ' bits. Projected error rate = '...
            num2str(totErr/num, '%10.1e') '. +++']);
    end
    BER(i) = totErr/num;
end

data = load("./Task 3/eps3_QPSK.mat");
BER_system_2 = data.BER;

data = load("./Task 1/uncoded_QPSK.mat");
BER_uncoded = data.BER;

semilogy(EbN0, BER, '-o', DisplayName="LDPC");
hold on;
semilogy(EbN0, BER_system_2, 'o', DisplayName="System 2");
semilogy(EbN0, BER_uncoded, '--', DisplayName="Uncoded QPSK")

% Plot formatting
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
legend show;
ylim([1e-4 1]);
xlim([-1,8]);
title('BER vs E_b/N_0 using QPSK');
hold off;

% ======================================================================= %
% End
% ======================================================================= %