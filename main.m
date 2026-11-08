% ======================================================================= %
% SSY125 Project
% ======================================================================= %
clc
clear

% ======================================================================= %
% Simulation Options
% ======================================================================= %
N = 3e5;  % simulate N bits each transmission (one block)
maxNumErrs = 100; % get at least 100 bit errors (more is better)
maxNum = 1e6; % OR stop if maxNum bits have been simulated
EbN0 = -1:8; % power efficiency range in dB

% ======================================================================= %
% Other Options
% ======================================================================= %

constellation_option = 'AMPM'; % Simply set to 'QPSK' or 'BPSK' or 'AMPM'
encoder_option = 4; % Set to i for encoder eps_i or 0 for uncoded case.
hard_decoding = false; % Set to true for hard decoding and false for soft decoding.

% ======================================================================= %
% Our own functions
% ======================================================================= %
function encodedBits = cnEncoder(inputBits, g1, g2)
    % inputBits: A vector of input bits (1xN) g1, g2: Generator polynomials
    % represented in binary. Leftmost bit is the D^0 term!
    
    % MATLAB automatically zero pads the information bits when using conv(). So
    % the shift register is automatically initilised to 00 state and will
    % automatically be terminated in the 00 state too.
    
    % Get the coded bits through the binary convolution
    c_1 = mod(conv(inputBits, g1),2);
    c_2 = mod(conv(inputBits, g2),2);
    
    % Reshape to interweave the encoded bitvectors
    encodedBits = reshape([c_1; c_2], 1, []);
    encodedBits = encodedBits';
end

function encoded_bits = eps4(input_bits)
    % Function to encode using Rate-2/3 convolutional encoder (eps_4)
    % input_bits: 1xN matrix
    % encoded_bits: 1xN interweaved bitstream
    
    % Reshape into a N/2 x 2 matrix
    % each row is [u^(1)_i, u^(2)_i]
    input_pairs = reshape(input_bits.',2,[]);
    input_pairs = input_pairs';

    % Number of input bit-pairs
    num_pairs = size(input_pairs,1);

    % We know that c^(2)_i and c^(3)_i are just copies of u^(1)_i and
    % u^(2)_i respectively, so we can assign those.
    code_bits = zeros(num_pairs, 3);
    code_bits(:,2) = input_pairs(:,1);
    code_bits(:,3) = input_pairs(:,2);

    % Initialize shift registers (3 memory elements)
    state = [0, 0, 0];

    % Encoding process
    for i = 1:num_pairs

        % Current input bits
        u1 = input_pairs(i, 1);
        u2 = input_pairs(i, 2);

        % Compute c^(1)_i based on the connection diagram
        c1 = state(3);

        % Store output
        code_bits(i, 1) = c1;

        % Update the states/shift register
        state = [state(3), mod(state(1) + u2, 2), mod(state(2) + u1, 2)];
    end

    % Transform the output into a single long bitstream:
    encoded_bits = reshape(code_bits.', 1, []);
end

function trellis = create_trellis(constraint_length, generators)
    % Calculate the degree of genrator polynomials.
    degree = max(constraint_length) - 1;
    nr_of_states = 2^degree;
    nr_of_inputs = 2; % As this is just for binary structures.

    % Initialise matrices to keep track of state transitions and outputs
    % Both matrices store the values in decimal base.
    nextStates = zeros(nr_of_states,nr_of_inputs);
    outputs = zeros(nr_of_states,nr_of_inputs);

    for state = 0:nr_of_states-1
        for input_bit = 0:nr_of_inputs-1
            % A lot is going on here on the next line.
            % bitshift(state, 1) can take state as a decimal value, for
            % example state = 2 and converts it implicitly to (1 0). Then
            % shifts it by one to the left, so (1 0 0). 
            % Then we add the output bit, for example 1. This gives us
            % (1 0 1). By using 2^degree - 1 which would be (1 1) in this
            % example it will output a new bit-sequence with ones only
            % where both (1 0 1) and (1 1) has ones. So (0 0 1). This
            % ensures our output is always within the bounds of our number
            % of states, we read it as (0 1). 
            next_state = bitand(bitshift(state, 1) + input_bit, 2^degree - 1);
            % Now we construct the matrix at the next state index using
            % this calculated state from above.
            nextStates(state+1, input_bit+1) = next_state;

            % ---------------
            % Compute outputs:
            % ---------------
            % Start initializing a vector of output bits. (Will later be
            % converted to decimal and stored in outputs matrix).
            % We have one output bit for every generator polynomial.
            output_bits = zeros(1,size(generators,1));
            % Loop through all generator polynomials
            for g_idx = 1:size(generators,1)
                g = fliplr(generators(g_idx,:)); % Our current generator polynomial
                % The generator polynomial should be loaded in binary
                % for example 111. This defines which of the memory bits we
                % want to use. 
                % We create the shift register. In state 10 for example we
                % will get 100 from bitshift. Then our input bit is also a
                % part of the register in some sense. So we add that and
                % get 101 (if we input 1). 
                reg = dec2bin(bitshift(state,1) + input_bit, constraint_length); % The shift register. Converted to binary strings

                % dec2bin outputs strings instead of arrays so this will
                % convert to binary arrays.
                reg = reg - '0';

                % The dot product will extract the shift-register bits that
                % are relevant using the generator polynomial.
                dot_result = mod(dot(g, reg), 2); % Perform mod-2 dot-product

                % The mod-2 dot product returns, well, the sum of the bits and thus
                % the output bit for this generator polynomial.
                output_bits(g_idx) = dot_result; 
            end
            % Now we have our output bits, the same number of bits as
            % number of generator polynomials. Now lets convert them into
            % decimal representation and store in outputs matrix.
            output_bits = num2str(output_bits);
            outputs(state + 1, input_bit + 1) = bin2dec(output_bits);
         end
     end
    % Construct the actual trellis structure for output.
    trellis.numInputSymbols = nr_of_inputs;
    trellis.numOutputSymbols = 2^size(generators,1);
    trellis.numStates = nr_of_states;
    trellis.nextStates = nextStates;
    trellis.outputs = outputs;


end

function decodedBits = viterbi(receivedBits, trellis)
    % recievedBits is a 1D matrix of bits. Could also change to take in the
    % symbols directly, for example pairs of bits for QPSK. 
    % trellis is the Trellis structure from the other function.
    
    % Fetch initial parameters
    bpsymb = log2(trellis.numOutputSymbols);
    numSteps = length(receivedBits)/bpsymb;
    numStates = trellis.numStates;
    numInputs = trellis.numInputSymbols;
    
    % Group the recieved bits into packets of size bpsymb and transpose
    % rS is now a matrix where every row is a "symbol"
    receivedSymbols = reshape(receivedBits, bpsymb, []).';
    
    % The trellis.outputs is a decimal matrix, so to compute Hamming distances
    % later on we need to convert to binary again. 
    % So, we extract the binary values:
    outputsBinary = de2bi(trellis.outputs, bpsymb, 'left-msb');
    % Now this is a (numStates*numInputs) x bpsymb matrix. Which means we are
    % mixing up states and input bits. To fix this, we reshape.
    % This new matrix will have dimensions:
    % numStates x numInputs x bpsymb
    possibleOutputsMatrix = reshape(outputsBinary, numStates, numInputs, bpsymb);
    
    
    % Expand our recieved "symbols" with some dummy-dimensions to facilitate
    % matrix operations and calculations. Will have dimensions:
    % nr of steps x 1 x 1 x bpsymb
    receivedExpanded = reshape(receivedSymbols, size(receivedSymbols,1), 1, 1, bpsymb);
    
    % Expand our possible output symbols also with a dummy-dimension so it fits
    % with receivedExpanded for matrix calculations. Will have dimensions:
    % 1 x nr of states x nr of inputs x bpsymb
    outputsExpanded = reshape(possibleOutputsMatrix, 1, size(possibleOutputsMatrix,1), size(possibleOutputsMatrix,2), bpsymb);
    
    % Mod-2 sum all received symbols with the possible output symbols. Leaves
    % ones wherever they are different. Sum up all ones across this 4:th
    % dimension. The 4:th dimension, bpsymb, is the number of bits being
    % processed at once for decoding. So for example for a rate 1/2 code the
    % 4:th dimension is only bpsymb = 2. 
    % This leaves us a matrix of dimension:
    % nr of steps x nr of states x nr of inputs
    % Which contains the Hamming distances between what was received and what
    % could have been output from the encoder. 
    hammingDistances = sum(mod(receivedExpanded + outputsExpanded, 2), 4);
    
    
    % Initialize accumulated path metrics
    pathMetrics = inf(numStates, numSteps + 1);
    pathMetrics(1, 1) = 0; 
    
    % Initialize the survivors matrix
    survivors = zeros(numStates, numSteps);
    
    for step = 1:numSteps
            currentMetrics = pathMetrics(:, step); % Metrics at the current step
            nextMetrics = inf(numStates, 1); % Initialize next step metric
    
            for state = 1:numStates
                for input = 1:numInputs
                    nextState = trellis.nextStates(state, input) + 1; % The +1 is to comply with matlab indexing
                    hammingCost = hammingDistances(step, state, input);
                    candidateMetric = currentMetrics(state) + hammingCost;
    
                    % Update path metrics and survivors
                    if candidateMetric < nextMetrics(nextState)
                        nextMetrics(nextState) = candidateMetric;
                        survivors(nextState, step) = state - 1; % The -1 is to switch from index to "state number"
                    end
                end
            end
    
            pathMetrics(:, step + 1) = nextMetrics; % Update metrics
    end
    
    % Traceback to find the best path
    decodedBits = zeros(1, numSteps);
    [~, finalState] = min(pathMetrics(:, end)); % Start from the state with the lowest metric
    for step = numSteps:-1:1
        previousState = survivors(finalState, step);
        inputBit = find(trellis.nextStates(previousState + 1, :) == (finalState - 1)) - 1;
        decodedBits(step) = inputBit; % Recover the input bit
        finalState = previousState + 1; % Move to the previous state
    end
end

function decodedBits = viterbi_soft(receivedSymbols, trellis, constellation)
    % Fetch initial parameters
    bpsymb = log2(length(constellation));
    rate = log2(trellis.numOutputSymbols)/bpsymb;
    numSteps = length(receivedSymbols)/rate;
    numStates = trellis.numStates;
    numInputs = trellis.numInputSymbols;
    
    if rate>1
        outputsBinary = de2bi(trellis.outputs, rate, 'left-msb');
    else
        outputsBinary = trellis.outputs;
    end
    possibleOutputsMatrix = reshape(outputsBinary, numStates, numInputs, rate);
    % Convert trellis outputs to constellation points
    possibleSymbols = constellation(possibleOutputsMatrix + 1);
    
    % Expand received symbols for vectorized distance computation
    receivedExpanded = reshape(receivedSymbols,rate,[]);
    if rate>1
        receivedExpanded = receivedExpanded';
    end
    receivedExpanded = reshape(receivedExpanded, numSteps, 1, 1, rate);
    
    outputsExpanded = reshape(possibleSymbols, 1, size(possibleSymbols,1), size(possibleSymbols,2),rate);
    
    % Compute Euclidean distances in the complex plane
    euclideanDistances = sum(abs(receivedExpanded - outputsExpanded).^2,4);
    
    % Initialize accumulated path metrics
    pathMetrics = inf(numStates, numSteps + 1);
    pathMetrics(1, 1) = 0; 
    
    % Initialize the survivors matrix
    survivors = zeros(numStates, numSteps);
    
    for step = 1:numSteps
            currentMetrics = pathMetrics(:, step); % Metrics at the current step
            nextMetrics = inf(numStates, 1); % Initialize next step metric
    
            for state = 1:numStates
                for input = 1:numInputs
                    nextState = trellis.nextStates(state, input) + 1; % The +1 is to comply with matlab indexing
                    % Only difference between hard-input and soft-input is that
                    % we use euclidean distance as a cost metric.
                    euclideanCost = euclideanDistances(step, state, input);
                    candidateMetric = currentMetrics(state) + euclideanCost;
    
                    % Update path metrics and survivors
                    if candidateMetric < nextMetrics(nextState)
                        nextMetrics(nextState) = candidateMetric;
                        survivors(nextState, step) = state - 1; % The -1 is to switch from index to "state number"
                    end
                end
            end
    
            pathMetrics(:, step + 1) = nextMetrics; % Update metrics
    end
    
    % Traceback to find the best path
    numInputBits = log2(trellis.numInputSymbols);
    decodedBits = zeros(1,numInputBits*numSteps);
    [~, finalState] = min(pathMetrics(:, end)); % Start from the state with the lowest metric
    for step = numSteps:-1:1
        previousState = survivors(finalState, step);
        inputBit_dec = find(trellis.nextStates(previousState + 1, :) == (finalState - 1)) - 1;
        inputBits = dec2bin(inputBit_dec, numInputBits) - '0'; % Convert from decimal to binary
        decodedBits(step*numInputBits-numInputBits+1:1:step*numInputBits) = inputBits; % Prepend the input bit(s)
        finalState = previousState + 1; % Move to the previous state
    end
end

% ======================================================================= %
% Simulation Chain
% ======================================================================= %
BER = zeros(1, length(EbN0)); % pre-allocate a vector for BER results

for i = 1:length(EbN0)
    totErr = 0;  % Number of errors observed
    num = 0; % Number of bits processed

    while((totErr < maxNumErrs) && (num < maxNum))
        % =====================================================================
        % % Begin processing one block of information
        % =====================================================================
        % % [SRC] generate N information bits

        info_bits = randi([0,1],1,N);

        % [ENC] convolutional encoder
        if encoder_option == 0
            code_rate = 1;
            tx_bits = info_bits; % No encoding

        elseif encoder_option == 1
            code_rate = 1/2;
            g1 = [1,0,1];
            g2 = [1,1,1];

            tx_bits = cnEncoder(info_bits,g1,g2); % Encoding using eps_1
            % The notation for the generator polynomials is that the
            % leftmost bit is D^0.

            trellis = create_trellis(3, [g1;g2]);
        elseif encoder_option == 2
            code_rate = 1/2;
            g1 = [1,0,1,1,1];
            g2 = [1,0,1,1,0];

            tx_bits = cnEncoder(info_bits,g1,g2);
            trellis = create_trellis(5, [g1;g2]);
        elseif encoder_option == 3
            code_rate = 1/2;
            g1 = [1,0,0,1,1];
            g2 = [1,1,0,1,1];

            tx_bits = cnEncoder(info_bits,g1,g2);
            trellis = create_trellis(5, [g1;g2]);

        elseif encoder_option == 4
            code_rate = 2/3;
            tx_bits = eps4(info_bits);
            
            trellis.numInputSymbols = 4;
            trellis.numOutputSymbols = 8;
            trellis.numStates = 8;

            trellis.nextStates = [0	2 1	3;
            4	6	5	7;
            1	3	0	2;
            5	7	4	6;
            2	0	3	1;
            6	4	7	5;
            3	1	2	0;
            7	5	6	4];

            trellis.outputs = [0 1 2 3;
            4	5	6	7;
            0	1	2	3;
            4	5	6	7;
            0	1	2	3;
            4	5	6	7;
            0	1	2	3;
            4	5	6	7];
         end

        % [MOD] symbol mapper
        switch constellation_option
            case 'QPSK'
                const = [(1 + 1i), (1 - 1i), (-1 + 1i), (-1 - 1i)]*(1/sqrt(2));
                bpsymb = 2;

            case 'BPSK'
                const = [1, -1];
                bpsymb = 1;

            case 'AMPM'
                % Check if the input length is a multiple of 3
                if mod(length(tx_bits), 3) ~= 0
                    tx_bits = [tx_bits, zeros(1, 3 - mod(length(tx_bits), 3))]; % if not, pad with zeros
                end
                % 000 -> (1,-1), 001 -> (-3,3), 010 -> (1,3), 011 -> (-3,1), 
                % 100 -> (3,-3), 101 -> (-1,1), 110 -> (3,1), 111 -> (-1,-3), 
                const = [(1 - 1i), (-3 + 3i), (1 + 3i), (-3 - 1i), (3 - 3i), (-1 + 1i), (3 + 1i), (-1 - 3i)] / sqrt(10);
                
                bpsymb = 3;
        end
        
        % Reshape and "send" the bits in packages of bpsymb
        m = reshape(tx_bits, bpsymb, [])';
        m_idx = bi2de(m, 'left-msb')'+1;
        x = const(m_idx);

        % [CHA] add Gaussian noise
        %signal_power = mean(abs(x).^2);
        E_s = mean(abs(const).^2);
        EbN0_linear = 10^(EbN0(i) / 10);
        noise_power = E_s / (code_rate * bpsymb * EbN0_linear);

        N_noise = length(x);
        n = sqrt(noise_power / 2) * (randn(1,N_noise) + 1i*randn(1,N_noise));

        y = x + n;

        if hard_decoding == true
            % [HR] Hard Receiver
             
            % Compute Euclidean distance to each constellation point
            dist = abs(repmat(y.', 1, length(const)) - repmat(const, length(y), 1));
            % Get the closest constellation points
            [~, idx] = min(dist, [], 2);

            % Convert symbol indices back to binary
            m_received = de2bi(idx-1, bpsymb, 'left-msb');
            
            received_bits = reshape(m_received.', [], 1);
            if encoder_option ~= 0
                decoded_bits = viterbi(received_bits, trellis);
            end
        elseif hard_decoding == false
            % [SR] Soft Receiver
            decoded_bits = viterbi_soft(y, trellis, const);
        end

        % =====================================================================
        % % End processing one block of information
        % =====================================================================
        % %

        if encoder_option == 0
            BitErrs = sum(xor(received_bits',info_bits));
        else
            BitErrs = sum(xor(decoded_bits(1:N), info_bits));
        end
        totErr = totErr + BitErrs;
        num = num + N;

        disp(['+++ ' num2str(totErr) '/' num2str(maxNumErrs) ' errors. '...
            num2str(num) '/' num2str(maxNum) ' bits. Projected error rate = '...
            num2str(totErr/num, '%10.1e') '. +++']);

    end
    BER(i) = totErr/num;
end

%save('./Task 3/eps4_AMPM.mat',"BER");
EbN0_linear = 10.^(EbN0/10);
% Calculate the upper bound
N=20;
spect = distspec(trellis, N);
dfree = spect.dfree;
weights = spect.weight;
Pb_soft_bound = zeros(size(EbN0));
for i = 1:length(EbN0)
    d_values = (dfree:dfree+N-1);
    Q_values = qfunc(sqrt(2 .* d_values .* code_rate .* EbN0_linear(i)));
    Pb_soft_bound(i) = sum(weights .* Q_values);
end
semilogy(EbN0, Pb_soft_bound, '-*', 'DisplayName', 'Soft Decoder Upper Bound');
hold on;
% Calculate the theoretical BER for the uncoded case.
BER_theory = qfunc(sqrt(bpsymb.*EbN0_linear));

semilogy(EbN0, BER_theory, '--', 'DisplayName', 'Uncoded QPSK (theoretical)');

semilogy(EbN0, BER, '-o', DisplayName=join(['Encoder', num2str(encoder_option)]));

% Plot formatting
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
legend show;
ylim([1e-4 1]);
title(join(['BER vs E_b/N_0 using ', constellation_option]));
hold off;
% ======================================================================= %
% End
% ======================================================================= %