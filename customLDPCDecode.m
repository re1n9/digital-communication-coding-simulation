function [decoded_bits, success, iter_count] = stableLDPCDecode(llr, H, max_iterations)
    [M, N] = size(H);
    numCodewords = size(llr, 2);
    decoded_bits = zeros(N, numCodewords);
    success = false(1, numCodewords);
    iter_count = zeros(1, numCodewords);
    
    % Pre-compute sparse matrix structure
    [check_indices, var_indices] = find(H);
    num_edges = length(check_indices);
    nnz_H = nnz(H);
    
    % Create lookup tables for connected nodes
    var_to_check_table = cell(N, 1);
    check_to_var_table = cell(M, 1);
    for v = 1:N
        var_to_check_table{v} = find(H(:,v));
    end
    for c = 1:M
        check_to_var_table{c} = find(H(c,:));
    end
    
    % Process one codeword at a time
    for cw = 1:numCodewords
        % Get current codeword's LLRs
        current_llr = llr(:, cw);
        
        % Initialize messages
        var_to_check = zeros(M, N);
        check_to_var = zeros(M, N);
        
        % Add damping factor for stability
        alpha = 0.85;  % Damping factor between 0 and 1
        
        % Main iteration loop
        for iter = 1:max_iterations
            % Store previous messages for damping
            prev_check_to_var = check_to_var;
            
            % Step 1: Variable to check messages
            for v = 1:N
                connected_checks = var_to_check_table{v};
                if ~isempty(connected_checks)
                    % Calculate total LLR for this variable node
                    total_llr = current_llr(v);
                    for c = connected_checks'
                        total_llr = total_llr + check_to_var(c,v);
                    end
                    
                    % Update messages to each connected check node
                    for c = connected_checks'
                        % Exclude the incoming message from the target check node
                        var_to_check(c,v) = total_llr - check_to_var(c,v);
                    end
                end
            end
            
            % Step 2: Check to variable messages using offset min-sum
            offset = 0.3;  % Offset parameter
            new_check_to_var = zeros(M, N);
            
            for c = 1:M
                connected_vars = check_to_var_table{c};
                if ~isempty(connected_vars)
                    % Get messages for this check node
                    messages = var_to_check(c, connected_vars);
                    signs = sign(messages);
                    magnitudes = abs(messages);
                    
                    % Calculate product of signs
                    sign_product = prod(signs);
                    
                    for i = 1:length(connected_vars)
                        v = connected_vars(i);
                        
                        % Find minimum magnitude excluding current variable
                        other_magnitudes = magnitudes;
                        other_magnitudes(i) = inf;
                        min_magnitude = min(other_magnitudes);
                        
                        % Apply offset min-sum
                        if min_magnitude > offset
                            min_magnitude = min_magnitude - offset;
                        end
                        
                        % Compute final message
                        final_sign = sign_product * signs(i);
                        new_check_to_var(c,v) = final_sign * min_magnitude;
                    end
                end
            end
            
            % Apply damping
            check_to_var = alpha * new_check_to_var + (1 - alpha) * prev_check_to_var;
            
            % Clip extremely large values for numerical stability
            max_llr = 25;
            check_to_var = min(max(check_to_var, -max_llr), max_llr);
            
            % Make hard decision
            current_decoded = zeros(N, 1);
            for v = 1:N
                connected_checks = var_to_check_table{v};
                total_message = current_llr(v);
                if ~isempty(connected_checks)
                    total_message = total_message + sum(check_to_var(connected_checks,v));
                end
                current_decoded(v) = total_message < 0;
            end
            
            % Check syndrome
            syndrome = mod(H * current_decoded, 2);
            num_errors = sum(syndrome);
            fprintf('CW: %d Iter: %d syndrome: %d\n', cw, iter, num_errors);
            
            if num_errors == 0
                decoded_bits(:,cw) = current_decoded;
                success(cw) = true;
                iter_count(cw) = iter;
                break;
            end
            
            if iter == max_iterations
                decoded_bits(:,cw) = current_decoded;
                iter_count(cw) = max_iterations;
            end
        end
    end
    
    decoded_bits = decoded_bits(1:size(llr,1), :);
end