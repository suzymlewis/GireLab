function [turb_label] = get_turb(errors);
%%%%%%%%%%%%%%% OUTPUTS ARRAY WHERE EACH INDEX RESPRESENTS ITS CORRESPONSING TRIAL NUMBER AND ITS ENTRY IS THE TURBULENCE LABEL
    % INPUT (optional)
        % errors: Trial numbers for each blocks if they are abnormal
            % 1 x 8 array
                % errors(1) : errors(end-1) are the start index for each block
                % errors(end) is the final trial number
                    % e.g. [ 1 11 16 21 26 31 37 42]
                        % 1:11- block 1, always medium turbulence
                        % 11:16 21:26 and 31:37 - block 2 4 6
                            % High/Low as indicated by user
                        % 16:21 26:31 and 37:42 - block 3 5 7
                            % Opposite 2 4 6 indicated label choice
    % OUTPUT
        % turb_label: m x 1 array
            % Indices correspond to trial number
                % e.g. turb_label(2) is the turbulence condition for trial 2
            % Conditions:
                % 0 : Low
                % 1 : Medium
                % 2 : High
%%%%%%%%%%%%%%%
if nargin == 0
    warning('Please note this function assumes there are no errors in user input unless optional argument "errors" is used')
    i = input('Are trials 10-15 high(1) or low(2): ');
    turb_label(1:10) = 1;
    if i == 1
        % HIGH-LOW 
        turb_label([11:15 21:25 31:35]) = 2;
        turb_label([16:20 26:30 36:40]) = 0;
    elseif i == 2
        % LOW-HIGH
        turb_label([11:15 21:25 31:35]) = 0;
        turb_label([16:20 26:30 36:40]) = 2;
    else
        error('Input not acceptable turbulence level')
    end
end
if nargin == 1
    if size(errors)~= [8 1] | size(errors)~= [1 8]
        error ('"Errors" input array needs to be 8 x 1 array.')
    end
    turb_label(errors(1):errors(2)) = 1;
    i = input('Are block 2 trials high(1) or low(2): ');
    if i == 1
        % HIGH-LOW
        turb_label([errors(2):errors(3) errors(4):errors(5) errors(6):errors(7)]) = 2;
        turb_label([errors(3):errors(4) errors(5):errors(6) errors(7):errors(8)]) = 0;
    elseif i == 2
        % LOW-HIGH
        turb_label([errors(2):errors(3) errors(4):errors(5) errors(6):errors(7)]) = 0;
        turb_label([errors(3):errors(4) errors(5):errors(6) errors(7):errors(8)]) = 2;
    else
        error('Input not acceptable turbulence level')
    end
end

