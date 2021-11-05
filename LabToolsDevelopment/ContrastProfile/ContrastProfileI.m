function [CP, MP_PosNeg_NEW, MP_PosPos_NEW, T_Pos_NEW] = ContrastProfileI(T_Neg, T_Pos, t_Pos, MP_PosNeg, MP_PosPos, m)
    %%%Ryan Mercer
    %%%Demo Version
    %%%
    %%%Requirements: statistics and machine learning (zscore)
    %%%
    %%% Example 1 input
    % m = 1000;
    % sampleFrequency = 20;
    % Ninitial = 10000;
    % Nfinal = Ninitial + ceil(9*24*60*60*sampleFrequency);
    % T_Neg = randn(Ninitial,1);
    % T_Pos = randn(Nfinal,1);
    % t_Pos = randn();
    % MP_PosNeg = randn(Nfinal-m+1,1);
    % MP_PosPos = randn(Nfinal-m+1,1);

    %%% Example 2 input
%     m = 1000;
%     sampleFrequency = 1;
%     Ninitial = 10000;
%     Nfinal = Ninitial + ceil(10*30*24*60*60*sampleFrequency);
%     T_Neg = randn(Ninitial,1);
%     T_Pos = randn(Nfinal,1);
%     t_Pos = randn();
%     MP_PosNeg = randn(Nfinal-m+1,1);
%     MP_PosPos = randn(Nfinal-m+1,1);
    
    dataOrientation = 0; %0 for row, 1 for column
    %Change to row vector for internal use
    if size(T_Pos,1) == 1 
        %save the data orientation for matching output format to input
        %give orientation priority to positiveTS
        dataOrientation = 1; 
        T_Pos = T_Pos';
    end
    if size(T_Neg,1) == 1 
        %give orientation priority to positiveTS, do not save negative
        %orientation if different than positive
        T_Neg = T_Neg';
    end
    if size(MP_PosNeg,1) == 1 
        %give orientation priority to positiveTS, do not save negative
        %orientation if different than positive
        MP_PosNeg = MP_PosNeg';
    end
    if size(MP_PosPos,1) == 1 
        %give orientation priority to positiveTS, do not save negative
        %orientation if different than positive
        MP_PosPos = MP_PosPos';
    end
    
    
    
    T_Pos_NEW = [T_Pos; t_Pos]; %line 1, ContrastProfileI Algorithm 
    startIndex = length(T_Pos_NEW) - m + 1;
    endIndex = length(T_Pos_NEW);
    NEW = T_Pos_NEW(startIndex:endIndex); %line 2, ContrastProfileI Algorithm
    
    
    %%% Update MP_PosNeg, the matrix profile AB-join
    massK = pow2(nextpow2(sqrt(length(T_Neg))))*4;
    DP_PosNeg = [];
    while isempty(DP_PosNeg)
        DP_PosNeg = MASS_V4(T_Neg, NEW, massK);%line 3, ContrastProfileI Algorithm 
        massK = massK * 2;
    end
    MP_PosNeg_NEW = [MP_PosNeg; min(DP_PosNeg)];%line 4, ContrastProfileI Algorithm
    
    
    %%% Update MP_PosPos, the matrix profile self-join
    massK = pow2(nextpow2(sqrt(length(T_Pos)))+1);
    DP_PosPos = [];
    while isempty(DP_PosPos)
        DP_PosPos = MASS_V4(T_Pos, NEW, massK); %line 5, ContrastProfileI Algorithm
        massK = massK * 2;
    end
    DP_PosPos(length(MP_PosPos)+1:end) = []; %take self-join exclusion zone into account
    MP_PosPos = min(MP_PosPos, DP_PosPos); %line 6, ContrastProfileI Algorithm
    MP_PosPos_NEW = [MP_PosPos; min(DP_PosPos)]; %line 7, ContrastProfileI Algorithm
    

    %%% Contrast Profile
    CP = MP_PosNeg_NEW - MP_PosPos_NEW; %line 8, ContrastProfileI Algorithm
    
    %%% return in the same vector orientation received
    if dataOrientation == 1
       CP = CP';
       MP_PosNeg_NEW = MP_PosNeg_NEW';
       MP_PosPos_NEW = MP_PosPos_NEW';
       T_Pos_NEW = T_Pos_NEW';
    end
    
end