m = 300;
sampleFrequency = 1;
Ninitial = 10000;
Nfinal = Ninitial + ceil(10*30*24*60*60*sampleFrequency);
T_Neg = randn(Ninitial,1);
T_Pos = randn(Nfinal,1);
t_Pos = randn();
MP_PosNeg = randn(Nfinal-m+1,1);
MP_PosPos = randn(Nfinal-m+1,1);


tic;
[CP, MP_PosNeg, MP_PosPos, T_Pos] = ContrastProfileI(T_Neg, T_Pos, randn(), MP_PosNeg, MP_PosPos, m);
toc;

%%%Run the following to get the brute force time
% m = 300;
% sampleFrequency = 1;
% Ninitial = 10000;
% Nfinal = Ninitial + ceil(6*30*24*60*60*sampleFrequency);
% T_Neg = randn(Ninitial,1);
% T_Pos = randn(Nfinal,1);
% 
% tic;
% ContrastProfile(T_Neg, T_Pos, m);
% toc;