function [MFCC_Coef] = melfcc_test()
close all
clear all
clc
usedSR =  22050;
coeffNum = 12;
hopTime = 0.016;%0.020;
winTime = 0.032;%0.025;


[d, sr, duration] = readAudio(usedSR);

% Convert to MFCCs very close to those genrated by feacalc -sr 22050 -nyq 8000 -dith -hpf -opf htk -delta 0 -plp no -dom cep -com yes -frq mel -filt tri -win 32 -step 16 -cep 20
[MFCC_Coef,aspc] = melfcc(d*3.3752, sr, 'maxfreq', 10000, 'numcep', coeffNum, 'nbands', 22, 'fbtype', 'fcmel', 'dcttype', 1, 'usecmp', 1, 'wintime', winTime, 'hoptime', hopTime, 'preemph', 0, 'dither', 1);
% .. then convert the cepstra back to audio (same options)

drawPlot(MFCC_Coef,d, hopTime,duration);

end