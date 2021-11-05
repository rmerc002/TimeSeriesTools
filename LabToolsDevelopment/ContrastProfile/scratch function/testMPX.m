subLength = 4000;
noiseLvl = 0;
ts = [];
ts = [ts; rand(2*subLength,1)];
ts = [ts; sin(linspace(0,5*2*pi,subLength))' + noiseLvl*rand(subLength,1)];
ts = [ts; rand(2*subLength,1)];
ts = [ts; sin(linspace(0,5*2*pi,subLength))'+ noiseLvl*rand(subLength,1)];
ts = [ts; rand(2*subLength,1)];

ts2 = [];
ts2 = [ts2; rand(2*subLength,1)];
ts2 = [ts2; sin(linspace(0,5*2*pi,subLength))' + noiseLvl*rand(subLength,1)];
ts2 = [ts2; rand(2*subLength,1)];
ts2 = [ts2; sin(linspace(0,5*2*pi,subLength))'+ noiseLvl*rand(subLength,1)];
ts2 = [ts2; rand(2*subLength,1)];

[mp,mpi] = mpx_ryan(ts,subLength,subLength);
figure; plot(real(mp));
ylim([-0.5, 1.5]);
title(sprintf('subLength = %d',subLength)); 
