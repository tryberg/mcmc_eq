
%in
phaseFile = 'Hood_syn.pha';
staFile = 'stations.dat';
vLon = -121.6950;
vLat = 45.3740;

%out
mcmcPhaFile = '../picks.mcmc';
mcmcSrcFile='../quakes.dat';   % mcmc src file

%%
pha2mcmc(phaseFile,staFile,mcmcPhaFile,mcmcSrcFile,vLat,vLon)