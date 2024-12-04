
%in
phaseFile = 'Hood_syn.pha';
staFile = 'stations.sta';
vLon = -121.6950;
vLat = 45.3740;

%out dir
outDir = '../';

%%
pha2mcmc(phaseFile,staFile,vLat,vLon,outDir)