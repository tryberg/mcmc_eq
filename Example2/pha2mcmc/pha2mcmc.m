function pha2mcmc(phaseFile,staFile,vlat,vlon,outDir)

%{

convert HypoDD phase data format for use in mcmc_eq codes, 
J. Pesicek: July 2020

this file cannot have any unused stations and must start numbering at 0,
see mcmc_eq manual

https://github.com/tryberg/mcmc_eq

phaseFile:  hypoDD phase data file
staFile:    hypoDD station.dat file
vlat, vlon: coordinates for cartesian transformation

requires Mapping Toolbox for geodetic2enu

%}

staNames = getStaNamesFromPhaFile(phaseFile);
% sd = readtable(staFile);

fid = fopen(staFile);
C = textscan(fid,'%s %f %f %d\n');
lon = C{3}; lat = C{2};
Elevation = C{4};
dep = -double(Elevation)/1000;
sname = C{1};
sd = table(sname,lon,lat,Elevation,dep, 'VariableNames', {'Station','Longitude','Latitude','Elevation','Depth'});

[~,ia,~] = intersect(sd.Station,staNames);
[C,IA] = setdiff(sd.Station,staNames);

if ~isempty(IA)
    warning('removing unused stations')
    disp(C')
end
sd = sd(ia,:);

npw = zeros(numel(sd.Station),5);
nsw = npw;

%% grid origin
E = wgs84Ellipsoid('kilometer');
% sz = -sd.Elevation./1000;
[xEast,yNorth,~] = geodetic2enu(sd.Latitude,sd.Longitude,0,vlat,vlon,0,E);

%%
mcmcPhaFile=fullfile(outDir,'picks.mcmc');   
mcmcEqFile=fullfile(outDir,'quakes.dat');   
mcmcStaFile=fullfile(outDir,'stations.dat');  

%%
fid1=fopen(phaseFile,'r');
fid2=fopen(mcmcPhaFile,'w');
fid3=fopen(mcmcEqFile,'w');

%%
e1=0;
ids = [];
if fid1 ~= -1
    line = fgetl(fid1);
    while 1

        if ~ischar(line)
            break;
        end

        if line(1) == '#'
            %             disp('  !!  EVENT LINE  !!')
            e1 = e1 + 1;

            [yr,mo,dy,hr,mn,sc,lat,lon,dep,~,~,~,~,~] = readEventLine(line);

            line = fgetl(fid1);

            phaCt = 0; pCt = 0; sCt = 0; clear pha

            % need to instead read in all phase lines into matrix
            while line(1) ~= '#' & line ~= -1

                is = find(strcmpi(strtrim(line(1:6)),sd.Station), 1);

                if ~isempty(is)

                    phaCt = phaCt + 1;

                    [sta,tt,wgt,phas] = readPhaLine(line);

                    pha(phaCt,1) = {sta};     %sta
                    pha(phaCt,2) = {num2str(tt)};    %tt
                    pha(phaCt,3) = {num2str(wgt)};   %wgt
                    pha(phaCt,4) = {phas};   %pha

                    if strcmp(pha(phaCt,4),'P') || strcmp(pha(phaCt,4),'Pg')
                        pCt = pCt + 1;
                    elseif strcmp(pha(phaCt,4),'S') || strcmp(pha(phaCt,4),'Sg')
                        sCt = sCt + 1;
                    else
                        error(['FATAL: Bad phase: ',pha{phaCt,4}])
                    end
                else
                    error(['phase at missing station: ',strtrim(line(1:6))])
                end

                line = fgetl(fid1);

            end

            % sort phases by P, then S
            if phaCt > 0
                pha = sortrows(pha,[4,1]);% this may be needed by mcmc
            end
            %% mcmc format
            fecha = sprintf('%04d%02d%02d%02d%02d%05.2f',yr,mo,dy,hr,mn,sc);
            reftime=datetime(fecha,'InputFormat','yyyyMMddHHmmss.SSS');
            reftime = char(string(reftime,'yyyyMMddHHmmss.SSS'));

            [x(e1),y(e1),~] = geodetic2enu(lat,lon,0,vlat,vlon,0,E);

            zi = dep; % z(e1)
            idi = e1-1;
            ids = [ids; idi];
            fprintf(fid3,'%3d %8.3f %8.3f %8.3f %s %8.3f\n',idi,x(e1),y(e1),zi,reftime,0);

            fprintf(fid2,'# %d %d %d %s\n',...
                idi,pCt,sCt,reftime);

            for i=1:phaCt %

                wgt = str2double(pha{i,3});
                qual = pickWeight2quality(wgt);
                if qual > 3
                    warning('bad pick qual value > 3, setting to 3')
                    qual = 3;
                end

                sii = find(strcmpi(sd.Station,pha{i,1}));
                if numel(sii) > 1
                    warning(['duplicate stations, should not occur, fix it: ',char(sd.Station(sii(1)))])
                    sii = sii(1); %NOTE: fix later?
                end

                zi = sd.Depth(sii); %
                %                 zi = zUp(sii);
                fprintf(fid2,'%4s %03d %1s %8.3f %8.3f %8.3f %8.3f %d\n',...
                    pha{i,1},sii-1,pha{i,4},xEast(sii),yNorth(sii),zi,str2double(pha{i,2}),qual);

                j = find(strcmp(pha(i,1),sd.Station));
                if strcmp(pha(i,4),'P')
                    npw(j,qual+1) = npw(j,qual+1) + 1 ;
                else
                    nsw(j,qual+1) = nsw(j,qual+1) + 1 ;
                end

            end

        end

    end
end
fclose(fid1);
fclose(fid2);
fclose(fid3);
ls(mcmcPhaFile)
ls(mcmcEqFile)

if max(ids) ~= numel(ids)-1
    error('FATAL')
end
%% qual hit count
I = sum(npw,2) == 0;
if any(I)
    disp(sd(I,:).Station')
    error('station with zero P picks')
end
I = sum(nsw,2) == 0;
if any(I)
    disp(sd(I,:).Station')
    warning('station with zero S picks')
end

%% stations in cartesian
fid4=fopen(mcmcStaFile,'w'); %only needed for synthetics and disp_error script

for i=1:numel(xEast)
    % station number    x   y   z   Pcorr   Scorr (added later)
    fprintf(fid4,'%3d %8.3f %8.3f %8.3f %8.3f %8.3f\n',i-1,xEast(i),yNorth(i),sd.Depth(i),0,0);
end
fclose(fid4);
ls(mcmcStaFile)

end
%%
function staNames = getStaNamesFromPhaFile(phaseFile)

% Read the file into a cell array of strings
fileContents = fileread(phaseFile);
lines = strsplit(fileContents, '\n');

% Filter out lines that start with '#'
filteredLines = lines(~startsWith(lines, '#'));

% Extract the first column from the remaining lines
firstColumn = cellfun(@(x) strtok(x), filteredLines, 'UniformOutput', false);

% Sort and remove duplicates
uniqueSorted = unique(sort(firstColumn));

% Display the result
% disp(uniqueSorted);

staNames = uniqueSorted(2:end)'; % remove first empty cell, not sure why its there

end
%%
function [yr,mo,dy,hr,mn,sc,lat,lon,dep,mag,eh,ez,rms,id] = readEventLine(line)

C = textscan(line,'%s %d %d %d %d %d %f %f %f %f %f %f %f %f %s');
yr = C{2}; mo = C{3}; dy = C{4}; hr = C{5}; mn = C{6}; sc = C{7};
lat = C{8}; lon = C{9}; dep = C{10}; mag = C{11};
eh = C{12}; ez = C{13}; rms = C{14}; id = C{15};

end
%%
function [sta,tt,wgt,phas] = readPhaLine(line)

C = textscan(line,'%s %f %f %c');
sta = char(C{1});
tt = C{2}; phas = C{4}; wgt = C{3};

if numel(phas) > 2
    error('bad line read, try again')
end

end
%%
function  qual = pickWeight2quality(wgt)
%
if      wgt <= 1.0 && wgt > 0.5
    qual = 0;
elseif  wgt <= 0.5 && wgt > 0.2
    qual = 1;
elseif  wgt <= 0.2 && wgt > 0.1
    qual = 2;
elseif  wgt <= 0.1 && wgt > 0.05
    qual = 3;
elseif  wgt < 0 % hypoDD neg flag to keep regardless
    warning(['neg pick wgt ',num2str(wgt),' set to qual 0: flagged or error?'])
    qual = 0;
elseif wgt == 0
    qual = 4;
else
    qual = 4;
    disp(['pick wgt ',num2str(wgt),' set to qual 4 (bad for mcmc)'])
end

end