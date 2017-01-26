%~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~%
%                                                                 %
%           Projeto de Pesquisa em Mudancas Climaticas            %
%            Instituto Nacional de Pesquisas Espaciais            %
%              Ocª.Leilane Gonçalves dos Passos                   %
%                         2016-2017                               %
%                                                                 %
%~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~%
%       PROGRAMA PARA CALCULAR O RMSD ROMS/REYNOLDS_SST           %
%                   DATA: 19/12/2016                              %
%~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~%
%                                                                 %
% authors: Leilane/Leonardo                                       %
%~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~%
%%
clear all; close all; warning off; clc;

load '~/ROTINAS/MATLAB/RMSD/rmsd_min_max_CCSM4.mat'
load '~/ROTINAS/MATLAB/RMSD/rmsdT_CCSM4.mat'

% Definindo numero de dias em cada mês para fazer a media dos campos Reynolds
MM  = [31 28 31 30 31 30 31 31 30 31 30 31];  

%-------------------------------------------------------------------------%
%Carrega lat e lon do arquivo REYNOLDS
%-------------------------------------------------------------------------%

rey      = netcdf.open('~/ANALISES/REYNOLDS_SST/sst_REYNOLDS_19880101.nc','NC_NOWRITE');    
it       = netcdf.inqVarID(rey,'LONN303_N40');
lon_rey  = netcdf.getVar(rey,it,'double');
it       = netcdf.inqVarID(rey,'LAT101_360');
lat_rey  = netcdf.getVar(rey,it,'double');

% Cria o grid para Reynolds
[Xlon_rey,Ylat_rey] = meshgrid(lon_rey,lat_rey);   

%-------------------------------------------------------------------------%
% Lê variaveis de lat e lon da minha região do GRIDfile ROMS
%-------------------------------------------------------------------------%

gfile     = netcdf.open('~/ROMS/GRD_files/sao12_LEILANE_grd.nc','NC_NOWRITE');
it_r      = netcdf.inqVarID(gfile,'lon_rho');
lon_rho   = netcdf.getVar(gfile,it_r);
it_r      = netcdf.inqVarID(gfile,'lat_rho');
lat_rho   = netcdf.getVar(gfile,it_r);

%-------------------------------------------------------------------------%
%Carrega lat e lon e tempo do arquivo do CCSM4
%-------------------------------------------------------------------------%

nOGCM1 = netcdf.open('~/RESULTADOS/CMIP5/CCSM4/historical/ocean/raw/thetao/thetao_Omon_CCSM4_historical_r6i1p1_198001-198912.nc','NC_NOWRITE');
nOGCM2 = netcdf.open('~/RESULTADOS/CMIP5/CCSM4/historical/ocean/raw/thetao/thetao_Omon_CCSM4_historical_r6i1p1_199001-199912.nc','NC_NOWRITE');

%  Get the lat and lon arrays
varOGCM = netcdf.inqVarID(nOGCM1,'lat');
lat     = netcdf.getVar(nOGCM1,varOGCM,'double');
varOGCM = netcdf.inqVarID(nOGCM1,'lon');
lon     = netcdf.getVar(nOGCM1,varOGCM,'double');
varOGCM = netcdf.inqVarID(nOGCM1,'time');
tOGCM1   = netcdf.getVar(nOGCM1,varOGCM,'double');
varOGCM = netcdf.inqVarID(nOGCM2,'time');
tOGCM2   = netcdf.getVar(nOGCM2,varOGCM,'double');

% varOGCM = netcdf.inqVarID(nOGCM,'lev');
% lev   = netcdf.getVar(nOGCM,varOGCM,'double');

% Convert time in noLeapDataVec
time_OGCM1 = noLeapDateVec(tOGCM1);
time_OGCM2 = noLeapDateVec(tOGCM2);

% Convert longitude from (0-360) to (-180 to 180)
lon1=rem((lon+180),360)-180;

% Create big Xlon and YLAT arrays
[Xlat,Ylon] = meshgrid(lat(1,:),lon1(:,1));

% Encontra os inds do CCSM4 do grid correspondente a minha gridfile
[ila_min,dist] = near(Xlat(1,:),min(min(lat_rho))); %Pega um limite um pouco maior
[ila_max,dist] = near(Xlat(1,:),max(max(lat_rho)));

fc = find(Ylon(:,1)>=min(min(lon_rho)) & Ylon(:,1)<=max(max(lon_rho)));

%-------------------------------------------------------------------------%
% Separa a regiao de interesse do CCSM4
%-------------------------------------------------------------------------%

lat_c  = Xlat(fc,ila_min:ila_max);
lon_c  = Ylon(fc,ila_min:ila_max);

% Cria um contador e matrizes para erro medio quadratico
n0 = 0; 
rmsd_sst = zeros(size(Xlon_rey));

%%
%=========================================================================%
% LOOP YEARLY
%=========================================================================%
for y=1987:1991
    
    %=====================================================================%
    % LOOP MONTHLY 
    %=====================================================================%    
    for m = 1:12      

        n1 = 0;
        mean_rey = zeros(size(Xlon_rey));
        
        %=================================================================%
        % LOOP DAYLY para calculo da media mensal Reynolds
        %=================================================================%        
        for d = 1:MM(m)
            n1  = n1 + 1;
            
            % Abre o arquivo do mes e ano desejado
            rey     = netcdf.open(strcat('~/ANALISES/REYNOLDS_SST/sst_REYNOLDS_',num2str(y),sprintf('%02g',m),sprintf('%02g',d),'.nc'),'NC_NOWRITE');
            it      = netcdf.inqVarID(rey,'SST_REYNOLDS');
            sst_rey = netcdf.getVar(rey,it)';       
    
            % Cria uma mascara para o arquivo de Rey
            f          = find(sst_rey<=-100);
            sst_rey(f) = nan;
            
            % Faz um somatorio dos dados de sst para posterior calculo da
            % media
            mean_rey = mean_rey + sst_rey;
        end 
        %=================================================================%
        % END LOOP DAYLY 
        %=================================================================%
        
        % Calculo da media mensal de sst-Rey
        mean_rey = mean_rey/n1;
      
        %-----------------------------------------------------------------%
        %Carrega variavel do arquivo CCSM4
        %-----------------------------------------------------------------%
             
        % Get variable ID of the field variable requested, given its name.
        varOGCM1 = netcdf.inqVarID(nOGCM1,'thetao');
        varOGCM2 = netcdf.inqVarID(nOGCM2,'thetao');
        
        if y <= 1989
            % Find year indices for CCSM4
            iOGCM1 = find(time_OGCM1(:,1)==y & time_OGCM1(:,2)==m);
            vard = netcdf.getVar(nOGCM1,varOGCM1,[0 0 0 (iOGCM1-1)],[320 384 1 1],'double'); % Pega lat e lon da camada 1(superficie=5m) tempo (indice encontrado)
        else
            % Find year indices for CCSM4
            iOGCM2 = find(time_OGCM2(:,1)==y & time_OGCM2(:,2)==m);
            vard = netcdf.getVar(nOGCM2,varOGCM2,[0 0 0 (iOGCM2-1)],[320 384 1 1],'double'); % Pega lat e lon da camada 1(superficie=5m) tempo (indice encontrado)
        end
        
        % Find Land Mask
        ii = find(vard >= 1e+20);
        vard(ii) = nan;
        
        % Converted temp units from 'Kelvin' to 'Celsius';
        vard = vard - 273.15;       
      
        % Separa a variavel na regiao de interesse (GRIDfile)
        thetao = vard(fc,ila_min:ila_max);
        
        %-----------------------------------------------------------------%
        
        % INTERPOLA dados Reynolds para a grade do CCSM4 (para a região de interesse(GRIDfile) já separada)
%         mean_rey_interp = griddata(Xlon_rey,Ylat_rey,mean_rey,lon_c,lat_c,'cubic');
 
        % INTERPOLA dados Reynolds para a grade do CCSM4 (para a região de interesse(GRIDfile) já separada)
        thetao_interp = griddata(lon_c,lat_c,thetao,Xlon_rey,Ylat_rey,'cubic');
        
        % CALCULA o Erro medio quadratico para todos os experimentos
        n0             = n0+1;
        rmsd_sst      = rmsd_sst + (thetao_interp - mean_rey).^2;
        rmsd_sstT(n0) = nanrmse(thetao_interp,mean_rey);
        
        % Apaga a media sst-Rey a SST do CCSM4
        clear mean_rey thetao_interp
        
        disp(strcat(num2str(y),sprintf('%02g',m),'--OK!'))
        
        % Cria um vetor de tempo
        time(n0) = datenum(strcat(sprintf('%02g',m),'-',num2str(y)),'mm-yyyy');
    end
    %=====================================================================%    
    % END LOOP MONTHLY
    %=====================================================================%    
end
%=========================================================================%
% END LOOP YEARLY
%=========================================================================%
%%
% Plotando

% Termina o calculo de RMSD de acordo com o numero de meses (n)
rmsd_sst = sqrt(rmsd_sst/n0);

figure(1)
contourf(Xlon_rey,Ylat_rey,rmsd_sst); shading flat; colorbar; grid on;
set(gca, 'YLim', [min(min(lat_rho)), max(max(lat_rho))]);
set(gca, 'XLim', [min(min(lon_rho)), max(max(lon_rho))]);
set(gca, 'CLim', [rmsd_sst_min_CCSM4, rmsd_sst_max_CCSM4]);
title('RMSD CCSM4 em relação à OISST');
cblabel('Erro Médio Quadrático');
saveas(gcf,['~/RESULTADOS/COMPARACOES/CCSM4_OISST/RMSD_CCSM4_OISST.png'])

figure(2)
plot(time,rmsd_sstT,'k');hold on;
datetick('x','mmmyy');
set(gca, 'XLim', [min(time), max(time)]);
grid on;
ylabel('RMSD');
xlabel('Data da Simulação (mmmyy)')
title('Erro Médio Quadrático');
saveas(gcf,['~/RESULTADOS/COMPARACOES/CCSM4_OISST/RMSD_CCSM4_OISST_total.png'])
