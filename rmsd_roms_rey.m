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

% Abre apenas um arquivo para pegar lon e lat de Rey
rey      = netcdf.open('~/ANALISES/REYNOLDS_SST/sst_REYNOLDS_19880101.nc','NC_NOWRITE');    
it       = netcdf.inqVarID(rey,'LONN303_N40');
lon_rey  = netcdf.getVar(rey,it,'double');
it       = netcdf.inqVarID(rey,'LAT101_360');
lat_rey  = netcdf.getVar(rey,it,'double');

% Cria o grid para Reynolds
[Xlon,Ylat] = meshgrid(lon_rey,lat_rey);                                    

% Cria um contador e matrizes para erro medio quadratico
n0 = 0; 
rmsd_sst4 = zeros(size(Xlon));
rmsd_sst5 = zeros(size(Xlon));
rmsd_sst6 = zeros(size(Xlon));

% Loop Anual
for y=1990:1991
    
    % Loop mensal 
    for m = 1:12      

        n1 = 0;
        mean_rey = zeros(size(Xlon));
        
        % Loop diario para calculo da media mensa Reynolds
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
        
        % Calculo da media mensal de sst-Rey
        mean_rey = mean_rey/n1;

        % Lê variaveis comuns a todos os experimentos
        gfile     = netcdf.open('~/ROMS/GRD_files/sao12_LEILANE_grd.nc','NC_NOWRITE');
        it_r      = netcdf.inqVarID(gfile,'lon_rho');
        lon_rho   = netcdf.getVar(gfile,it_r);
        it_r      = netcdf.inqVarID(gfile,'lat_rho');
        lat_rho   = netcdf.getVar(gfile,it_r);

        % Abre resultados do ROMS para cada Experimento
        avg4       = netcdf.open(strcat('~/RESULTADOS/004lei/004lei_SST_sao12tbc_LEILANE_avg_',num2str(y),sprintf('%02g',m),'16.nc'),'NC_NOWRITE');
        avg5       = netcdf.open(strcat('~/RESULTADOS/005lei/005lei_SST_sao12tbc_LEILANE_avg_',num2str(y),sprintf('%02g',m),'16.nc'),'NC_NOWRITE');
        avg6       = netcdf.open(strcat('~/RESULTADOS/006lei/006lei_SST_sao12tbc_LEILANE_avg_',num2str(y),sprintf('%02g',m),'16.nc'),'NC_NOWRITE');
        
        it_r4      = netcdf.inqVarID(avg4,'temp');
        temp_roms4 = netcdf.getVar(avg4,it_r4,'double');
        it_r5      = netcdf.inqVarID(avg5,'temp');
        temp_roms5 = netcdf.getVar(avg5,it_r5,'double');
        it_r6      = netcdf.inqVarID(avg6,'temp');
        temp_roms6 = netcdf.getVar(avg6,it_r6,'double');
%         temp_roms = netcdf.getVar(avg,it_r,[0 0 49 0],[656 732 1 1],'double');
        
        % Cria mascara nas variáveis
        i4             = find(temp_roms4>=100);
        temp_roms4(i4) = nan;
        i5             = find(temp_roms5>=100);
        temp_roms5(i5) = nan;
        i6             = find(temp_roms6>=100);
        temp_roms6(i6) = nan;

        % Interpola o ROMS para os pontos de grade Reynolds
        sst_roms4 = griddata(lon_rho,lat_rho,temp_roms4,Xlon,Ylat,'cubic');
        sst_roms5 = griddata(lon_rho,lat_rho,temp_roms5,Xlon,Ylat,'cubic');
        sst_roms6 = griddata(lon_rho,lat_rho,temp_roms6,Xlon,Ylat,'cubic');
 
        % Calcula o Erro medio quadratico para todos os experimentos
        n0             = n0+1;
        rmsd_sst4      = rmsd_sst4 + (sst_roms4 - mean_rey).^2;
        rmsd_sstT4(n0) = nanrmse(sst_roms4,mean_rey);
        rmsd_sst5      = rmsd_sst5 + (sst_roms5 - mean_rey).^2;
        rmsd_sstT5(n0) = nanrmse(sst_roms5,mean_rey);
        rmsd_sst6      = rmsd_sst6 + (sst_roms6 - mean_rey).^2;
        rmsd_sstT6(n0) = nanrmse(sst_roms6,mean_rey);
        
        % Apaga a media sst-Rey a SST dos experimentos
        clear mean_rey sst_roms4 sst_roms5 sst_roms6
        disp(strcat(num2str(y),sprintf('%02g',m),'--OK!'))
        
        % Cria um vetor de tempo
        time(n0) = datenum(strcat(sprintf('%02g',m),'-',num2str(y)),'mm-yyyy');
    end
end

%%
% Plotando

% Termina o calculo de RMSD de acordo com o numero de meses (n)
rmsd_sst4 = sqrt(rmsd_sst4/n0);
rmsd_sst5 = sqrt(rmsd_sst5/n0);
rmsd_sst6 = sqrt(rmsd_sst6/n0);

figure(1)
contourf(Xlon,Ylat,rmsd_sst4); shading flat; colorbar; grid on;
set(gca, 'YLim', [min(min(lat_rho)), max(max(lat_rho))]);
set(gca, 'XLim', [min(min(lon_rho)), max(max(lon_rho))]);
set(gca, 'CLim', [rmsd_sst_min_CCSM4, rmsd_sst_max_CCSM4]);
title('Exp: 004lei - CC: Dani');
cblabel('Erro Médio Quadrático');
saveas(gcf,['~/RESULTADOS/COMPARACOES/004lei_005lei_006lei/004lei.png'])

figure(2)
contourf(Xlon,Ylat,rmsd_sst5); shading flat; colorbar; grid on;
set(gca, 'YLim', [min(min(lat_rho)), max(max(lat_rho))]);
set(gca, 'XLim', [min(min(lon_rho)), max(max(lon_rho))]);
set(gca, 'CLim', [rmsd_sst_min_CCSM4, rmsd_sst_max_CCSM4]);
title('Exp: 005lei - CC: Leo');
cblabel('Erro Médio Quadrático');
saveas(gcf,['~/RESULTADOS/COMPARACOES/004lei_005lei_006lei/005lei.png'])

figure(3)
contourf(Xlon,Ylat,rmsd_sst6); shading flat; colorbar; grid on;
set(gca, 'YLim', [min(min(lat_rho)), max(max(lat_rho))]);
set(gca, 'XLim', [min(min(lon_rho)), max(max(lon_rho))]);
set(gca, 'CLim', [rmsd_sst_min_CCSM4, rmsd_sst_max_CCSM4]);
title('Exp: 006lei - CC: Kate');
cblabel('Erro Médio Quadrático');
saveas(gcf,['~/RESULTADOS/COMPARACOES/004lei_005lei_006lei/006lei.png'])

figure(4)
plot(time_CCSM4,rmsd_sstT_CCSM4,'m');hold on;
plot(time,rmsd_sstT4,'k');
plot(time,rmsd_sstT5,'r');
plot(time,rmsd_sstT6,'b');
datetick('x','mmmyy');
set(gca, 'XLim', [min(time_CCSM4), max(time_CCSM4)]);
grid on;
ylabel('RMSD');
xlabel('Data da Simulação (mmmyy)')
title('Erro Médio Quadrático');
legend('CCSM4','CC: Dani','CC: Leo','CC: Kate')
saveas(gcf,['~/RESULTADOS/COMPARACOES/004lei_005lei_006lei/RMSD.png'])
