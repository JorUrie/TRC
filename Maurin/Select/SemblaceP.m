%% Limpiando comandos, variables e imagenes hechas anteriormente
clear
close
clc
%% Los datos fueron descargados en:
% Manchas solares: https://www.sidc.be/silso/datafiles
% Campo Magnético: http://wso.stanford.edu/#CoronalModels
Dok = readtable('Data-CM.csv', 'Range', 'C2:C784', 'ReadRowNames', false); % Data
%% Llamando los archivos csv creados por los programas en Python
% Manchas solares o Campo Magnético Py
% Aqui van las Phi que se han hecho en Python
PhiCF = importdata('PhiCF-CM-AATB.csv'); 
PhiCD = importdata('PhiCD-CM-AATB.csv'); 
% Vector de las fechas
% Con solo uno es necesario, deben tener el mismo tamaño
Date = importdata('FechaCD-CM.csv');
%% Conviertiendo de table a double 
% Solo es para los datos de Manchas y Campo magnético
Data1 = table2array(Dok);
%% Transponiendo vectores
Fecha = Date.';
PhiCFr = PhiCF';
PhiCDr = PhiCD';
Data = Data1.';
%% Manchas Solares o Campo Magnético
% Se hace la semblanza y se hace la diferencia
SCF = semblance(Fecha, Data, PhiCFr, 100); % Campo de Fuerza
SCD = semblance(Fecha, Data, PhiCDr, 100); % Campo de Fuerza
ST = SCF-SCD;
STi = SCD-SCF;
%% Haciendo las graficas
subplot(3,1,1); imagesc(SCD, [-1 1]); axis xy; axis tight; title('Semblance Convection-Diffusion AATB'); ylabel('Wavelength'); 
colormap(jet(256));
colorbar
subplot(3,1,2); imagesc(SCF,[-1 1]); axis xy; axis tight; title('Semblance Force Field AATB'); ylabel('Wavelength'); 
colormap(jet(256));
colorbar
subplot(3,1,3); imagesc(ST,[-1 1]); axis xy; axis tight; title('Diff. FF-CD'); ylabel('Wavelength'); 
colormap(jet(256));
colorbar
%%