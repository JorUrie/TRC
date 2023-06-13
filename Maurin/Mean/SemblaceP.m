%% Limpiando pantalla, variables e imagenes hechas anteriormente
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
PhiCF = importdata('PhiCF-CM-P.csv'); 
PhiCD = importdata('PhiCD-CM-P.csv'); 
% Vector de las fechas
% Con solo uno es necesario, deben tener el mismo tamaño
Date = importdata('FechaCDCM.csv');
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
figure(1)
subplot(3,1,1); imagesc(SCD,[-1 1]); axis xy; axis tight; title('Semblance Convection-Diffusion Magnetic field'); ylabel('Wavelength');
colormap(jet(256));
colorbar
subplot(3,1,2); imagesc(SCF,[-1 1]); axis xy; axis tight; title('Semblance Force Field Magnetic field'); ylabel('Wavelength'); 
colormap(jet(256));
colorbar
subplot(3,1,3); imagesc(ST,[-1 1]); axis xy; axis tight; title('Diff. FF-CD'); ylabel('Wavelength'); 
colormap(jet(256));
colorbar
%% Haciendo grafica para ver los valores iguales a -1, 0, 1
% Diferencia entre Campo de Fuerza y Conveccion-Difusion
STm1 = find(ST<0);
STM1 = find(ST>0);
m1T = length(STm1);
M1T = length(STM1);
STT = m1T+M1T;
m1t = m1T/STT; M1t= M1T/STT;
TR = [m1t, M1t];
lblT = {'<0','>0'};
%% Haciendo las graficas
figure(2)
plot(1); pie(TR); title('Diff. FF-CD'); lgdT = legend(lblT, 'Location', 'bestoutside');
