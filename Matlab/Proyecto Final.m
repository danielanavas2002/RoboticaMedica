% ------------------------------------------------------------------- 
% PROYECTO FINAL
% ROBÓTICA MÉDICA - 2025
% 
% Aury Campaneros & Daniela Navas 
% ------------------------------------------------------------------
close all; clear; clc;
%% FASE 1 - CENTROIDE Y TRAYECTORIA EN SLICER
% 1. Carga de información de archivos para uso en la App "Medical Image Labeler"
filepath = '4.000000-Ax T1 MP SPGR-PRE GAD-85514';
dataFolder = fullfile(filepath);

% Extraer la información general de la imagen
medVol = medicalVolume(dataFolder);

% 2. Definición de Centro de Masa (centroide del tumor)
load('groundTruthMed.mat')

% Obtención de Voxeles que conforman el volumen definido en la App
% Obtenemos la tabla con la información de todas las "Lables" definidas en
% las imágenes. 
labelDefs = gTruthMed.LabelDefinitions;

% Obtenemos la dirección del archivo en donde se guarda la información de
% las "Labels" generadas
labelFilePath = gTruthMed.LabelData;

% Obtenemos el volumen del set de imágenes utilizando "niftiread" o
% "dicomreadVolume" según el formato de las imágenes que estemos analizando
labelVolume = niftiread(labelFilePath); 

% Especificamos el nombre de la "Label" que queremos analizar
labelName = 'Tumor';  

% Seleccionamos el ID asignado a la "Label" que deseamos analizar
% La instrucción "strcmp" sirve para realizar una comparación de texto
% entonces tomaremos en cuenta solo las etiquetas que correspondan al
% nombre de la región definida anteriormente 
labelID = labelDefs.PixelLabelID(strcmp(labelDefs.Name, labelName));

% Utilizando la función "find" haremos un arreglo con los datos de
% "labelVolume" que cumplan con el índice la región en cuestion, lo que nos
% proporcionará la localización (coordenadas intrínsicas) de los Voxels 
% correspndientes
RegionVoxels = find(labelVolume == labelID);

% Generación de Matriz de Voxels 
% Luego encontraremos la localización en el formato (fila, columna, slice)
% de cada uno de los Voxels encontrados en el inciso anterior. Esto se
% registra desde la perspectiva transversal.
[rows, cols, slices] = ind2sub(size(labelVolume), RegionVoxels);

% Finalmente guardaremos la información de todos los Voxels en una sola
% matriz
voxelLocation = ([rows, cols, slices]);

% Cálculo de centroide (centro de masa) de la región

% Para calcular el centroide sacaremos el promedio de cada una de las 
% coordenadas de los Voxels, para esto utilizaremos la siguiente notación

% rows = i, cols = j, slices = k
cent_i = mean(rows);
cent_j = mean(cols);
cent_k = mean(slices);

centroid_voxels = [cent_i, cent_j, cent_k];

% 3. Conversión a Coordenadas del Paciente (en mm)
% Para esta parte necesitamos leer las propiedades del Nifti o Dicom
% original, en donde se encuentra la trasformación utilizada al generar las
% imágenes. para esto debemos realizar lo siguiente: 

% 1. Leer las propiedades del set de imágenes original
dicom_prop = dicominfo('4.000000-Ax T1 MP SPGR-PRE GAD-85514\1-001.dcm');

% Definición de Variables extrayendo información del Archivo DICOM para 
% construir matriz de transformación
ori = double(dicom_prop.ImageOrientationPatient(:));
B = ori(1:3);
C = ori(4:6);
BxC = cross(B,C);

pos = dicom_prop.ImagePositionPatient(:);
pos = [-pos(1,1); -pos(2,1); pos(3,1)];

pix = double(dicom_prop.PixelSpacing(:));
delta_x = -pix(1,1);
delta_y = -pix(2,1);

delta_z = dicom_prop.SpacingBetweenSlices;

T_1 = [ B  C BxC pos;
        0  0  0   1 ]; 

T_2 = [delta_x    0       0       0   ;
          0    delta_y    0       0   ;
          0       0    delta_z    0   ;
          0       0       0       1   ];

T = T_1 * T_2;

matriz_trans_3d = affinetform3d(T);

% Crear el objeto de referencia espacial tridimencional, el cual se
% calcula a partir de las dimenciones originales (en voxels) de las 
% imágenes y la matriz de transoformación
dicomFiles = dir(fullfile(filepath, '*.dcm')); % lista de archivos DICOM
numSlices = length(dicomFiles);                 % número real de slices
imageSize = double([dicom_prop.Rows, dicom_prop.Columns, numSlices]);

ref3d = medicalref3d(imageSize, matriz_trans_3d);

% Aplicamos la función Intrinsic (Voxels) a World (Paciente)
% Para este paso debemos de tomar en cuenta que existen ciertas diferencias
% entre el sistema de coordinadas de Matlab con 3D Slicer y 
[x_paciente, y_paciente, z_paciente] = intrinsicToWorld(ref3d, cent_j, cent_i, cent_k); 
centroid_mm = [x_paciente, y_paciente, z_paciente];

% 4.Gráfica 3D
% Leer el archivo JSON
dataTumor = jsondecode(fileread('Tumor.mrk.json'));
controlPointsTumor = dataTumor.markups(1).controlPoints;
Tumor = cell2mat(arrayfun(@(p) p.position(:)', controlPointsTumor, 'UniformOutput', false));

dataTrack1 = jsondecode(fileread('Track1.mrk.json'));
controlPointsTrack1 = dataTrack1.markups(1).controlPoints;
pt1_tray = cell2mat(arrayfun(@(p) p.position(:)', controlPointsTrack1, 'UniformOutput', false));

dataTrack2 = jsondecode(fileread('Track2.mrk.json'));
controlPointsTrack2 = dataTrack2.markups(1).controlPoints;
pt2_tray = cell2mat(arrayfun(@(p) p.position(:)', controlPointsTrack2, 'UniformOutput', false));

dataTrack3 = jsondecode(fileread('Track3.mrk.json'));
controlPointsTrack3 = dataTrack3.markups(1).controlPoints;
pt3_tray = cell2mat(arrayfun(@(p) p.position(:)', controlPointsTrack3, 'UniformOutput', false));

dataTrack4 = jsondecode(fileread('Track4.mrk.json'));
controlPointsTrack4 = dataTrack4.markups(1).controlPoints;
pt4_tray = cell2mat(arrayfun(@(p) p.position(:)', controlPointsTrack4, 'UniformOutput', false));

dataTrack5 = jsondecode(fileread('Track5.mrk.json'));
controlPointsTrack5 = dataTrack5.markups(1).controlPoints;
pt5_tray = cell2mat(arrayfun(@(p) p.position(:)', controlPointsTrack5, 'UniformOutput', false));

% Definir Trayectoria
trayectoria = [pt1_tray; pt2_tray; pt3_tray; pt4_tray; pt5_tray];

% Graficar Trayectoria
figure(1);
plot3(trayectoria(:,1), trayectoria(:,2), trayectoria(:,3), '-o', 'LineWidth', 2);
hold on;
scatter3(-centroid_mm(1), -centroid_mm(2), centroid_mm(3), 100, 'r', 'filled');
title('Centroide y trayectoria en coordenadas del paciente');
xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
grid on; xlim([-25 -15]);
legend('Trayectoria', 'Tumor');

%% FASE 2 - MAPEO DE TC Y MOCAP
% 1. Importe como Numeric Matriz una fila del archivo CSV
load("data.mat"); % cargue aquí el archivo .mat de los datos importados

% Hacer un reshape de la data para que cada fila corresponda a las
% coordenadas (x,y,z) de los markers.
marker_data = reshape(seccion11BE1, 3, []);

% Intercambiar coordenadas y con z, cambiar el signo para y.
marker_swap = marker_data([1 3 2], :);
marker_swap(2, :) = -marker_swap(2, :);

% Visualizar markers y determinar cuáles son los de interés
close all;
figure(2);
scatter3(marker_swap(1,:), marker_swap(2,:), marker_swap(3,:), 50, 'filled');
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Markers');
grid on; axis vis3d; hold on;
numPoints = size(marker_swap, 2); % número de columnas (puntos)
for i = 1:numPoints
    text(marker_swap(1,i), marker_swap(2,i), marker_swap(3,i), sprintf('%d', i), ...
        'FontSize', 8, 'FontWeight', 'bold', ...
        'Color', 'k', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom');
end
hold off;

% MARKERS DE INTERES: 9 10 11 12
% Defina el conjunto de puntos deseados
data = marker_swap(:, [9 10 11 12]);

% Visualice los puntos específicos dentro del conjunto deseado
figure(3);
scatter3(data(1,:), data(2,:), data(3,:), 50, 'filled', 'MarkerFaceColor',[0.4940 0.1840 0.5560]);    
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Markers - Proyecto');
grid on; axis vis3d; 

% 2. Importe los puntos de 3D Slicer en formato json
json1 = jsondecode(fileread('Marker1.mrk.json'));
controlPointsjson1 = json1.markups(1).controlPoints;

json2 = jsondecode(fileread('Marker2.mrk.json'));
controlPointsjson2 = json2.markups(1).controlPoints;

json3 = jsondecode(fileread('Marker3.mrk.json'));
controlPointsjson3 = json3.markups(1).controlPoints;

json4 = jsondecode(fileread('Marker4.mrk.json'));
controlPointsjson4 = json4.markups(1).controlPoints;

% Realice la correspondencia de los puntos
% Puntos 3D Slicer (json)
punto1_slicer = cell2mat(arrayfun(@(p) p.position(:)', controlPointsjson1, 'UniformOutput', false));
punto2_slicer = cell2mat(arrayfun(@(p) p.position(:)', controlPointsjson2, 'UniformOutput', false));
punto3_slicer = cell2mat(arrayfun(@(p) p.position(:)', controlPointsjson3, 'UniformOutput', false));
punto4_slicer = cell2mat(arrayfun(@(p) p.position(:)', controlPointsjson4, 'UniformOutput', false));

% Puntos del sistema de captura de movimiento (MOCAP)
punto1_mocap = data(:, 1);
punto2_mocap = data(:, 2);
punto3_mocap = data(:, 3);
punto4_mocap = data(:, 4);

puntos_slicer = [punto1_slicer; punto2_slicer; punto3_slicer; punto4_slicer];
puntos_mocap = [punto1_mocap'; punto2_mocap'; punto3_mocap'; punto4_mocap'];

% Graficar Markers en Slicer y MoCap
figure(4);
hold on;
scatter3(puntos_slicer(:,1), puntos_slicer(:,2), puntos_slicer(:,3), 50, 'filled', 'MarkerFaceColor',[0.4660 0.6740 0.1880]);
scatter3(puntos_mocap(:,1), puntos_mocap(:,2), puntos_mocap(:,3), 50, 'filled', 'MarkerFaceColor',[0.4940 0.1840 0.5560]);
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Markers - Slicer & MoCap');
legend(["Slicer", "MoCap"], 'Location','best');
grid on; axis equal;             
axis vis3d; view(45, 25);           
rotate3d on;            

% 3. Se obtiene la transformación mediante registro
movingPoints = puntos_slicer; 
fixedPoints = puntos_mocap; 
tform = fitgeotform3d(movingPoints, fixedPoints, "affine");

% 4. Se aplica la transformación a los puntos de slicer
mocapT_slicer = tform.A; %matriz de transformación

save('mocapT_slicer.mat', 'mocapT_slicer');

% Agregar un 1 a cada punto (para convertir a coordenadas homogéneas)
punto1_h = [puntos_slicer(1,:) 1]';
punto2_h = [puntos_slicer(2,:) 1]';
punto3_h = [puntos_slicer(3,:) 1]';
punto4_h = [puntos_slicer(4,:) 1]';

% aplique la matriz a cada uno de los puntos de 3D slicer
punto1_trans_h = mocapT_slicer * punto1_h;
punto2_trans_h = mocapT_slicer * punto2_h;
punto3_trans_h = mocapT_slicer * punto3_h;
punto4_trans_h = mocapT_slicer * punto4_h;

% Quitar el 1 (regresar a coordenadas cartesianas)
punto1_trans = punto1_trans_h(1:3)';
punto2_trans = punto2_trans_h(1:3)';
punto3_trans = punto3_trans_h(1:3)';
punto4_trans = punto4_trans_h(1:3)';

data_trans = [punto1_trans; punto2_trans; punto3_trans; punto4_trans]; % agrupe en un vector la data final

% Verifique la transformación graficando los dos conjuntos de puntos y como se sobreponen
figure(5);
hold on;
scatter3(puntos_mocap(:,1), puntos_mocap(:,2), puntos_mocap(:,3), 100, 'o', 'MarkerEdgeColor', [0.4940 0.1840 0.5560], 'LineWidth', 1.5);
scatter3(data_trans(:,1), data_trans(:,2), data_trans(:,3), 40, 'filled', 'MarkerFaceColor', [0.3010 0.7450 0.9330]);
legend({'MOCAP (Referencia)', 'Slicer Transformado'}, 'Location','best');
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Verificación del Registro 3D');
grid on; axis equal; 
axis vis3d; view(45,25);
hold off;

%% FASE 3 - DEFINICIÓN DE ROBOT Y TRAYECTORIA PARA EL PROCEDIMIENTO
% Se carga el modelo cinemático del robot. El robot UR5 se encuentra dentro
% de los modelos contenidos por defecto en la Robotics Toolbox
mdl_ur5

% El modelo aparece como ur5 en el workspace y también crea dos 
% configuraciones básicas qr y qv
q0 = zeros(1, 7);

%% 2. Mapear centride del Tumor y pasos de la trayectoria a posición en
% MoCap

% Agregar un 1 a cada punto (para convertir a coordenadas homogéneas)
puntoTumor_h = [Tumor 1]';

puntoTray1_h = [pt1_tray 1]';
puntoTray2_h = [pt2_tray 1]';
puntoTray3_h = [pt3_tray 1]';
puntoTray4_h = [pt4_tray 1]';
puntoTray5_h = [pt5_tray 1]';

% aplique la matriz a cada uno de los puntos de 3D slicer
puntoTumor_trans_h = mocapT_slicer * puntoTumor_h;
puntoTray1_trans_h = mocapT_slicer * puntoTray1_h;
puntoTray2_trans_h = mocapT_slicer * puntoTray2_h;
puntoTray3_trans_h = mocapT_slicer * puntoTray3_h;
puntoTray4_trans_h = mocapT_slicer * puntoTray4_h;
puntoTray5_trans_h = mocapT_slicer * puntoTray5_h;

% Quitar el 1 (regresar a coordenadas cartesianas)
puntoTumor_trans = puntoTumor_trans_h(1:3)';
puntoTray1_trans = puntoTray1_trans_h(1:3)';
puntoTray2_trans = puntoTray2_trans_h(1:3)';
puntoTray3_trans = puntoTray3_trans_h(1:3)';
puntoTray4_trans = puntoTray4_trans_h(1:3)';
puntoTray5_trans = puntoTray5_trans_h(1:3)';

Tray_trans = [puntoTray1_trans; puntoTray2_trans; puntoTray3_trans; puntoTray4_trans; puntoTray5_trans]; % agrupe en un vector la data final

figure(6);
hold on;
scatter3(puntos_mocap(:,1), puntos_mocap(:,2), puntos_mocap(:,3), 100, 'o', 'MarkerEdgeColor', [0.4940 0.1840 0.5560], 'LineWidth', 1.5);
scatter3(data_trans(:,1), data_trans(:,2), data_trans(:,3), 40, 'filled', 'MarkerFaceColor', [0.3010 0.7450 0.9330]);
% scatter3(Tray_trans(:,1), Tray_trans(:,2), Tray_trans(:,3), 40, 'filled', 'MarkerFaceColor', [0.4660 0.6740 0.1880]);
plot3(Tray_trans(:,1), Tray_trans(:,2), Tray_trans(:,3), '-o', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 1.5);
scatter3(puntoTumor_trans(1,1), puntoTumor_trans(1,2), puntoTumor_trans(1,3), 40, 'filled', 'MarkerFaceColor', [0.6350 0.0780 0.1840]);
legend({'MOCAP (Referencia)', 'Slicer Transformado', 'Trayectoria', 'Centroide Tumor'}, 'Location','best');
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Mapeo de Trayectoria y Centoide Tumor a MoCap');
grid on; axis equal; 
axis vis3d; view(45,25);
hold off;
%%


%% Inciso 3
obs1 = [0.7; 0.25; 0.9];
obs2 = [0.3; -0.25; 0.9];
figure('units', 'normalized', 'outerposition', [0 0 1 1]);
plot_sphere(pos1, 0.05, 'b');
plot_sphere(pos2, 0.05, 'b');
plot_sphere(obs1, 0.05, 'r');
plot_sphere(obs2, 0.05, 'r');

via1_3 = [0.7; 0.15; 0.5]; % Definir via  

Tvia1_3 = transl(via1_3); % Convertir en transformaciones homogeneas 

Ts1_3 = ctraj(T1, Tvia1_3, 100);   % tramo 1
Ts2_3 = ctraj(Tvia1_3, T2 * troty(90), 100);   % tramo 2
TsT_3 = cat(3, Ts1_3, Ts2_3);  % concatenar trayectorias

Q2 = panda.ikine(TsT_3, 'mask', [1 1 1 0 1 0]);

panda.plot(Q2, 'zoom', 1.5, 'delay', 0.05);

%% Inciso 4 
head = [0.8; 0.8; 0.2];
pos3 = [0.5; 0.5; 0.5];
figure('units', 'normalized', 'outerposition', [0 0 1 1]);
plot_sphere(head, 0.5, 'r', 0.5);
plot_sphere(pos3, 0.05, 'g');
hold on;
Td1 = transl(0.5, 0.5, 0.5) * trotx(-135) * trotz(-45);
Td2 = Td1 * transl(0 , 0, 0.2);
trplot(Td1, 'frame', 'A', 'color', 'g');
trplot(Td2, 'frame', 'B', 'color', 'c');
hold off;

T0_4 = panda.fkine(q0).T; % Extraer matriz homogenea para iniciar traslación                 

Ts1_4 = ctraj(T0_4, Td1, 100);
Ts2_4 = ctraj(Td1, Td2, 100);

Q3_1 = panda.ikine(Ts1_4, 'mask', [1 1 1 0 0 0], 'q0', q0);             
Q3_2 = panda.ikine(Ts2_4, 'mask', [1 1 1 0 0 0], 'q0', Q3_1(end,:));   

Q3 = [Q3_1; Q3_2]; % Concatenar Trayectorias

panda.plot(Q3, 'zoom', 1.5, 'delay', 0.05);
%% Inciso 5
pos1_5 = [0; 0.5; 0.3]; % Definir puntos de Trayectorias
pos2_5 = [0; 0.5; 1.1];
pos3_5 = [0.2; 0.5; 1.08];
pos4_5 = [0.3; 0.5; 0.99]; 
pos5_5 = [0.38; 0.5; 0.9];
pos6_5 = [0.4; 0.5; 0.8];
pos7_5 = [0.41; 0.5; 0.7];
pos8_5 = [0.4; 0.5; 0.6];
pos9_5 = [0.38; 0.5; 0.5];
pos10_5 = [0.3; 0.5; 0.4];
pos11_5 = [0.3; 0.5; 0.4];

T1_5 = transl(pos1_5); % Convetir a matricez homogeneas 
T2_5 = transl(pos2_5);
T3_5 = transl(pos3_5);
T4_5 = transl(pos4_5);
T5_5 = transl(pos5_5);
T6_5 = transl(pos6_5);
T7_5 = transl(pos7_5);
T8_5 = transl(pos8_5);
T9_5 = transl(pos9_5);
T10_5 = transl(pos10_5);
T11_5 = transl(pos11_5);

Ts1_5 = ctraj(T1_5, T2_5, 100); % Tramos del Recorrido
Ts2_5 = ctraj(T2_5, T3_5, 50);
Ts3_5 = ctraj(T3_5, T4_5, 50);
Ts4_5 = ctraj(T4_5, T5_5, 50);
Ts5_5 = ctraj(T5_5, T6_5, 50);
Ts6_5 = ctraj(T6_5, T7_5, 50);
Ts7_5 = ctraj(T7_5, T8_5, 50);
Ts8_5 = ctraj(T8_5, T9_5, 50);
Ts9_5 = ctraj(T9_5, T10_5, 50);
Ts10_5 = ctraj(T10_5, T11_5, 50);
Ts11_5 = ctraj(T11_5, T1_5, 50);

Q4_1 = panda.ikine(Ts1_5, 'mask', [1 1 1 0 0 0], 'q0', q0); % Trayectorias
Q4_2 = panda.ikine(Ts2_5, 'mask', [1 1 1 0 0 0], 'q0', Q4_1(end,:));
Q4_3 = panda.ikine(Ts3_5, 'mask', [1 1 1 0 0 0], 'q0', Q4_2(end,:));
Q4_4 = panda.ikine(Ts4_5, 'mask', [1 1 1 0 0 0], 'q0', Q4_3(end,:));
Q4_5 = panda.ikine(Ts5_5, 'mask', [1 1 1 0 0 0], 'q0', Q4_4(end,:));
Q4_6 = panda.ikine(Ts6_5, 'mask', [1 1 1 0 0 0], 'q0', Q4_5(end,:));
Q4_7 = panda.ikine(Ts7_5, 'mask', [1 1 1 0 0 0], 'q0', Q4_6(end,:));
Q4_8 = panda.ikine(Ts8_5, 'mask', [1 1 1 0 0 0], 'q0', Q4_7(end,:));
Q4_9 = panda.ikine(Ts9_5, 'mask', [1 1 1 0 0 0], 'q0', Q4_8(end,:));
Q4_10 = panda.ikine(Ts10_5, 'mask', [1 1 1 0 0 0], 'q0', Q4_9(end,:));
Q4_11 = panda.ikine(Ts11_5, 'mask', [1 1 1 0 0 0], 'q0', Q4_10(end,:));

Q4 = [Q4_1; Q4_2; Q4_3; Q4_4; Q4_5; Q4_6; Q4_7; Q4_8; Q4_9; Q4_10; Q4_11]; % Concatenar Trayectorias

figure; % Figura
panda.plot(Q4, 'zoom', 1.5, 'delay', 0.005);
view(0, 10); % Ajustar vista

traj = []; % Vector para extraer puntos del recorrido
for i = 1:size(Ts1_5,3)
    traj = [traj; transl(Ts1_5(:,:,i))'];
end
for i = 1:size(Ts2_5,3)
    traj = [traj; transl(Ts2_5(:,:,i))'];
end
for i = 1:size(Ts3_5,3)
    traj = [traj; transl(Ts3_5(:,:,i))'];
end
for i = 1:size(Ts4_5,3)
    traj = [traj; transl(Ts4_5(:,:,i))'];
end
for i = 1:size(Ts5_5,3)
    traj = [traj; transl(Ts5_5(:,:,i))'];
end
for i = 1:size(Ts6_5,3)
    traj = [traj; transl(Ts6_5(:,:,i))'];
end
for i = 1:size(Ts7_5,3)
    traj = [traj; transl(Ts7_5(:,:,i))'];
end
for i = 1:size(Ts8_5,3)
    traj = [traj; transl(Ts8_5(:,:,i))'];
end
for i = 1:size(Ts9_5,3)
    traj = [traj; transl(Ts9_5(:,:,i))'];
end
for i = 1:size(Ts10_5,3)
    traj = [traj; transl(Ts10_5(:,:,i))'];
end
for i = 1:size(Ts11_5,3)
    traj = [traj; transl(Ts11_5(:,:,i))'];
end

hold on; % Graficar puntos
plot3(traj(:,1), traj(:,2), traj(:,3), 'b-', 'LineWidth', 2);
grid on; axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');












%% FASE 4 - EJECUTAR TRAYECTORIA 