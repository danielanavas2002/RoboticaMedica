% ------------------------------------------------------------------- 
% PROYECTO FINAL
% ROBÓTICA MÉDICA - 2025
% 
% Aury Campaneros & Daniela Navas 
% ------------------------------------------------------------------
close all; clear; clc;
%% Fase 1 - Segmentación, markups y trayectoria
% 1. Carga de información de archivos para uso en la App "Medical Image Labeler"
filepath = 'C:\Users\Daniela\Desktop\Robótica Médica\RoboticaMedica\Matlab\4.000000-Ax T1 MP SPGR-PRE GAD-85514';
dataFolder = fullfile(filepath);

% Extraer la información general de la imagen
medVol = medicalVolume(dataFolder);

% 2. Definición de Centro de Masa (centroide del tumor)
load('C:\Users\Daniela\Desktop\Robótica Médica\RoboticaMedica\Matlab\groundTruthMed.mat')

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
dicom_prop = dicominfo('C:\Users\Daniela\Desktop\Robótica Médica\RoboticaMedica\Matlab\4.000000-Ax T1 MP SPGR-PRE GAD-85514\1-001.dcm');

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
dataTumor = jsondecode(fileread('C:\Users\Daniela\Desktop\Robótica Médica\Proyecto\Slicer\Tumor.mrk.json'));
controlPointsTumor = dataTumor.markups(1).controlPoints;
Tumor = cell2mat(arrayfun(@(p) p.position(:)', controlPointsTumor, 'UniformOutput', false));

dataTrack1 = jsondecode(fileread('C:\Users\Daniela\Desktop\Robótica Médica\Proyecto\Slicer\Track1.mrk.json'));
controlPointsTrack1 = dataTrack1.markups(1).controlPoints;
pt1_tray = cell2mat(arrayfun(@(p) p.position(:)', controlPointsTrack1, 'UniformOutput', false));

dataTrack2 = jsondecode(fileread('C:\Users\Daniela\Desktop\Robótica Médica\Proyecto\Slicer\Track2.mrk.json'));
controlPointsTrack2 = dataTrack2.markups(1).controlPoints;
pt2_tray = cell2mat(arrayfun(@(p) p.position(:)', controlPointsTrack2, 'UniformOutput', false));

dataTrack3 = jsondecode(fileread('C:\Users\Daniela\Desktop\Robótica Médica\Proyecto\Slicer\Track3.mrk.json'));
controlPointsTrack3 = dataTrack3.markups(1).controlPoints;
pt3_tray = cell2mat(arrayfun(@(p) p.position(:)', controlPointsTrack3, 'UniformOutput', false));

dataTrack4 = jsondecode(fileread('C:\Users\Daniela\Desktop\Robótica Médica\Proyecto\Slicer\Track4.mrk.json'));
controlPointsTrack4 = dataTrack4.markups(1).controlPoints;
pt4_tray = cell2mat(arrayfun(@(p) p.position(:)', controlPointsTrack4, 'UniformOutput', false));

dataTrack5 = jsondecode(fileread('C:\Users\Daniela\Desktop\Robótica Médica\Proyecto\Slicer\Track5.mrk.json'));
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

%% Fase 2 - Transformación entre sistemas de coordenadas mediante registro
% 1. Importe como Numeric Matriz una fila del archivo CSV
load("C:\Users\Daniela\Desktop\Robótica Médica\Proyecto\Matlab\data.mat"); % cargue aquí el archivo .mat de los datos importados

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
json1 = jsondecode(fileread('C:\Users\Daniela\Desktop\Robótica Médica\Laboratorio\Lab10\Marker1.mrk.json'));
controlPointsjson1 = json1.markups(1).controlPoints;

json2 = jsondecode(fileread('C:\Users\Daniela\Desktop\Robótica Médica\Laboratorio\Lab10\Marker2.mrk.json'));
controlPointsjson2 = json2.markups(1).controlPoints;

json3 = jsondecode(fileread('C:\Users\Daniela\Desktop\Robótica Médica\Laboratorio\Lab10\Marker3.mrk.json'));
controlPointsjson3 = json3.markups(1).controlPoints;

json4 = jsondecode(fileread('C:\Users\Daniela\Desktop\Robótica Médica\Laboratorio\Lab10\Marker4.mrk.json'));
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

% Mapear centroide del tumor y trayectoria a sistema MoCap
puntoTumor = Tumor(1, :);   % [x y z] en Slicer
puntoTray1 = pt1_tray(1, :);
puntoTray2 = pt2_tray(1, :);
puntoTray3 = pt3_tray(1, :);
puntoTray4 = pt4_tray(1, :);
puntoTray5 = pt5_tray(1, :);

% Pasar a homogéneas
puntoTumor_h  = [puntoTumor  1]';
puntoTray1_h  = [puntoTray1  1]';
puntoTray2_h  = [puntoTray2  1]';
puntoTray3_h  = [puntoTray3  1]';
puntoTray4_h  = [puntoTray4  1]';
puntoTray5_h  = [puntoTray5  1]';

% Aplicar la transformación Slicer -> MoCap
puntoTumor_trans_h = mocapT_slicer * puntoTumor_h;
puntoTray1_trans_h = mocapT_slicer * puntoTray1_h;
puntoTray2_trans_h = mocapT_slicer * puntoTray2_h;
puntoTray3_trans_h = mocapT_slicer * puntoTray3_h;
puntoTray4_trans_h = mocapT_slicer * puntoTray4_h;
puntoTray5_trans_h = mocapT_slicer * puntoTray5_h;

% Quitar el 1 (volver a cartesiano)
puntoTumor_trans = puntoTumor_trans_h(1:3)';
puntoTray1_trans = puntoTray1_trans_h(1:3)';
puntoTray2_trans = puntoTray2_trans_h(1:3)';
puntoTray3_trans = puntoTray3_trans_h(1:3)';
puntoTray4_trans = puntoTray4_trans_h(1:3)';
puntoTray5_trans = puntoTray5_trans_h(1:3)';

% Trayectoria en sistema MoCap
Tray_trans = [puntoTray1_trans;
              puntoTray2_trans;
              puntoTray3_trans;
              puntoTray4_trans;
              puntoTray5_trans];

% Si ves que las coordenadas están en mm, descomenta esto:
% puntoTumor_trans = puntoTumor_trans / 1000;
% Tray_trans       = Tray_trans / 1000;

% Gráfica: Verificación Tumor + Trayectoria en Sistema MoCap
figure(6); clf; hold on;
scatter3(puntos_mocap(:,1), puntos_mocap(:,2), puntos_mocap(:,3), 100, 'o', 'MarkerEdgeColor', [0.4940 0.1840 0.5560], 'LineWidth', 1.5);
scatter3(data_trans(:,1), data_trans(:,2), data_trans(:,3), 40, 'filled', 'MarkerFaceColor', [0.3010 0.7450 0.9330]);
plot3(Tray_trans(:,1), Tray_trans(:,2), Tray_trans(:,3), '-o', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 1.5, 'MarkerFaceColor', [0.4660 0.6740 0.1880]);
scatter3(puntoTumor_trans(1), puntoTumor_trans(2), puntoTumor_trans(3), 120, 'filled', 'MarkerFaceColor', [0.6350 0.0780 0.1840]);
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Tumor y Trayectoria en Coordenadas del Sistema MoCap');
legend({'MoCap (Referencia)', 'Markers Slicer', 'Trayectoria', 'Centroide Tumor'},'Location', 'best');
grid on; axis equal; 
axis vis3d; view(45, 25);
hold off;
%% Fase 3 - Visualización de la cabeza y patología mediante la Robotics Toolbox 
head_center = mean(data, 2) - [40; 0; 0];   

% Convertir mm → m
head_center_m = head_center / 1000;
puntoTumor_trans_m = puntoTumor_trans / 1000;
Tray_trans_m = Tray_trans / 1000;
data_trans_m = data_trans / 1000;       
puntos_mocap_m = puntos_mocap / 1000;   
escala = 3;

head_center_scaled = head_center_m * escala;
puntoTumor_scaled  = puntoTumor_trans_m * escala;
Tray_scaled        = Tray_trans_m * escala;

radio_cabeza = 0.13 * escala;
radio_tumor  = 0.01 * escala;
tam_mesa     = 0.25 * escala;
altura_mesa  = head_center_scaled(3) - radio_cabeza - 0.05;

[Xm,Ym] = meshgrid(-tam_mesa:0.02*escala:tam_mesa);
Zm = ones(size(Xm))*altura_mesa;

figure(7); clf; hold on; grid on; axis equal; axis vis3d;

surf(Xm, Ym, Zm,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none','FaceAlpha',0.35);

plot_sphere(head_center_scaled, radio_cabeza, 'b', 0.25);
plot_sphere(puntoTumor_scaled', radio_tumor, 'r', 1.0);
plot3(Tray_scaled(:,1),Tray_scaled(:,2),Tray_scaled(:,3),'-og','LineWidth',2);

xlabel('X'); ylabel('Y'); zlabel('Z');
title('Escena completa escalada (escala = 3)');
xlim([head_center_scaled(1)-1*escala, head_center_scaled(1)+1*escala]);
ylim([head_center_scaled(2)-1*escala, head_center_scaled(2)+1*escala]);
zlim([altura_mesa-0.2*escala, head_center_scaled(3)+0.6*escala]);

hold off;

%% Fase 4 — Instalación del manipulador serial dentro de la simulación
mdl_panda
q0 = zeros(1,7);

offset_base = [
    head_center_scaled(1) - 0.35*escala
    head_center_scaled(2) - 0.20*escala
    altura_mesa
];
panda.base = transl(offset_base);
panda.tool = transl(0,0,0.10*escala) * trotx(pi);

figure(8); clf; hold on; grid on; axis equal; axis vis3d;

surf(Xm, Ym, Zm,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none','FaceAlpha',0.35);

plot_sphere(head_center_scaled, radio_cabeza, 'b', 0.25);
plot_sphere(puntoTumor_scaled', radio_tumor, 'r', 1.0);
plot3(Tray_scaled(:,1),Tray_scaled(:,2),Tray_scaled(:,3),'-og','LineWidth',2);

panda.plot(q0,'view',[45 25],'nojoints','noname');

xlim([head_center_scaled(1)-1*escala, head_center_scaled(1)+1*escala]);
ylim([head_center_scaled(2)-1*escala, head_center_scaled(2)+1*escala]);
zlim([altura_mesa-0.2*escala, head_center_scaled(3)+0.6*escala]);

title('Robot Panda con mundo escalado');
xlabel('X'); ylabel('Y'); zlabel('Z');
hold off;
%% Fase 5 — Movimiento del robot solo siguiendo la trayectoria verde

%% Fase 3

%% Fase 4

%% Fase 5
