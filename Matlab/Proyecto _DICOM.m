% =========================================================================
% BE3027 - LABORATORIO 9: Procesamiento de imágenes médicas utilizando el
%          "Medical Imaging Toolbox" de MATLAB
% -------------------------------------------------------------------------
% Ver las instrucciones en la guía adjunta
% =========================================================================
%% 6. Repetición del proceso con archivos DICOM
% 2. Carga de información de archivos para uso en la App "Medical Image Labeler"
% UTILIZANDO ARCHIVOS NIFTI

% Guarde en la variable "filepath" la dirección de la carpeta en donde se
% encuentra la carpeta con el conjunto de imágenes médicas
% Nota: la dirección no debe incluir el la carpeta "MedicalVolumeDIfTIData"
% sino que llegar hasta la ubicación anterior
filepath = 'C:\Users\Daniela\Desktop\Robótica Médica\Proyecto\Matlab\4.000000-Ax T1 MP SPGR-PRE GAD-85514';
dataFolder = fullfile(filepath);

% La función "medicalVolume" detecta de forma automática el formato de los
% archivos y extrae la información general de la imagen

medVol = medicalVolume(dataFolder);

%% 3. Definición de Centro de Masa (centroide del tumor)
load('C:\Users\Daniela\Desktop\Robótica Médica\Proyecto\Matlab\groundTruthMed.mat')

% 3.(a) Obtención de Voxeles que conforman el volumen definido en la App

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
labelName = 'Tumor';  %NOTA: este nombre debe coincidir con el nombre utilizado en el app

% Seleccionamos el ID asignado a la "Label" que deseamos analizar (si solo
% hay una label definida sería la 1)

% la instrucción "strcmp" sirve para realizar una comparación de texto
% entonces tomaremos en cuenta solo las etiquetas que correspondan al
% nombre de la región definida anteriormente 
labelID = labelDefs.PixelLabelID(strcmp(labelDefs.Name, labelName));

% Utilizando la función "find" haremos un arreglo con los datos de
% "labelVolume" que cumplan con el índice la región en cuestion, lo que nos
% proporcionará la localización (coordenadas intrínsicas) de los Voxels 
% correspndientes

RegionVoxels = find(labelVolume == labelID);

%% 3.(b) Generación de Matriz de Voxels 
% Luego encontraremos la localización en el formato (fila, columna, slice)
% de cada uno de los Voxels encontrados en el inciso anterior. Esto se
% registra desde la perspectiva transversal.
[rows, cols, slices] = ind2sub(size(labelVolume), RegionVoxels);

% Finalmente guardaremos la información de todos los Voxels en una sola
% matriz
voxelLocation = ([rows, cols, slices]);

%% 3.(c) Cálculo de centroide (centro de masa) de la región

% Para calcular el centroide sacaremos el promedio de cada una de las 
% coordenadas de los Voxels, para esto utilizaremos la siguiente notación

% rows = i, cols = j, slices = k
cent_i = mean(rows);
cent_j = mean(cols);
cent_k = mean(slices);

centroid_voxels = [cent_i, cent_j, cent_k];

%% 3.(d) Conversión a Coordenadas del Paciente (en mm)
% Para esta parte necesitamos leer las propiedades del Nifti o Dicom
% original, en donde se encuentra la trasformación utilizada al generar las
% imágenes. para esto debemos realizar lo siguiente: 

% 1. Leer las propiedades del set de imágenes original utilizando la
% función niftiinfo (para archivos DICOM la instrucción es dicominfo)
% nifti_prop = niftiinfo(labelFilePath);

dicom_prop = dicominfo('C:\Users\Daniela\Desktop\Robótica Médica\Proyecto\Matlab\4.000000-Ax T1 MP SPGR-PRE GAD-85514\1-001.dcm');

% Para este paso es necesario utilizar la transformación T. Hay dos
% elementos de alta importancia. La primera es la función affinetform3d, la
% cual se encarga de "ajustar" la matriz de transofrmación en un formato 
% adecuado para representa cómo se deben mover, rotar, escalar o 
% distorsionar puntos o imágenes en un espacio 3D.

% Luego está el T, esta es la matriz de transformación utilizada al generar
% las imágenes, sin embargo debe utilizarse su transpuesta debido a que el
% toobox maneja una orientación distinta para esta matriz. 
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
% delta_z = dicom_prop.SliceThickness;

T_1 = [ B  C BxC pos;
        0  0  0   1 ]; 

T_2 = [delta_x    0       0       0   ;
          0    delta_y    0       0   ;
          0       0    delta_z    0   ;
          0       0       0       1   ];

T = T_1 * T_2;

matriz_trans_3d = affinetform3d(T);

% 2. Creamos el objeto de referencia espacial tridimencional, el cual se
% calcula a partir de las dimenciones originales (en voxels) de las 
% imágenes y la matriz de transofrmaicón

dicomFiles = dir(fullfile(filepath, '*.dcm')); % lista de archivos DICOM
numSlices = length(dicomFiles);                 % número real de slices
imageSize = double([dicom_prop.Rows, dicom_prop.Columns, numSlices]);

ref3d = medicalref3d(imageSize, matriz_trans_3d);

% 3. Aplicamos la función Intrinsic (Voxels) a World (Paciente)
% Para este paso debemos de tomar en cuenta que existen ciertas diferencias
% entre el sistema de coordinadas de Matlab con 3D Slicer y 
[x_paciente, y_paciente, z_paciente] = intrinsicToWorld(ref3d, cent_j, cent_i, cent_k); %Ingrese aquí los parámetros de la función

% centroid_mm = [x_paciente, -y_paciente, z_paciente-161]; %cree un solo elemento con las coordenadas del centroide
centroid_mm = [x_paciente, y_paciente, z_paciente];

%% 5.Gráfica 3D
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

trayectoria = [pt1_tray; pt2_tray; pt3_tray; pt4_tray; pt5_tray];

figure;
plot3(trayectoria(:,1), trayectoria(:,2), trayectoria(:,3), '-o', 'LineWidth', 2);
hold on;
scatter3(-centroid_mm(1), -centroid_mm(2), centroid_mm(3), 100, 'r', 'filled');
title('Centroide y trayectoria en coordenadas del paciente');
xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
grid on; xlim([-25 -15]);
legend('Trayectoria', 'Tumor');

%% 3.(d) Conversión a Coordenadas del Paciente (en mm)

% Para esta parte necesitamos leer las propiedades del Nifti o Dicom
% original, en donde se encuentra la trasformación utilizada al generar las
% imágenes. para esto debemos realizar lo siguiente: 

% 1. Leer las propiedades del set de imágenes original utilizando la
% función niftiinfo (para archivos DICOM la instrucción es dicominfo)
dicom_prop = niftiinfo(labelFilePath);

% Para este paso es necesario utilizar la transformación T. Hay dos
% elementos de alta importancia. La primera es la función affinetform3d, la
% cual se encarga de "ajustar" la matriz de transofrmación en un formato 
% adecuado para representa cómo se deben mover, rotar, escalar o 
% distorsionar puntos o imágenes en un espacio 3D.

% Luego está el T, esta es la matriz de transformación utilizada al generar
% las imágenes, sin embargo debe utilizarse su transpuesta debido a que el
% toobox maneja una orientación distinta para esta matriz. 

%matriz_trans_3d = affinetform3d(dicom_prop.Transform.T');
matriz_trans_3d = affinetform3d(dicom_prop.Transform.T');

% 2. Creamos el objeto de referencia espacial tridimencional, el cual se
% calcula a partir de las dimenciones originales (en voxels) de las 
% imágenes y la matriz de transofrmaicón
%ref3d = medicalref3d(dicom_prop.ImageSize(1:3), matriz_trans_3d);
ref3d = medicalref3d(dicom_prop.ImageSize(1:3), matriz_trans_3d);

% 3. Aplicamos la función Intrinsic (Voxels) a World (Paciente)
% Para este paso debemos de tomar en cuenta que existen ciertas diferencias
% entre el sistema de coordinadas de Matlab con 3D Slicer y 

[x_paciente, y_paciente, z_paciente] = intrinsicToWorld(ref3d, cent_j, cent_i, cent_k); %Ingrese aquí los parámetros de la función

centroid_mm = [-x_paciente, y_paciente, z_paciente]; %cree un solo elemento con las coordenadas del centroide