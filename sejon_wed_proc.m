%% Load data
clear; close all; clc;

load('D:\2022_Spring\research\학위논문연구\dataset\SejongDataset\sc1_data.mat');
load('D:\2022_Spring\research\학위논문연구\dataset\SejongDataset\sc2_data.mat');
load('D:\2022_Spring\research\학위논문연구\dataset\SejongDataset\sc3_data.mat');
load('D:\2022_Spring\research\학위논문연구\dataset\SejongDataset\sc4_data.mat');
load('D:\2022_Spring\research\학위논문연구\dataset\SejongDataset\sc5_data.mat');
% sc3_data1 = load('D:\2022_Spring\research\학위논문연구\dataset\SejongDataset\sc3_data1.mat');
% sc3_data2 = load('D:\2022_Spring\research\학위논문연구\dataset\SejongDataset\sc3_data2.mat');
% sc3_data3 = load('D:\2022_Spring\research\학위논문연구\dataset\SejongDataset\sc3_data3.mat');
% sc3_data4 = load('D:\2022_Spring\research\학위논문연구\dataset\SejongDataset\sc3_data4.mat');
% sc4_data1 = load('D:\2022_Spring\research\학위논문연구\dataset\SejongDataset\sc4_data1.mat');
% sc4_data2 = load('D:\2022_Spring\research\학위논문연구\dataset\SejongDataset\sc4_data2.mat');
% sc5_data1 = load('D:\2022_Spring\research\학위논문연구\dataset\SejongDataset\sc5_data1.mat');
% sc5_data2 = load('D:\2022_Spring\research\학위논문연구\dataset\SejongDataset\sc5_data2.mat');
% sc5_data3 = load('D:\2022_Spring\research\학위논문연구\dataset\SejongDataset\sc5_data3.mat');

load('D:\2022_Spring\research\학위논문연구\dataset\SejongDataset\ppk_data.mat')
% 
% 
sc1_data.gnssppk = ppk_data.sc1;
sc2_data.gnssppk = ppk_data.sc2;

% % Merge sliced data sc3
% 
% sc3_data1.lane.t = sc3_data1.lane.t';
% sc3_data2.lane.t = sc3_data2.lane.t';
% sc3_data3.lane.t = sc3_data3.lane.t';
% sc3_data4.lane.t = sc3_data4.lane.t';
% sc3_data1.imu.t = sc3_data1.imu.t';
% sc3_data2.imu.t = sc3_data2.imu.t';
% sc3_data3.imu.t = sc3_data3.imu.t';
% sc3_data4.imu.t = sc3_data4.imu.t';
% 
% sc3_data1.can = FullTranspose(sc3_data1.can);
% sc3_data2.can = FullTranspose(sc3_data2.can);
% sc3_data3.can = FullTranspose(sc3_data3.can);
% sc3_data4.can = FullTranspose(sc3_data4.can);
% sc3_data1.gnss = FullTranspose(sc3_data1.gnss);
% sc3_data2.gnss = FullTranspose(sc3_data2.gnss);
% sc3_data3.gnss = FullTranspose(sc3_data3.gnss);
% sc3_data4.gnss = FullTranspose(sc3_data4.gnss);
% 
% sc3_data = struct();
% sc3_data.lane = struct();
% sc3_data.can = struct();
% sc3_data.gnss = struct();
% sc3_data.imu = struct();
% 
% sc3_data.lane = CatStructFields(1,sc3_data1.lane,sc3_data2.lane,sc3_data3.lane,sc3_data4.lane);
% sc3_data.can = CatStructFields(1,sc3_data1.can,sc3_data2.can,sc3_data3.can,sc3_data4.can);
% sc3_data.gnss = CatStructFields(1,sc3_data1.gnss,sc3_data2.gnss,sc3_data3.gnss,sc3_data4.gnss);
% sc3_data.imu = CatStructFields(1,sc3_data1.imu,sc3_data2.imu,sc3_data3.imu,sc3_data4.imu);
sc3_data.gnssppk = ppk_data.sc3;
% 
% 
% % Merge sliced data sc4
% 
% sc4_data1.lane.t = sc4_data1.lane.t';
% sc4_data2.lane.t = sc4_data2.lane.t';
% sc4_data1.imu.t = sc4_data1.imu.t';
% sc4_data2.imu.t = sc4_data2.imu.t';
% sc4_data1.can = FullTranspose(sc4_data1.can);
% sc4_data2.can = FullTranspose(sc4_data2.can);
% sc4_data1.gnss = FullTranspose(sc4_data1.gnss);
% sc4_data2.gnss = FullTranspose(sc4_data2.gnss);
% 
% sc4_data = struct();
% sc4_data.lane = struct();
% sc4_data.can = struct();
% sc4_data.gnss = struct();
% sc4_data.imu = struct();
% 
% sc4_data.lane = CatStructFields(1,sc4_data1.lane,sc4_data2.lane);
% sc4_data.can = CatStructFields(1,sc4_data1.can,sc4_data2.can);
% sc4_data.gnss = CatStructFields(1,sc4_data1.gnss,sc4_data2.gnss);
% sc4_data.imu = CatStructFields(1,sc4_data1.imu,sc4_data2.imu);
sc4_data.gnssppk = ppk_data.sc4;
% % 
% % Merge sliced data sc5
% 
% sc5_data1.lane.t = sc5_data1.lane.t';
% sc5_data2.lane.t = sc5_data2.lane.t';
% sc5_data3.lane.t = sc5_data3.lane.t';
% sc5_data1.imu.t = sc5_data1.imu.t';
% sc5_data2.imu.t = sc5_data2.imu.t';
% sc5_data3.imu.t = sc5_data3.imu.t';
% sc5_data1.can = FullTranspose(sc5_data1.can);
% sc5_data2.can = FullTranspose(sc5_data2.can);
% sc5_data3.can = FullTranspose(sc5_data3.can);
% sc5_data1.gnss = FullTranspose(sc5_data1.gnss);
% sc5_data2.gnss = FullTranspose(sc5_data2.gnss);
% sc5_data3.gnss = FullTranspose(sc5_data3.gnss);
% 
% sc5_data = struct();
% sc5_data.lane = struct();
% sc5_data.can = struct();
% sc5_data.gnss = struct();
% sc5_data.imu = struct();
% 
% sc5_data.lane = CatStructFields(1,sc5_data1.lane,sc5_data2.lane,sc5_data3.lane);
% sc5_data.can = CatStructFields(1,sc5_data1.can,sc5_data2.can,sc5_data3.can);
% sc5_data.gnss = CatStructFields(1,sc5_data1.gnss,sc5_data2.gnss,sc5_data3.gnss);
% sc5_data.imu = CatStructFields(1,sc5_data1.imu,sc5_data2.imu,sc5_data3.imu);
sc5_data.gnssppk = ppk_data.sc5;

%% Associate timestamp

% reference: imu, lane
% need to change: can, gnssrtk

% n = length(sc2_data.lane.t);
% sc1_data.can.idxs = zeros(1,n); % Effective CAN idxs
% for i=1:n
%     [~,idx] = min((sc1_data.lane.t(i) - sc1_data.can.t).^2);
%     sc1_data.can.idxs(i) = idx;
% end
% sc1_data.can.aligned_ws_rl = sc1_data.can.ws_rl(sc1_data.can.idxs);
% sc1_data.can.aligned_ws_rr = sc1_data.can.ws_rr(sc1_data.can.idxs);
% sc1_data.can.aligned_ws_lB = sc1_data.can.leftBlinker(sc1_data.can.idxs);
% sc1_data.can.aligned_ws_rB = sc1_data.can.rightBlinker(sc1_data.can.idxs);


n = length(sc2_data.gnssppk.t);
sc2_data.gnssppk.idxs = zeros(1,n);
for i=1:n
    [~,idx] = min((sc2_data.lane.t - sc2_data.gnssppk.t(i)).^2);
    sc2_data.gnssppk.idxs(i) = idx;
end
%% Check liveLocationKalman accuracy
ref = ppk_data.ref;
doi = sc2_data;
N = length(doi.imu.lat);
N2 = length(doi.gnssppk.lat);
pos = zeros(N,2);
pos_ppk = zeros(N2,2);
for i=1:N
    pos(i,:) = geo_to_lin(doi.imu.lat(i),doi.imu.lon(i),ref);
end

for i=1:N2
    pos_ppk(i,:) = geo_to_lin(doi.gnssppk.lat(i),doi.gnssppk.lon(i),ref);
end
figure(1);
plot(pos(120:end,1),pos(120:end,2),'b.'); hold on; axis equal; grid on; 
plot(pos_ppk(:,1),pos_ppk(:,2),'r.');

for i=1:100:length(sc2_data.gnssppk.idxs)
    plot([pos_ppk(i,1) pos(sc2_data.gnssppk.idxs(i),1)],[pos_ppk(i,2) pos(sc2_data.gnssppk.idxs(i),2)],'k--');
end
legend('liveLocationKalman','PPK')

figure(2);
plot(doi.imu.lat,'r-'); hold on; grid on;
plot(doi.imu.lon,'b-');
legend('Latitude','Longitude')

figure(3);
plot(doi.gnssppk.cov(1,:))
%% Concatenation Function
function S = CatStructFields(dim, varargin)
    F = cellfun(@fieldnames,varargin,'uni',0);
    assert(isequal(F{:}),'All structures must have the same field names.')
    T = [varargin{:}];
    S = struct();
    F = F{1};
    for k = 1:numel(F)
        S.(F{k}) = cat(dim,T.(F{k}));
    end
end
%% Field Transpose
function S = FullTranspose(F)
    S = struct();
    fields = fieldnames(F);
    
    for i=1:length(fields)
        name = fields{i};
        S.(name) = F.(name)';
    end
end

%% Latitude/Longitude to TM Coord. Converter

function lpos=geo_to_lin(PosLat, PosLon, PosRef)

% PosLat : Latitude in Degrees
% PosLon : Longitude in Degrees
% PosRef : Origin of Local coordinate
% lpos : [ localX, localY ] in meter, East-North Coordinate

% Convert Geographic coordinate into Linear coordinate with WGS84
% Ellipsoid model constants (actual values here are for WGS84

R0 = 6378137.0;
E=1/298.257223563;

deltaLon = PosLon - PosRef(2);
deltaLat = PosLat - PosRef(1);

Rn = R0*(1-E^2)/((1-E^2 * (sind(PosRef(1))^2))^(3/2));
Re = R0/((1-E^2 *(sind(PosRef(1))^2))^(1/2));
Ra = Re*Rn / sqrt( Re^2*sind(PosRef(1))^2 + Rn^2*cosd(PosRef(1))^2 );

localX = Ra * cosd(PosRef(1)) * deltaLon * pi/180;
localY = Rn * deltaLat * pi/180;

lpos = [ localX, localY ];

% display('geo_to_lin.m >> DONE')
end
