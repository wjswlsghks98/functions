%%
clear; close all; clc;
sc1 = load('/home/acl/JINHWAN/20220803SJDATA/ppk/20220803_021724.mat');
sc2 = load('/home/acl/JINHWAN/20220803SJDATA/ppk/20220803_023138.mat');
sc3 = load('/home/acl/JINHWAN/20220803SJDATA/ppk/20220803_024836.mat');
sc4 = load('/home/acl/JINHWAN/20220803SJDATA/ppk/20220803_032522.mat');
sc5 = load('/home/acl/JINHWAN/20220803SJDATA/ppk/20220803_035409.mat');

sc1 = table2cell(sc1.sc1.data);
sc2 = table2cell(sc2.sc2.data);
sc3 = table2cell(sc3.sc3.data);
sc4 = table2cell(sc4.sc4.data);
sc5 = table2cell(sc5.sc5.data);

%% 
[m,~] = size(sc1);
for i=1:m
    hms_ = split(sc1{i,2},":");
    hms = double(hms_(strlength(hms_)>0));
    sc1{i,16} = posixtime(datetime(2022,8,3,hms(1),hms(2),hms(3))) + sc1{i,14}; 
end
[m,~] = size(sc2);
for i=1:m
    hms_ = split(sc2{i,2},":");
    hms = double(hms_(strlength(hms_)>0));
    sc2{i,16} = posixtime(datetime(2022,8,3,hms(1),hms(2),hms(3))) + sc2{i,14}; 
end
[m,~] = size(sc3);
for i=1:m
    hms_ = split(sc3{i,2},":");
    hms = double(hms_(strlength(hms_)>0));
    sc3{i,16} = posixtime(datetime(2022,8,3,hms(1),hms(2),hms(3))) + sc3{i,14}; 
end
[m,~] = size(sc4);
for i=1:m
    hms_ = split(sc4{i,2},":");
    hms = double(hms_(strlength(hms_)>0));
    sc4{i,16} = posixtime(datetime(2022,8,3,hms(1),hms(2),hms(3))) + sc4{i,14}; 
end
[m,~] = size(sc5);
for i=1:m
    hms_ = split(sc5{i,2},":");
    hms = double(hms_(strlength(hms_)>0));
    sc5{i,16} = posixtime(datetime(2022,8,3,hms(1),hms(2),hms(3))) + sc5{i,14}; 
end

%% Save GNSS data
ppk_data = struct();
ppk_data.ref = [36.522213381, 127.303189872]; % at readme.txt (POS file header)
ppk_data.sc1 = struct();
ppk_data.sc2 = struct();
ppk_data.sc3 = struct();
ppk_data.sc4 = struct();
ppk_data.sc5 = struct();

%%

% GNSS timestamp (unix)
ppk_data.sc1.t = cell2mat(sc1(:,16));
ppk_data.sc2.t = cell2mat(sc2(:,16));
ppk_data.sc3.t = cell2mat(sc3(:,16));
ppk_data.sc4.t = cell2mat(sc4(:,16));
ppk_data.sc5.t = cell2mat(sc5(:,16));

% GNSS lat,lon,alt
ppk_data.sc1.lat = cell2mat(sc1(:,3));
ppk_data.sc1.lon = cell2mat(sc1(:,4));
ppk_data.sc1.alt = cell2mat(sc1(:,5));
ppk_data.sc1.num_sat = cell2mat(sc1(:,7)); % Number of satelites

ppk_data.sc2.lat = cell2mat(sc2(:,3));
ppk_data.sc2.lon = cell2mat(sc2(:,4));
ppk_data.sc2.alt = cell2mat(sc2(:,5));
ppk_data.sc2.num_sat = cell2mat(sc2(:,7)); % Number of satelites

ppk_data.sc3.lat = cell2mat(sc3(:,3));
ppk_data.sc3.lon = cell2mat(sc3(:,4));
ppk_data.sc3.alt = cell2mat(sc3(:,5));
ppk_data.sc3.num_sat = cell2mat(sc3(:,7)); % Number of satelites

ppk_data.sc4.lat = cell2mat(sc4(:,3));
ppk_data.sc4.lon = cell2mat(sc4(:,4));
ppk_data.sc4.alt = cell2mat(sc4(:,5));
ppk_data.sc4.num_sat = cell2mat(sc4(:,7)); % Number of satelites

ppk_data.sc5.lat = cell2mat(sc5(:,3));
ppk_data.sc5.lon = cell2mat(sc5(:,4));
ppk_data.sc5.alt = cell2mat(sc5(:,5));
ppk_data.sc5.num_sat = cell2mat(sc5(:,7)); % Number of satelites

% GNSS covariance 
% N E  U NE EU UN
% y x  z xy xz yz
% 8 9 10 11 12 13

% sc1
m = size(sc1,1);
ppk_data.sc1.cov = zeros(9,m);
ppk_data.sc1.cov(1,:) = cell2mat(sc1(:,9)).^2';
ppk_data.sc1.cov(2,:) = cell2mat(sc1(:,11)).^2';
ppk_data.sc1.cov(3,:) = cell2mat(sc1(:,12)).^2';

ppk_data.sc1.cov(5,:) = cell2mat(sc1(:,8)).^2';
ppk_data.sc1.cov(6,:) = cell2mat(sc1(:,13)).^2';

ppk_data.sc1.cov(9,:) = cell2mat(sc1(:,10)).^2';

negidxs2 = cell2mat(sc1(:,11)) < 0;
negidxs3 = cell2mat(sc1(:,12)) < 0;
negidxs6 = cell2mat(sc1(:,13)) < 0;

ppk_data.sc1.cov(2,negidxs2) = -ppk_data.sc1.cov(2,negidxs2);
ppk_data.sc1.cov(3,negidxs3) = -ppk_data.sc1.cov(3,negidxs3);
ppk_data.sc1.cov(6,negidxs6) = -ppk_data.sc1.cov(6,negidxs6);

ppk_data.sc1.cov(4,:) = ppk_data.sc1.cov(2,:);
ppk_data.sc1.cov(7,:) = ppk_data.sc1.cov(3,:);
ppk_data.sc1.cov(8,:) = ppk_data.sc1.cov(6,:);

% sc2
m = size(sc2,1);
ppk_data.sc2.cov = zeros(9,m);
ppk_data.sc2.cov(1,:) = cell2mat(sc2(:,9)).^2';
ppk_data.sc2.cov(2,:) = cell2mat(sc2(:,11)).^2';
ppk_data.sc2.cov(3,:) = cell2mat(sc2(:,12)).^2';

ppk_data.sc2.cov(5,:) = cell2mat(sc2(:,8)).^2';
ppk_data.sc2.cov(6,:) = cell2mat(sc2(:,13)).^2';

ppk_data.sc2.cov(9,:) = cell2mat(sc2(:,10)).^2';

negidxs2 = cell2mat(sc2(:,11)) < 0;
negidxs3 = cell2mat(sc2(:,12)) < 0;
negidxs6 = cell2mat(sc2(:,13)) < 0;

ppk_data.sc2.cov(2,negidxs2) = -ppk_data.sc2.cov(2,negidxs2);
ppk_data.sc2.cov(3,negidxs3) = -ppk_data.sc2.cov(3,negidxs3);
ppk_data.sc2.cov(6,negidxs6) = -ppk_data.sc2.cov(6,negidxs6);

ppk_data.sc2.cov(4,:) = ppk_data.sc2.cov(2,:);
ppk_data.sc2.cov(7,:) = ppk_data.sc2.cov(3,:);
ppk_data.sc2.cov(8,:) = ppk_data.sc2.cov(6,:);

% sc3
m = size(sc3,1);
ppk_data.sc3.cov = zeros(9,m);
ppk_data.sc3.cov(1,:) = cell2mat(sc3(:,9)).^2';
ppk_data.sc3.cov(2,:) = cell2mat(sc3(:,11)).^2';
ppk_data.sc3.cov(3,:) = cell2mat(sc3(:,12)).^2';

ppk_data.sc3.cov(5,:) = cell2mat(sc3(:,8)).^2';
ppk_data.sc3.cov(6,:) = cell2mat(sc3(:,13)).^2';

ppk_data.sc3.cov(9,:) = cell2mat(sc3(:,10)).^2';

negidxs2 = cell2mat(sc3(:,11)) < 0;
negidxs3 = cell2mat(sc3(:,12)) < 0;
negidxs6 = cell2mat(sc3(:,13)) < 0;

ppk_data.sc3.cov(2,negidxs2) = -ppk_data.sc3.cov(2,negidxs2);
ppk_data.sc3.cov(3,negidxs3) = -ppk_data.sc3.cov(3,negidxs3);
ppk_data.sc3.cov(6,negidxs6) = -ppk_data.sc3.cov(6,negidxs6);

ppk_data.sc3.cov(4,:) = ppk_data.sc3.cov(2,:);
ppk_data.sc3.cov(7,:) = ppk_data.sc3.cov(3,:);
ppk_data.sc3.cov(8,:) = ppk_data.sc3.cov(6,:);

% sc4
m = size(sc4,1);
ppk_data.sc4.cov = zeros(9,m);
ppk_data.sc4.cov(1,:) = cell2mat(sc4(:,9)).^2';
ppk_data.sc4.cov(2,:) = cell2mat(sc4(:,11)).^2';
ppk_data.sc4.cov(3,:) = cell2mat(sc4(:,12)).^2';

ppk_data.sc4.cov(5,:) = cell2mat(sc4(:,8)).^2';
ppk_data.sc4.cov(6,:) = cell2mat(sc4(:,13)).^2';

ppk_data.sc4.cov(9,:) = cell2mat(sc4(:,10)).^2';

negidxs2 = cell2mat(sc4(:,11)) < 0;
negidxs3 = cell2mat(sc4(:,12)) < 0;
negidxs6 = cell2mat(sc4(:,13)) < 0;

ppk_data.sc4.cov(2,negidxs2) = -ppk_data.sc4.cov(2,negidxs2);
ppk_data.sc4.cov(3,negidxs3) = -ppk_data.sc4.cov(3,negidxs3);
ppk_data.sc4.cov(6,negidxs6) = -ppk_data.sc4.cov(6,negidxs6);

ppk_data.sc4.cov(4,:) = ppk_data.sc4.cov(2,:);
ppk_data.sc4.cov(7,:) = ppk_data.sc4.cov(3,:);
ppk_data.sc4.cov(8,:) = ppk_data.sc4.cov(6,:);

% sc5

m = size(sc5,1);
ppk_data.sc5.cov = zeros(9,m);
ppk_data.sc5.cov(1,:) = cell2mat(sc5(:,9)).^2';
ppk_data.sc5.cov(2,:) = cell2mat(sc5(:,11)).^2';
ppk_data.sc5.cov(3,:) = cell2mat(sc5(:,12)).^2';

ppk_data.sc5.cov(5,:) = cell2mat(sc5(:,8)).^2';
ppk_data.sc5.cov(6,:) = cell2mat(sc5(:,13)).^2';

ppk_data.sc5.cov(9,:) = cell2mat(sc5(:,10)).^2';

negidxs2 = cell2mat(sc5(:,11)) < 0;
negidxs3 = cell2mat(sc5(:,12)) < 0;
negidxs6 = cell2mat(sc5(:,13)) < 0;

ppk_data.sc5.cov(2,negidxs2) = -ppk_data.sc5.cov(2,negidxs2);
ppk_data.sc5.cov(3,negidxs3) = -ppk_data.sc5.cov(3,negidxs3);
ppk_data.sc5.cov(6,negidxs6) = -ppk_data.sc5.cov(6,negidxs6);

ppk_data.sc5.cov(4,:) = ppk_data.sc5.cov(2,:);
ppk_data.sc5.cov(7,:) = ppk_data.sc5.cov(3,:);
ppk_data.sc5.cov(8,:) = ppk_data.sc5.cov(6,:);

