

bag = rosbag("D:\2021-03-21-13-06-01.bag");


output = struct();

%% 

temptopic = "/ublox/navpvt"
tempbag = select(bag, 'Topic' , temptopic);
tempmsg = readMessages(tempbag, 'DataFormat', 'struct');
templen = length(tempmsg);
msgname = 'ubloxvel';
output.(msgname) = struct();
output.(msgname).vx = zeros(1,templen);
output.(msgname).vy = zeros(1,templen);
output.(msgname).vz = zeros(1,templen);

for i = 1:length(tempmsg)
    msg = tempmsg{i};
    output.(msgname).vx(i) = msg.VelE;
    output.(msgname).vy(i) = msg.VelN;
    output.(msgname).vz(i) = msg.VelD;
end
%% 

temptopic = "/ublox/fix"
tempbag = select(bag, 'Topic' , temptopic);
tempmsg = readMessages(tempbag, 'DataFormat', 'struct');
templen = length(tempmsg);
msgname = 'ubloxgps';
output.(msgname) = struct();
output.(msgname).t = zeros(1,templen);
output.(msgname).lat = zeros(1,templen);
output.(msgname).lon = zeros(1,templen);
output.(msgname).alt = zeros(1,templen);

for i = 1:length(tempmsg)
    msg = tempmsg{i};
    output.(msgname).t(i) = double(msg.Header.Stamp.Sec) + double(msg.Header.Stamp.Nsec)*1e-9;
    output.(msgname).lat(i) = msg.Latitude;
    output.(msgname).lon(i) = msg.Longitude;
    output.(msgname).alt(i) = msg.Altitude;
end


%%

temptopic = "/ublox/fix_velocity"
tempbag = select(bag, 'Topic' , temptopic);
tempmsg = readMessages(tempbag, 'DataFormat', 'struct');
templen = length(tempmsg);
msgname = 'ubloxvel';
output.(msgname) = struct();
output.(msgname).vx = zeros(1,templen);
output.(msgname).vy = zeros(1,templen);
output.(msgname).vz = zeros(1,templen);
output.(msgname).vxstd = zeros(1,templen);
output.(msgname).vystd = zeros(1,templen);
output.(msgname).vzstd = zeros(1,templen);

for i = 1:length(tempmsg)
    msg = tempmsg{i};
    output.(msgname).t(i) = double(msg.Header.Stamp.Sec) + double(msg.Header.Stamp.Nsec)*1e-9;
    output.(msgname).vx(i) = msg.Twist.Twist.Linear.X;
    output.(msgname).vy(i) = msg.Twist.Twist.Linear.Y;
    output.(msgname).vz(i) = msg.Twist.Twist.Linear.Z;
    output.(msgname).vxstd(i) = msg.Twist.Covariance(1);
    output.(msgname).vystd(i) = msg.Twist.Covariance(8);
    output.(msgname).vzstd(i) = msg.Twist.Covariance(15);
end

%% 

temptopic = "/RT/gps/vel"
tempbag = select(bag, 'Topic' , temptopic);
tempmsg = readMessages(tempbag, 'DataFormat', 'struct');
templen = length(tempmsg);
msgname = 'rtvel';
output.(msgname) = struct();
output.(msgname).t = zeros(1,templen);
output.(msgname).vx = zeros(1,templen);
output.(msgname).vy = zeros(1,templen);
output.(msgname).vz = zeros(1,templen);

for i = 1:length(tempmsg)
    msg = tempmsg{i};
    output.(msgname).t(i) = double(msg.Header.Stamp.Sec) + double(msg.Header.Stamp.Nsec)*1e-9;
    output.(msgname).vx(i) = msg.Twist.Twist.Linear.X;
    output.(msgname).vy(i) = msg.Twist.Twist.Linear.Y;
    output.(msgname).vz(i) = msg.Twist.Twist.Linear.Z;
end

%% 

temptopic = "/RT/gps/fix"
tempbag = select(bag, 'Topic' , temptopic);
tempmsg = readMessages(tempbag, 'DataFormat', 'struct');
templen = length(tempmsg);
msgname = 'rtgps';
output.(msgname) = struct();
output.(msgname).t = zeros(1,templen);
output.(msgname).lat = zeros(1,templen);
output.(msgname).lon = zeros(1,templen);
output.(msgname).alt = zeros(1,templen);

for i = 1:length(tempmsg)
    msg = tempmsg{i};
    output.(msgname).t(i) = double(msg.Header.Stamp.Sec) + double(msg.Header.Stamp.Nsec)*1e-9;
    output.(msgname).lat(i) = msg.Latitude;
    output.(msgname).lon(i) = msg.Longitude;
    output.(msgname).alt(i) = msg.Altitude;
end

%%

temptopic = "/RT/imu/data"
tempbag = select(bag, 'Topic' , temptopic);
tempmsg = readMessages(tempbag, 'DataFormat', 'struct');
templen = length(tempmsg);
msgname = 'rtimu';
output.(msgname) = struct();
output.(msgname).t = zeros(1,templen);
output.(msgname).ax = zeros(1,templen);
output.(msgname).ay = zeros(1,templen);
output.(msgname).az = zeros(1,templen);
output.(msgname).wx = zeros(1,templen);
output.(msgname).wy = zeros(1,templen);
output.(msgname).wz = zeros(1,templen);
output.(msgname).rx = zeros(1,templen);
output.(msgname).ry = zeros(1,templen);
output.(msgname).rz = zeros(1,templen);

for i = 1:length(tempmsg)
    msg = tempmsg{i};
    output.(msgname).t(i) = double(msg.Header.Stamp.Sec) + double(msg.Header.Stamp.Nsec)*1e-9;
    output.(msgname).ax(i) = msg.LinearAcceleration.X;
    output.(msgname).ay(i) = -msg.LinearAcceleration.Y;
    output.(msgname).az(i) = msg.LinearAcceleration.Z;
    output.(msgname).wx(i) = msg.AngularVelocity.X;
    output.(msgname).wy(i) = -msg.AngularVelocity.Y;
    output.(msgname).wz(i) = msg.AngularVelocity.Z;
    temp = quat2eul([msg.Orientation.W, msg.Orientation.X, msg.Orientation.Y, msg.Orientation.Z]);
    output.(msgname).rx(i) = temp(3);
    output.(msgname).ry(i) = -temp(2);
    output.(msgname).rz(i) = temp(1);
end

% %%
% 
% temptopic = "/openpilot/carState"
% tempbag = select(bag, 'Topic' , temptopic);
% tempmsg = readMessages(tempbag, 'DataFormat', 'struct');
% templen = length(tempmsg);
% msgname = 'opcarstate';
% output.(msgname) = struct();
% output.(msgname).t = zeros(1,templen);
% output.(msgname).delta = zeros(1,templen);
% output.(msgname).wsf = zeros(1,templen);
% output.(msgname).wsr = zeros(1,templen);
% output.(msgname).v = zeros(1,templen);
% for i = 1:length(tempmsg)
%     msg = tempmsg{i};
%     output.(msgname).t(i) = double(msg.Header.Stamp.Sec) + double(msg.Header.Stamp.Nsec)*1e-9;
%     output.(msgname).delta(i) = msg.SteeringAngle.Data;
%     output.(msgname).wsf(i) = (msg.WheelSpeeds.Fl.Data+msg.WheelSpeeds.Fr.Data)/2;
%     output.(msgname).wsr(i) = (msg.WheelSpeeds.Rl.Data+msg.WheelSpeeds.Rr.Data)/2;
%     output.(msgname).v(i) = msg.VEgo.Data;
% end
% 
%%

temptopic = "/xsens/filter/quaternion"
tempbag = select(bag, 'Topic' , temptopic);
tempmsg = readMessages(tempbag, 'DataFormat', 'struct');
templen = length(tempmsg);
msgname = 'xs';
output.(msgname) = struct();
output.(msgname).t = zeros(1,templen);
output.(msgname).rx = zeros(1,templen);
output.(msgname).ry = zeros(1,templen);
output.(msgname).rz = zeros(1,templen);
output.(msgname).vx = zeros(1,templen);
output.(msgname).vy = zeros(1,templen);
output.(msgname).vz = zeros(1,templen);
output.(msgname).ax = zeros(1,templen);
output.(msgname).ay = zeros(1,templen);
output.(msgname).az = zeros(1,templen);
output.(msgname).wx = zeros(1,templen);
output.(msgname).wy = zeros(1,templen);
output.(msgname).wz = zeros(1,templen);

for i = 1:length(tempmsg)
    msg = tempmsg{i};
    output.(msgname).t(i) = double(msg.Header.Stamp.Sec) + double(msg.Header.Stamp.Nsec)*1e-9;
    temp = quat2eul([msg.Quaternion.W, msg.Quaternion.X, msg.Quaternion.Y, msg.Quaternion.Z]);
    output.(msgname).rx(i) = temp(3);
    output.(msgname).ry(i) = temp(2);
    output.(msgname).rz(i) = temp(1);
end

temptopic = "/xsens/filter/velocity"
tempbag = select(bag, 'Topic' , temptopic);
tempmsg = readMessages(tempbag, 'DataFormat', 'struct');
for i = 1:length(tempmsg)
    msg = tempmsg{i};
    output.(msgname).vx(i) = msg.Vector.X;
    output.(msgname).vy(i) = msg.Vector.Y;
    output.(msgname).vz(i) = msg.Vector.Z;
end

temptopic = "/xsens/imu/acceleration"
tempbag = select(bag, 'Topic' , temptopic);
tempmsg = readMessages(tempbag, 'DataFormat', 'struct');
for i = 1:length(tempmsg)
    msg = tempmsg{i};
    output.(msgname).ax(i) = msg.Vector.X;
    output.(msgname).ay(i) = msg.Vector.Y;
    output.(msgname).az(i) = msg.Vector.Z;
end

temptopic = "/xsens/imu/angular_velocity"
tempbag = select(bag, 'Topic' , temptopic);
tempmsg = readMessages(tempbag, 'DataFormat', 'struct');
for i = 1:length(tempmsg)
    msg = tempmsg{i};
    output.(msgname).wx(i) = msg.Vector.X;
    output.(msgname).wy(i) = msg.Vector.Y;
    output.(msgname).wz(i) = msg.Vector.Z;
end

%%

temptopic = "cam0/pub/image/compressed"
tempbag = select(bag,'Topic', temptopic);
tempmsg = readMessages(tempbag, 'DataFormat', 'struct');
templen = length(tempmsg);
msgname = 'cam';
output.(msgname) = struct();
output.(msgname).t = zeros(1,templen-1);
for i = 2:length(tempmsg)
    msg = tempmsg{i};
    output.(msgname).t(i-1) = double(msg.Header.Stamp.Sec) + double(msg.Header.Stamp.Nsec)*1e-9;
end

%%

output.cam.x = [0,0.1875,0.7500,1.6875,3,4.6875,6.750,9.1875,12,15.1875,18.75,22.6875,27,31.6875,36.75,42.1875,48,54.1875,60.75,67.6875,75,82.6875,90.75,99.1875,108,117.1875,126.75,136.6875,147,157.6875,168.75,180.1875,192];
output.cam.l = -csvread('ly.csv');
output.cam.ll = -csvread('lly.csv');
output.cam.le = -csvread('ley.csv');
output.cam.r = -csvread('ry.csv');
output.cam.rr = -csvread('rry.csv');
output.cam.re = -csvread('rey.csv');

% Data Interpolation
output.cam.x_inter = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100];
output.cam.l_inter = zeros(size(output.cam.l,1),11);
output.cam.r_inter = zeros(size(output.cam.r,1),11);

output.cam.l_inter(:,1) = output.cam.l(:,1); 
output.cam.r_inter(:,1) = output.cam.r(:,1); 
% for i=1:10
%     j=1;
%     while output.cam.x(j) < output.cam.x_inter(i+1)
%         j = j+1;
%     end
%     output.cam.l_inter(:,i+1) = output.cam.l(:,j-1) + (output.cam.l(:,j) - output.cam.l(:,j-1))./(output.cam.x(j) - output.cam.x(j-1)) .* (output.cam.x_inter(i+1) - output.cam.x(j-1));
%     output.cam.r_inter(:,i+1) = output.cam.r(:,j-1) + (output.cam.r(:,j) - output.cam.r(:,j-1))./(output.cam.x(j) - output.cam.x(j-1)) .* (output.cam.x_inter(i+1) - output.cam.x(j-1));
% end

for i = 1:numel(output.cam.t)
    output.cam.l_inter(i,:) = interp1(output.cam.x, output.cam.l(i,:), output.cam.x_inter);
    output.cam.r_inter(i,:) = interp1(output.cam.x, output.cam.r(i,:), output.cam.x_inter);
end



tempmsg = csvread('pose.csv');
output.cam.vx = tempmsg(:,1);
output.cam.vy = -tempmsg(:,2);
output.cam.vz = -tempmsg(:,3);
output.cam.wz = tempmsg(:,6);
output.cam.wy = -tempmsg(:,5);
output.cam.wx = -tempmsg(:,4);

lane = struct();
lane.lx = zeros(size(output.cam.l_inter));
lane.ly = zeros(size(output.cam.l_inter));
lane.rx = zeros(size(output.cam.l_inter));
lane.ry = zeros(size(output.cam.l_inter));
lane.posx = zeros(size(output.ubloxvel.t,2),1);
lane.posy = zeros(size(output.ubloxvel.t,2),1);
lane.posz = zeros(size(output.ubloxvel.t,2),1);
lane.imu.a = zeros(3,size(output.ubloxvel.t,2));
lane.imu.w = zeros(3,size(output.ubloxvel.t,2));
lane.xs.a = zeros(3,size(output.ubloxvel.t,2));
lane.xs.w = zeros(3,size(output.ubloxvel.t,2));
lane.xs.v = zeros(3,size(output.ubloxvel.t,2));
lane.xs.r = zeros(3,size(output.ubloxvel.t,2));
lane.xs.vcov = zeros(3,size(output.ubloxvel.t,2));
lane.xs.t = zeros(size(output.ubloxvel.t,2),1);
lane.eul.x = zeros(size(output.ubloxvel.t,2),1);
lane.eul.y = zeros(size(output.ubloxvel.t,2),1);
lane.eul.z = zeros(size(output.ubloxvel.t,2),1);
lane.vel.x = zeros(size(output.ubloxvel.t,2),1);
lane.vel.y = zeros(size(output.ubloxvel.t,2),1);
lane.vel.z = zeros(size(output.ubloxvel.t,2),1);
lane.t = zeros(size(output.ubloxvel.t,2),1);

%% Lane convert to global

lpos = geo_to_lin(output.rtgps.lat', output.rtgps.lon',[mean(output.rtgps.lat), mean(output.rtgps.lon)]);
lpos2 = geo_to_lin(output.ubloxgps.lat', output.ubloxgps.lon',[mean(output.rtgps.lat), mean(output.rtgps.lon)]);
output.ubloxgps.x = lpos2(:,1)';
output.ubloxgps.y = lpos2(:,2)';
output.rtgps.x = lpos(:,1)';
output.rtgps.y = lpos(:,2)';
for i = 1:length(output.cam.t)
[~,idx1] = min((output.cam.t(i) - output.rtgps.t).^2);
[~,idx2] = min((output.cam.t(i) - output.xs.t).^2);

lane.lx(i,:) = output.rtgps.x(idx1) + cos(output.rtimu.rz(idx1))*output.cam.x_inter - sin(output.rtimu.rz(idx1))*output.cam.l_inter(i,:);
lane.ly(i,:) = output.rtgps.y(idx1) + sin(output.rtimu.rz(idx1))*output.cam.x_inter + cos(output.rtimu.rz(idx1))*output.cam.l_inter(i,:);

lane.rx(i,:) = output.rtgps.x(idx1) + cos(output.rtimu.rz(idx1))*output.cam.x_inter - sin(output.rtimu.rz(idx1))*output.cam.r_inter(i,:);
lane.ry(i,:) = output.rtgps.y(idx1) + sin(output.rtimu.rz(idx1))*output.cam.x_inter + cos(output.rtimu.rz(idx1))*output.cam.r_inter(i,:);

lane.posx(i) = output.rtgps.x(idx1); lane.posy(i) = output.rtgps.y(idx1); lane.posz(i) = output.rtgps.alt(idx1);
lane.eul.x(i) = output.rtimu.rx(idx1); lane.eul.y(i) = output.rtimu.ry(idx1); lane.eul.z(i) = output.rtimu.rz(idx1);
lane.vel.x(i) = output.rtvel.vx(idx1); lane.vel.y(i) = output.rtvel.vy(idx1); lane.vel.z(i) =output.rtvel.vz(idx1);

lane.imu.a(:,i) = [output.rtimu.ax(idx1) output.rtimu.ay(idx1) output.rtimu.az(idx1)]' - [0 0 0]';
lane.imu.w(:,i) = [output.rtimu.wx(idx1) output.rtimu.wy(idx1) output.rtimu.wz(idx1)]' - [0 0 0]';

lane.xs.a(:,i) = [output.xs.ax(idx2) output.xs.ay(idx2) output.xs.az(idx2)]';
lane.xs.v(:,i) = [output.ubloxvel.vx(i) output.ubloxvel.vy(i) output.ubloxvel.vz(i)]';
lane.xs.w(:,i) = [output.xs.wx(idx2) output.xs.wy(idx2) output.xs.wz(idx2)]';
lane.xs.r(:,i) = [output.xs.rx(idx2) output.xs.ry(idx2) output.xs.rz(idx2)]';
lane.xs.vcov(:,i) = [output.ubloxvel.vxstd(i) output.ubloxvel.vystd(i) output.ubloxvel.vzstd(i)]';
lane.xs.t(i) = output.xs.t(idx2);
lane.t(i) = output.rtgps.t(idx1);
end