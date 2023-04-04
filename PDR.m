clc;clear;close all
load('matlab.mat');
%% 传感器数据
figure

t=data(:,1);
subplot(3,1,1)
hold on
grid on
box on
plot(t,data(:,2),Color='b',LineStyle='--',LineWidth=1)
plot(t,data(:,3),Color=[0 0.40 0],LineStyle=':',LineWidth=2)
plot(t,data(:,4),Color='r',LineStyle='-',LineWidth=1)
title("Gyro data",FontSize=20)
legend("Gx","Gy","Gz",Location="bestoutside",fontsize=20)
ylabel("(deg/s)",FontSize=20)
subplot(3,1,2)
hold on
grid on
box on
plot(t,data(:,5),Color='b',LineStyle='--',LineWidth=1)
plot(t,data(:,6),Color=[0 0.40 0],LineStyle=':',LineWidth=2)
plot(t,data(:,7),Color='r',LineStyle='-',LineWidth=1)
title("Accelerometer data",FontSize=20)
legend("Ax","Ay","Az",Location="bestoutside",fontsize=20)
ylabel("(m/s^2)",FontSize=20)
subplot(3,1,3)
hold on
grid on
box on
plot(t,data(:,8),Color='b',LineStyle='--',LineWidth=1)
plot(t,data(:,9),Color=[0 0.40 0],LineStyle=':',LineWidth=2)
plot(t,data(:,10),Color='r',LineStyle='-',LineWidth=1)
title("Magnetometers",FontSize=20)
legend("Mx","My","Mz",Location="bestoutside",fontsize=20)
ylabel("(Gauss)",FontSize=20)

%% 水平姿态角
roll=zeros(6829,1);
pitch=zeros(6829,1);
for i=1:6829
    roll(i)=atan2(-data(i,6),-data(i,7));
    pitch(i)=atan2(data(i,5),sqrt(data(i,6)^2+data(i,7)^2));
end
figure
box on
hold on
grid on
plot(data(:,1),roll*180/pi,LineStyle="-",Color="b",LineWidth=2);
plot(data(:,1),pitch*180/pi,LineStyle="--",Color="r",LineWidth=2);
title("Horizontal angles,without smoothing",fontsize=20)
legend("Roll","Pitch",fontsize=20,Location="bestoutside")
xlabel("Time(s)",fontsize=20)
ylabel("(deg)",fontsize=20)

%% 水平姿态角平滑
ave_roll=movmean(roll,20);
ave_pitch=movmean(pitch,20);
figure
box on
hold on
grid on
plot(data(:,1),ave_roll*180/pi,LineStyle="-",Color="b",LineWidth=3);
plot(data(:,1),ave_pitch*180/pi,LineStyle="--",Color="r",LineWidth=3);
title("Horizontal angles",fontsize=20)
legend("Roll","Pitch",fontsize=20,Location="bestoutside")
xlabel("Time(s)",fontsize=20)
ylabel("(deg)",fontsize=20)


%% 磁强计数据计算航向角
Psi=zeros(6829,1);
D=14.11*pi/180;
for i=1:6829
    m_x_n1=data(i,8)*cos(ave_pitch(i))+data(i,9)*sin(ave_pitch(i))*sin(ave_roll(i))+data(i,10)*cos(ave_roll(i))*sin(ave_pitch(i));
    m_y_n1=data(i,9)*cos(ave_roll(i))-data(i,10)*sin(ave_roll(i));
    Psi(i)=-atan2(m_y_n1,m_x_n1)+D;
    Psi(i)=Psi(i)*180/pi;
    if(Psi(i)>=180)
        Psi(i)=Psi(i)-360;
    elseif(Psi(i)<-180)
        Psi(i)=Psi(i)+360;
    end
end
figure
grid on
box on
plot(data(:,1),Psi,LineStyle="-",Color="b",LineWidth=3);
title("Heading from magnetometers",fontsize=20)
legend("\Psi_m",fontsize=20,Location="bestoutside")
xlabel("Time(s)",fontsize=20)
ylabel("(deg)",fontsize=20)

%% 陀螺数据计算航向角
errornum=0.3;%简单Z轴常值零偏补偿(deg/s)
data(:,4)=data(:,4)+errornum;
WD=zeros(6829,1);
Psi_g=-90*pi/180;
for i=1:6829
    WD(i)=Psi_g+[-sin(ave_roll(i)) sin(ave_roll(i))*cos(ave_pitch(i)) cos(ave_roll(i))*cos(ave_pitch(i))]*data(i,2:4)'*pi/180*0.05;
    WD(i)=WD(i)*180/pi;
    if(WD(i)>=180)
        WD(i)=WD(i)-360;
    elseif(WD(i)<-180)
        WD(i)=WD(i)+360;
    end
    Psi_g=WD(i,1)*pi/180;
    
end
figure
box on
grid on
hold on
plot(data(:,1),Psi,LineStyle="-",Color="b",LineWidth=3);
plot(data(:,1),WD,LineStyle="--",color="r",LineWidth=3);
title("Heading from magnetometers and gyro",fontsize=20)
legend("\Psi_m","\Psi_g",fontsize=20,Location="bestoutside")
xlabel("Time(s)",fontsize=20)
ylabel("(deg)",fontsize=20)

%% 脚步探测
num=9.2;
step=[];
accNorms=zeros(6829,1);
for i=1:6829
    accNorms(i,1)=norm(data(i,5:7));
end
mid=-accNorms;
[~,locs]=findpeaks(mid);
for i=1:length(locs)
    if accNorms(locs(i))<num
        step=[step;[locs(i)*0.05,accNorms(locs(i))]];
    end
end
figure
grid on
box on
hold on
plot(data(:,1),accNorms(:,1),LineStyle="-",Color="b",LineWidth=3);
scatter(step(:,1),step(:,2),"r",Marker="+",SizeData=100,LineWidth=2);
title("Steps",fontsize=20)
legend("Accel","Step",fontsize=20,Location="bestoutside")
xlabel("Time(s)",fontsize=20)
ylabel("Acceleration",fontsize=20)

%% PDR解算
steplong=0.7;
ENloca=zeros(length(step(:,1)),2);
for i=1:length(step(:,1))
    loca=int64(step(:,1)/0.05);
    j=(loca(i));
    ENloca(i+1,1)=ENloca(i,1)+steplong*sin(Psi(j)*pi/180);
    ENloca(i+1,2)=ENloca(i,2)+steplong*cos(Psi(j)*pi/180);
end
ENloca_G=zeros(length(step(:,1)),2);
for i=1:length(step(:,1))
    loca=int64(step(:,1)/0.05);
    j=(loca(i));
    ENloca_G(i+1,1)=ENloca_G(i,1)+steplong*sin(WD(j)*pi/180);
    ENloca_G(i+1,2)=ENloca_G(i,2)+steplong*cos(WD(j)*pi/180);
end
figure
hold on
box on
plot(ENloca(:,1),ENloca(:,2),LineStyle="-",Color="b",LineWidth=2);
plot(ENloca_G(:,1),ENloca_G(:,2),LineStyle="--",Color="r",LineWidth=2)
set(gca,'xtick',-200:50:100)
set(gca,'ytick',-200:50:50)
title("PDR using Heading from magnetometers and gyro",fontsize=20)
legend("P_m","P_g")
xlabel("East(m)",fontsize=20)
ylabel("North(m)",fontsize=20)
grid on
axis equal


% %% 对陀螺Z轴加入常值零偏补偿
% data2=data;
% errorNum=0.2;
% data2(:,4)=data2(:,4)-errorNum;
% % Z轴简单修正陀螺数据计算航向角
% WD2=zeros(6829,1);
% Psi_g=-90*pi/180;
% for i=1:6829
%     WD2(i)=Psi_g+[-sin(ave_roll(i)) sin(ave_roll(i))*cos(ave_pitch(i)) cos(ave_roll(i))*cos(ave_pitch(i))]*data2(i,2:4)'*pi/180*0.05;
%     WD2(i)=WD2(i)*180/pi;
%     if(WD2(i)>=180)
%         WD2(i)=WD2(i)-360;
%     elseif(WD2(i)<-180)
%         WD2(i)=WD2(i)+360;
%     end
%     Psi_g=WD2(i,1)*pi/180;
%     
% end
% figure
% grid on
% hold on
% plot(data(:,1),Psi,LineStyle="-",Color="b",LineWidth=2);
% plot(data(:,1),WD2,LineStyle="--",color="r",LineWidth=2);
% title("Heading from magnetometers and gyro",fontsize=20)
% legend("\Psi_m","\Psi_g",fontsize=20,Location="bestoutside")
% xlabel("Time(s)",fontsize=20)
% ylabel("(deg)",fontsize=20)
% % PDR解算
% steplong=0.7;
% ENloca=zeros(length(step(:,1)),2);
% for i=1:length(step(:,1))
%     loca=int64(step(:,1)/0.05);
%     j=(loca(i));
%     ENloca(i+1,1)=ENloca(i,1)+steplong*sin(Psi(j)*pi/180);
%     ENloca(i+1,2)=ENloca(i,2)+steplong*cos(Psi(j)*pi/180);
% end
% ENloca_G2=zeros(length(step(:,1)),2);
% for i=1:length(step(:,1))
%     loca=int64(step(:,1)/0.05);
%     j=(loca(i));
%     ENloca_G2(i+1,1)=ENloca_G2(i,1)+steplong*sin(WD2(j)*pi/180);
%     ENloca_G2(i+1,2)=ENloca_G2(i,2)+steplong*cos(WD2(j)*pi/180);
% end
% figure
% hold on
% plot(ENloca(:,1),ENloca(:,2),LineStyle="-",Color="b",LineWidth=2);
% plot(ENloca_G(:,1),ENloca_G2(:,2),LineStyle="--",Color="r",LineWidth=2)
% set(gca,'xtick',-200:50:100)
% set(gca,'ytick',-200:50:50)
% title("PDR using Heading from magnetometers and gyro",fontsize=20)
% legend("P_m","P_g")
% xlabel("East(m)")
% ylabel("North(m)")
% grid on
% axis equal