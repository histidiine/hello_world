function MJ_LT_Plots

MJ1 = load('Test_MJ_Aurelie.txt');
MJ2 = load('Elvira_MJ_1.txt');
MJ3 = load('Martina_MJ_1.txt');

LT1 = load('Elvira_LT_1.txt');
LT2 = load('Martina_LT_1.txt');
LT3 = load('Test_LT_Aurelie.txt');

VMJ1 = sqrt(MJ1(:,14).^2+MJ1(:,15).^2+MJ1(:,16).^2);
VMJ2 = sqrt(MJ2(:,14).^2+MJ2(:,15).^2+MJ2(:,16).^2);
VMJ3 = sqrt(MJ3(:,14).^2+MJ3(:,15).^2+MJ3(:,16).^2);

VLT1 = sqrt(LT1(:,14).^2+LT1(:,15).^2+LT1(:,16).^2);
VLT2 = sqrt(LT2(:,14).^2+LT2(:,15).^2+LT2(:,16).^2);
VLT3 = sqrt(LT3(:,14).^2+LT3(:,15).^2+LT3(:,16).^2);

VMJ1filt = sgolayfilt(VMJ1,1,21);
VMJ2filt = sgolayfilt(VMJ2,1,21);
VMJ3filt = sgolayfilt(VMJ3,1,21);

VLT1filt = sgolayfilt(VLT1,1,21);
VLT2filt = sgolayfilt(VLT2,1,21);
VLT3filt = sgolayfilt(VLT3,1,21);

timePoints=[];

figure;plot(VMJ1filt);
pause;
[x,y]= ginput(16); 
timePoints(:,end+1) = round(x);

figure;plot(VLT1filt);
pause;
[x,y]= ginput(16); 
timePoints(:,end+1) = round(x);


MJ_Mov1 = VMJ1filt(timePoints(1,1):timePoints(2,1)); %6
Target_Mov1 = MJ1(timePoints(1,1),33);
str1 = sprintf('Target %d',Target_Mov1);
figure;plot(MJ_Mov1),title(str1);
hold on;

LT_Mov3 = VLT1filt(timePoints(5,2):timePoints(6,2)); %6
Target_Mov3_LT = LT1(timePoints(5,2),33);
str3_LT = sprintf('Target %d',Target_Mov3_LT);
plot(LT_Mov3,'Color', 'r'),title(str3_LT);
legend('Minimum Jerk','Linear');

%

MJ_Mov2 = VMJ1filt(timePoints(3,1):timePoints(4,1)); %3
Target_Mov2 = MJ1(timePoints(3,1),33);
str2 = sprintf('Target %d',Target_Mov2);
figure;plot(MJ_Mov2),title(str2);
hold on;

LT_Mov8 = VLT1filt(timePoints(15,2):timePoints(16,2));%3
Target_Mov8_LT = LT1(timePoints(15,2),33);
str8_LT = sprintf('Target %d',Target_Mov8_LT);
plot(LT_Mov8,'Color', 'r'),title(str8_LT);
legend('Minimum Jerk','Linear');

% 

MJ_Mov3 = VMJ1filt(timePoints(5,1):timePoints(6,1)); %2
Target_Mov3 = MJ1(timePoints(5,1),33);
str3 = sprintf('Target %d',Target_Mov3);
figure;plot(MJ_Mov3),title(str3);
hold on;

LT_Mov5 = VLT1filt(timePoints(9,2):timePoints(10,2));%2
Target_Mov5_LT = LT1(timePoints(9,2),33);
str5_LT = sprintf('Target %d',Target_Mov5_LT);
plot(LT_Mov5,'Color', 'r'),title(str5_LT);
legend('Minimum Jerk','Linear');

%

MJ_Mov4 = VMJ1filt(timePoints(7,1):timePoints(8,1)); %8
Target_Mov4 = MJ1(timePoints(7,1),33);
str4 = sprintf('Target %d',Target_Mov4);
figure;plot(MJ_Mov4),title(str4);
hold on;

LT_Mov2 = VLT1filt(timePoints(3,2):timePoints(4,2)); %8
Target_Mov2_LT = LT1(timePoints(3,2),33);
str2_LT = sprintf('Target %d',Target_Mov2_LT);
plot(LT_Mov2,'Color', 'r'),title(str2_LT);
legend('Minimum Jerk','Linear');

%

MJ_Mov5 = VMJ1filt(timePoints(9,1):timePoints(10,1)); %4
Target_Mov5 = MJ1(timePoints(9,1),33);
str5 = sprintf('Target %d',Target_Mov5);
figure;plot(MJ_Mov5),title(str5);
hold on;

LT_Mov4 = VLT1filt(timePoints(7,2):timePoints(8,2));%4
Target_Mov4_LT = LT1(timePoints(7,2),33);
str4_LT = sprintf('Target %d',Target_Mov4_LT);
plot(LT_Mov4,'Color', 'r'),title(str4_LT);
legend('Minimum Jerk','Linear');

%

MJ_Mov6 = VMJ1filt(timePoints(11,1):timePoints(12,1)); %7
Target_Mov6 = MJ1(timePoints(11,1),33);
str6 = sprintf('Target %d',Target_Mov6);
figure;plot(MJ_Mov6),title(str6);
hold on;

LT_Mov6 = VLT1filt(timePoints(11,2):timePoints(12,2));%7
Target_Mov6_LT = LT1(timePoints(11,2),33);
str6_LT = sprintf('Target %d',Target_Mov6_LT);
plot(LT_Mov6,'Color', 'r'),title(str6_LT);
legend('Minimum Jerk','Linear');

%

MJ_Mov7 = VMJ1filt(timePoints(13,1):timePoints(14,1)); %1
Target_Mov7 = MJ1(timePoints(13,1),33);
str7 = sprintf('Target %d',Target_Mov7);
figure;plot(MJ_Mov7),title(str7);
hold on;

LT_Mov7 = VLT1filt(timePoints(13,2):timePoints(14,2));%1
Target_Mov7_LT = LT1(timePoints(13,2),33);
str7_LT = sprintf('Target %d',Target_Mov7_LT);
plot(LT_Mov7,'Color', 'r'),title(str7_LT);
legend('Minimum Jerk','Linear');

%

MJ_Mov8 = VMJ1filt(timePoints(15,1):timePoints(16,1)); %5
Target_Mov8 = MJ1(timePoints(15,1),33);
str8 = sprintf('Target %d',Target_Mov8);
figure;plot(MJ_Mov8),title(str8);
hold on;

LT_Mov1 = VLT1filt(timePoints(1,2):timePoints(2,2)); %5
Target_Mov1_LT = LT1(timePoints(1,2),33);
str1_LT = sprintf('Target %d',Target_Mov1_LT);
plot(LT_Mov1,'Color', 'r'),title(str1_LT);
legend('Minimum Jerk','Linear');




