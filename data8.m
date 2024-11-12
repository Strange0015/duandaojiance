clc; 
close all; 
clear all;
file_0=1;%第一个读取文件
file_1=2;%第二个读取文件    
file_end=5000;%最后一个读取文件
%%采集文件读取
folder = 'D:\研究生资料\断刀检测\数据\2024_11_7\z'; % 替换为你的文件夹路径
files = dir(fullfile(folder, '*.csv'));
files = files(~[files.isdir]);% 排除子文件夹
%%采集文件按文件名排序
sorted_files = sort_nat({files.name});
file0 = fullfile(folder, sorted_files{1,file_0});
%提取第一个.csv文件的铣削电流数据存入y
m0=csvread(file0,9,0);
% %分号分隔符
% m1 = importdata(file0, ';');
% m0=m1.data(1:50000,:);
y=m0(:,2);
%%开始遍历第二个文件到最后一个文件，前面已经提取第一个文件数据完毕
for i = file_1:file_end
    % length(sorted_files)遍历并处理CSV文件
    file = fullfile(folder, sorted_files{1,i}); % file = fullfile(folder, sorted_files(i).name); 
    % 在这里编写你的代码来处理CSV文件
    m=csvread(file,9,0);
%     m2 = importdata(file, ';');
%     m=m2.data(1:50000,:);
    y2=m(:,2);
    y=[y;y2];
    disp(['处理文件名称: ', file]);
end
%%设置采样频率5000Hz绘制时间电流图像/
Fs=5000;%采样频率
Ts=1/Fs;%采样时间0.001s
N=length(y);%实际电流个数
N1=N-1;
t=0:Ts:N1*Ts;
figure(1);
plot(t,y);
xlabel('X-采样时间');
ylabel('Y-实际电流');
title('原始数据');

%%中值滤波，绘制滤波后图像
y1=smooth(y,100);
figure(2);
plot(t,y1);
xlabel('X-采样时间/s');
ylabel('Y-滤波后实际电流/A');
title('滤波后图像');
% ylim([2.46 2.61]);
grid on;
%%铣削电流上升沿检测 
k=1;
p=1;
while k<N-5000
    data_1=y1(k:k+5000);%将电流数据每5000个为一组存入中间变量data_1
    rms_value(p,1)=rms(data_1);%分组的电流数据进行均方根计算，并存入rms数组
    k=k+5000;%下一个电流数据组初始位置
    p=p+1;%每个电流数据组的均方根数组位置
end
t_1=0:1:length(rms_value)-1;%分组后数据时间划分，时间间隔1s
figure(3);
plot(t_1,rms_value);%输出：时间-各分组均方根函数图像（各分组时间为每秒）
xlabel('X-时间/秒');
ylabel('Y-各分组均方根/A');
ylim([1.3 2.5]);
title('时间-各分组均方根函数图像（调试上升沿用）');
title('3号刀铣削图像');
B=diff(rms_value);%各个均方根组的差分值，即前值与后值做差，即每秒的斜率
upward_data=rms_value(B>=0.017);%筛选出差分值大于0.017的上升沿
% up_1=upward_data(upward_data>0.9);
time=find(B(2:end)>=0.017)+3;%得出上升沿后时间两秒的时间,单位秒
%去除相近上升沿数据，容差为9
tolerance=12;
time1=remove_near_duplicates(time,tolerance);%去除相近上升沿后时间
% max=10*lenth(file_end);
% 逻辑索引，找出矩阵中不大于阈值的元素
time1(time1>98330)=[];
k=length(time1);
%崩刃检测
    f_0=5000*time1(k-1);%倒数第二个上升沿时间
    f_1=5000*time1(k);%最后一个上升沿时间
    subArray_0=t(f_0:f_0+5*5000);%抓取时间t中的数据
    columnElements_0=y1(f_0:f_0+5*5000,1);%抓取y中的数据
    columnElements_0_1=y1(f_0+10*5000:f_0+12*5000,1);%抓取y中的数据
    up_rms_0=rms(columnElements_0);%倒数第二个上升沿的RMS值
    up_rms_1=rms(columnElements_0_1);%倒数第二个下降沿的RMS值
    subArray_1=t(f_1:f_1+20*5000);%抓取时间t中的数据
    columnElements_1=y1(f_1:f_1+20*5000,1);%抓取y中的数据
    broken=[subArray_1' columnElements_1];%最后一个上升沿时间和电流的矩阵，第一列为时间，第二列为对应电流值
    threshold =0.5*(up_rms_0+up_rms_1);%下降沿判断标准
    selected_values = broken(broken(:, 2) <threshold, 1);%筛选出小于判断标准的电流
    fprintf('刀具在%f秒停轴或崩刃\n', selected_values(1));%输出一个符合标准电流对应的时间即为崩刃时间
m=1;
while m<k+1
f=5000*time1(m);
    subArray=t(f:f+5*5000);%抓取时间t中的数据
    columnElements=y1(f:f+5*5000,1);%抓取y中的数据
    up_rms_0=rms(columnElements);
    %%如果上升沿电流数据高于2.45则说明铣削电流数据有效记录数据，否则记录的上一个铣削电流数据作为此槽数数据
    if up_rms_0>2.45    
    up_rms(m,1)=rms(columnElements);
    %切削电流与空转电流增量（空转电流取切削电流前5-3秒的数据做RMS值）
    subArray1=t(f-6*5000:f-4*5000);%抓取时间t中的数据
    columnElements1=y1(f-6*5000:f-4*5000,1);%抓取y中的数据
    up_rms4(m,1)=rms(columnElements1);
    increment(m,1)=up_rms(m,1)-up_rms4(m,1);
    m=m+1;
    else  
    up_rms(m,1)=up_rms(m-1,1);
    %切削电流与空转电流增量（空转电流取切削电流前5-3秒的数据做RMS值）
    f=5000*time1(m-1);
    subArray1=t(f-6*5000:f-4*5000);%抓取时间t中的数据
    columnElements1=y1(f-6*5000:f-4*5000,1);%抓取y中的数据
    up_rms4(m,1)=rms(columnElements1);
    increment(m,1)=up_rms(m,1)-up_rms4(m,1);
    m=m+1;
    end
end
data_3(:,1)=up_rms;
date_end=[time1,up_rms];
t1=1:1:length(time1); 
figure(4);
plot(t1,up_rms);
xlabel('X-槽数');
ylabel('Y-滤波后上升沿RMS值');
title('结果图像');
% ylim([2.46 2.61]);
grid on;
% 计算移动平均值
up_rms1 = smooth(up_rms, 100);
% 计算偏差
deviations = abs(up_rms - up_rms1);
% 过滤掉偏离阈值0.03以外的数据点停轴按0.05
up_rms2 = up_rms(deviations < 0.05);
% 去除平滑后数据末端的点
up_rms3 = up_rms2(1:end-1); 
% up_rms2=up_rms(inRange);
t2=1:1:length(up_rms3); 
figure(5);
plot(t2,up_rms3);
xlabel('X-槽数');
ylabel('Y-去除异常后上升沿RMS值');
title('结果图像');
% ylim([2.45 2.65]);
% ylim([-0.2 1]);
grid on;

%increment(m,1)=up_rms(m,1)-up_rms4(m,1);
increment4=[up_rms,up_rms4];
% 计算移动平均值
increment(increment>0.12)=[];
increment(increment<0.03)=[];
increment1 = smooth(increment, 100);
% 计算偏差
deviations = abs(increment - increment1);
% 过滤掉偏离阈值0.03以外的数据点停轴按0.5
increment2 = increment(deviations < 0.04);
% 去除平滑后数据末端的点
increment3 = increment2(1:end-1); 
% up_rms2=up_rms(inRange);
t3=1:1:length(increment); 
t4=1:1:length(increment3); 
figure(7);
plot(t3,increment);
xlabel('X-槽数');
ylabel('Y-增量RMS值');
title('增量结果图像');
figure(8);
plot(t4,increment3);
xlabel('X-槽数');
ylabel('Y-去除异常后增量RMS值');
title('增量结果图像');
% ylim([-0.03 0.08]);
grid on;

%增量阈值判断
k=1;
p=1;
N=length(increment3);
while k<N-10
    data_3=increment3(k:k+9);%将电流数据每1000个为一组存入中间变量data_1
    data_4(p,1)=rms(data_3);%分组的电流数据进行均方根计算，并存入rms数组
    k=k+10;%下一个电流数据组初始位置
    p=p+1;%每个电流数据组的均方根数组位置
end
t_2=1:1:length(data_4);%分组后数据时间划分，范围0-86s，时间间隔1s
figure(9);
plot(t_2,data_4);%输出：时间-各分组均方根函数图像 
xlabel('X-槽数组/s');
ylabel('Y-各分组均方根/A');
title('时间-各分组均方根函数图像（检测下降沿）');
F2=[t_2',data_4];
column_indices_0 = F2(:, 2) > 0.075; % 生成逻辑索引0.075
E3 = F2(column_indices_0, 1);
inRange2 =(E3 >= 0.075); 
result3 = E3(inRange2); % 符合条件的元素
if ~all(result3(:) == 0)
    fprintf('数据增量达到0.075在：%d0槽，刀具存在崩刃可能，请检查刀具状态\n', result3(1));
end
column_indices_1 = F2(:, 2) > 0.08; % 生成逻辑索引0.075
E4 = F2(column_indices_1, 1);
inRange3 = (E4 >= 0.08); 
result4 = E4(inRange3); % 符合条件的元素
if ~all(result4(:) == 0)
    fprintf('数据增量达到0.08在：%d0槽，刀具存在崩刃可能，请检查刀具状态\n', result4(1));
end
column_indices_2 = F2(:, 2) > 0.085; % 生成逻辑索引0.075
E5 = F2(column_indices_2, 1);
inRange4 = (E5 >= 0.085); 
result5 = E5(inRange4); % 符合条件的元素
if ~all(result5(:) == 0)
    fprintf('数据增量达到0.085在：%d0槽，刀具存在崩刃可能，请检查刀具状态\n', result5(1));
end

%2000槽以内缓慢上升监测
k=1;
p=1;
N=length(up_rms3);
while k<N-100
    data_1=up_rms3(k:k+99);%将电流数据每1000个为一组存入中间变量data_1
    data_2(p,1)=rms(data_1);%分组的电流数据进行均方根计算，并存入rms数组
    k=k+100;%下一个电流数据组初始位置
    p=p+1;%每个电流数据组的均方根数组位置
end
t_1=1:1:length(data_2);%分组后数据时间划分，范围0-86s，时间间隔1s
figure(6);
plot(t_1,data_2);%输出：时间-各分组均方根函数图像 
xlabel('X-槽数组/s');
ylabel('Y-各分组均方根/A');
title('时间-各分组均方根函数图像（检测下降沿）');
D=data_2;%筛选出差分值小于-0.012的下降沿
a1=2.45;%a1=mean(up_rms(1:300));%第一段均方根平均值 
D1=D-a1;%与前300组数据增量
F=[t_1',D1];

% 检查第二列是否满足0.06增量条件，并使用逻辑索引获取同行第一列的值
column_indices = F(:, 2) > 0.10; % 生成逻辑索引0.14
E2 = F(column_indices, 1); % 提取对应的第一列值
% 逻辑索引，找出矩阵中大于阈值的元素
inRange1 = (E2 >=10) & (E2 <= 20); % 逻辑索引
result1 = E2(inRange1); % 符合条件的元素
% 检查列向量是否不是零向量
if ~all(result1(:) == 0)
    fprintf('上升趋势数据增量达到0.10在：%d00槽，刀具存在崩刃的可能，请检查刀具状态\n', result1(1));  % 输出第一个数据
end
%2500槽内达到增量的下降趋势监测
C=diff(data_2);%各个均方根组的差分值，即前值与后值做差
D2=data_2(C<=-0.002);%筛选出差分值小于-0.012的下降沿
E=find(C(1:end)<=-0.002)+1;%得出下降沿槽数,单位秒
D3=D2-a1;%与前300组数据增量
F1=[E,D3];
% 检查第二列是否满足0.06增量条件，并使用逻辑索引获取同行第一列的值
column_indices = F1(:, 2) > 0.075; % 生成逻辑索引
E3 = F1(column_indices, 1); % 提取对应的第一列值
% 逻辑索引，找出矩阵中不大于阈值的元素
inRange1 = (E3 >= 17) & (E3 <= 23); % 逻辑索引17-23
result2 = E3(inRange1); % 符合条件的元素
if ~all(result2(:) == 0)
    fprintf('下降趋势数据增量达到0.075在：%d00槽，刀具存在崩刃可能，请检查刀具状态\n', result2(1));
end
% disp(result1); % 显示结果 


