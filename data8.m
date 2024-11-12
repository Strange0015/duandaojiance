clc; 
close all; 
clear all;
file_0=1;%��һ����ȡ�ļ�
file_1=2;%�ڶ�����ȡ�ļ�    
file_end=5000;%���һ����ȡ�ļ�
%%�ɼ��ļ���ȡ
folder = 'D:\�о�������\�ϵ����\����\2024_11_7\z'; % �滻Ϊ����ļ���·��
files = dir(fullfile(folder, '*.csv'));
files = files(~[files.isdir]);% �ų����ļ���
%%�ɼ��ļ����ļ�������
sorted_files = sort_nat({files.name});
file0 = fullfile(folder, sorted_files{1,file_0});
%��ȡ��һ��.csv�ļ���ϳ���������ݴ���y
m0=csvread(file0,9,0);
% %�ֺŷָ���
% m1 = importdata(file0, ';');
% m0=m1.data(1:50000,:);
y=m0(:,2);
%%��ʼ�����ڶ����ļ������һ���ļ���ǰ���Ѿ���ȡ��һ���ļ��������
for i = file_1:file_end
    % length(sorted_files)����������CSV�ļ�
    file = fullfile(folder, sorted_files{1,i}); % file = fullfile(folder, sorted_files(i).name); 
    % �������д��Ĵ���������CSV�ļ�
    m=csvread(file,9,0);
%     m2 = importdata(file, ';');
%     m=m2.data(1:50000,:);
    y2=m(:,2);
    y=[y;y2];
    disp(['�����ļ�����: ', file]);
end
%%���ò���Ƶ��5000Hz����ʱ�����ͼ��/
Fs=5000;%����Ƶ��
Ts=1/Fs;%����ʱ��0.001s
N=length(y);%ʵ�ʵ�������
N1=N-1;
t=0:Ts:N1*Ts;
figure(1);
plot(t,y);
xlabel('X-����ʱ��');
ylabel('Y-ʵ�ʵ���');
title('ԭʼ����');

%%��ֵ�˲��������˲���ͼ��
y1=smooth(y,100);
figure(2);
plot(t,y1);
xlabel('X-����ʱ��/s');
ylabel('Y-�˲���ʵ�ʵ���/A');
title('�˲���ͼ��');
% ylim([2.46 2.61]);
grid on;
%%ϳ�����������ؼ�� 
k=1;
p=1;
while k<N-5000
    data_1=y1(k:k+5000);%����������ÿ5000��Ϊһ������м����data_1
    rms_value(p,1)=rms(data_1);%����ĵ������ݽ��о��������㣬������rms����
    k=k+5000;%��һ�������������ʼλ��
    p=p+1;%ÿ������������ľ���������λ��
end
t_1=0:1:length(rms_value)-1;%���������ʱ�仮�֣�ʱ����1s
figure(3);
plot(t_1,rms_value);%�����ʱ��-���������������ͼ�񣨸�����ʱ��Ϊÿ�룩
xlabel('X-ʱ��/��');
ylabel('Y-�����������/A');
ylim([1.3 2.5]);
title('ʱ��-���������������ͼ�񣨵����������ã�');
title('3�ŵ�ϳ��ͼ��');
B=diff(rms_value);%������������Ĳ��ֵ����ǰֵ���ֵ�����ÿ���б��
upward_data=rms_value(B>=0.017);%ɸѡ�����ֵ����0.017��������
% up_1=upward_data(upward_data>0.9);
time=find(B(2:end)>=0.017)+3;%�ó������غ�ʱ�������ʱ��,��λ��
%ȥ��������������ݣ��ݲ�Ϊ9
tolerance=12;
time1=remove_near_duplicates(time,tolerance);%ȥ����������غ�ʱ��
% max=10*lenth(file_end);
% �߼��������ҳ������в�������ֵ��Ԫ��
time1(time1>98330)=[];
k=length(time1);
%���м��
    f_0=5000*time1(k-1);%�����ڶ���������ʱ��
    f_1=5000*time1(k);%���һ��������ʱ��
    subArray_0=t(f_0:f_0+5*5000);%ץȡʱ��t�е�����
    columnElements_0=y1(f_0:f_0+5*5000,1);%ץȡy�е�����
    columnElements_0_1=y1(f_0+10*5000:f_0+12*5000,1);%ץȡy�е�����
    up_rms_0=rms(columnElements_0);%�����ڶ��������ص�RMSֵ
    up_rms_1=rms(columnElements_0_1);%�����ڶ����½��ص�RMSֵ
    subArray_1=t(f_1:f_1+20*5000);%ץȡʱ��t�е�����
    columnElements_1=y1(f_1:f_1+20*5000,1);%ץȡy�е�����
    broken=[subArray_1' columnElements_1];%���һ��������ʱ��͵����ľ��󣬵�һ��Ϊʱ�䣬�ڶ���Ϊ��Ӧ����ֵ
    threshold =0.5*(up_rms_0+up_rms_1);%�½����жϱ�׼
    selected_values = broken(broken(:, 2) <threshold, 1);%ɸѡ��С���жϱ�׼�ĵ���
    fprintf('������%f��ͣ������\n', selected_values(1));%���һ�����ϱ�׼������Ӧ��ʱ�伴Ϊ����ʱ��
m=1;
while m<k+1
f=5000*time1(m);
    subArray=t(f:f+5*5000);%ץȡʱ��t�е�����
    columnElements=y1(f:f+5*5000,1);%ץȡy�е�����
    up_rms_0=rms(columnElements);
    %%��������ص������ݸ���2.45��˵��ϳ������������Ч��¼���ݣ������¼����һ��ϳ������������Ϊ�˲�������
    if up_rms_0>2.45    
    up_rms(m,1)=rms(columnElements);
    %�����������ת������������ת����ȡ��������ǰ5-3���������RMSֵ��
    subArray1=t(f-6*5000:f-4*5000);%ץȡʱ��t�е�����
    columnElements1=y1(f-6*5000:f-4*5000,1);%ץȡy�е�����
    up_rms4(m,1)=rms(columnElements1);
    increment(m,1)=up_rms(m,1)-up_rms4(m,1);
    m=m+1;
    else  
    up_rms(m,1)=up_rms(m-1,1);
    %�����������ת������������ת����ȡ��������ǰ5-3���������RMSֵ��
    f=5000*time1(m-1);
    subArray1=t(f-6*5000:f-4*5000);%ץȡʱ��t�е�����
    columnElements1=y1(f-6*5000:f-4*5000,1);%ץȡy�е�����
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
xlabel('X-����');
ylabel('Y-�˲���������RMSֵ');
title('���ͼ��');
% ylim([2.46 2.61]);
grid on;
% �����ƶ�ƽ��ֵ
up_rms1 = smooth(up_rms, 100);
% ����ƫ��
deviations = abs(up_rms - up_rms1);
% ���˵�ƫ����ֵ0.03��������ݵ�ͣ�ᰴ0.05
up_rms2 = up_rms(deviations < 0.05);
% ȥ��ƽ��������ĩ�˵ĵ�
up_rms3 = up_rms2(1:end-1); 
% up_rms2=up_rms(inRange);
t2=1:1:length(up_rms3); 
figure(5);
plot(t2,up_rms3);
xlabel('X-����');
ylabel('Y-ȥ���쳣��������RMSֵ');
title('���ͼ��');
% ylim([2.45 2.65]);
% ylim([-0.2 1]);
grid on;

%increment(m,1)=up_rms(m,1)-up_rms4(m,1);
increment4=[up_rms,up_rms4];
% �����ƶ�ƽ��ֵ
increment(increment>0.12)=[];
increment(increment<0.03)=[];
increment1 = smooth(increment, 100);
% ����ƫ��
deviations = abs(increment - increment1);
% ���˵�ƫ����ֵ0.03��������ݵ�ͣ�ᰴ0.5
increment2 = increment(deviations < 0.04);
% ȥ��ƽ��������ĩ�˵ĵ�
increment3 = increment2(1:end-1); 
% up_rms2=up_rms(inRange);
t3=1:1:length(increment); 
t4=1:1:length(increment3); 
figure(7);
plot(t3,increment);
xlabel('X-����');
ylabel('Y-����RMSֵ');
title('�������ͼ��');
figure(8);
plot(t4,increment3);
xlabel('X-����');
ylabel('Y-ȥ���쳣������RMSֵ');
title('�������ͼ��');
% ylim([-0.03 0.08]);
grid on;

%������ֵ�ж�
k=1;
p=1;
N=length(increment3);
while k<N-10
    data_3=increment3(k:k+9);%����������ÿ1000��Ϊһ������м����data_1
    data_4(p,1)=rms(data_3);%����ĵ������ݽ��о��������㣬������rms����
    k=k+10;%��һ�������������ʼλ��
    p=p+1;%ÿ������������ľ���������λ��
end
t_2=1:1:length(data_4);%���������ʱ�仮�֣���Χ0-86s��ʱ����1s
figure(9);
plot(t_2,data_4);%�����ʱ��-���������������ͼ�� 
xlabel('X-������/s');
ylabel('Y-�����������/A');
title('ʱ��-���������������ͼ�񣨼���½��أ�');
F2=[t_2',data_4];
column_indices_0 = F2(:, 2) > 0.075; % �����߼�����0.075
E3 = F2(column_indices_0, 1);
inRange2 =(E3 >= 0.075); 
result3 = E3(inRange2); % ����������Ԫ��
if ~all(result3(:) == 0)
    fprintf('���������ﵽ0.075�ڣ�%d0�ۣ����ߴ��ڱ��п��ܣ����鵶��״̬\n', result3(1));
end
column_indices_1 = F2(:, 2) > 0.08; % �����߼�����0.075
E4 = F2(column_indices_1, 1);
inRange3 = (E4 >= 0.08); 
result4 = E4(inRange3); % ����������Ԫ��
if ~all(result4(:) == 0)
    fprintf('���������ﵽ0.08�ڣ�%d0�ۣ����ߴ��ڱ��п��ܣ����鵶��״̬\n', result4(1));
end
column_indices_2 = F2(:, 2) > 0.085; % �����߼�����0.075
E5 = F2(column_indices_2, 1);
inRange4 = (E5 >= 0.085); 
result5 = E5(inRange4); % ����������Ԫ��
if ~all(result5(:) == 0)
    fprintf('���������ﵽ0.085�ڣ�%d0�ۣ����ߴ��ڱ��п��ܣ����鵶��״̬\n', result5(1));
end

%2000�����ڻ����������
k=1;
p=1;
N=length(up_rms3);
while k<N-100
    data_1=up_rms3(k:k+99);%����������ÿ1000��Ϊһ������м����data_1
    data_2(p,1)=rms(data_1);%����ĵ������ݽ��о��������㣬������rms����
    k=k+100;%��һ�������������ʼλ��
    p=p+1;%ÿ������������ľ���������λ��
end
t_1=1:1:length(data_2);%���������ʱ�仮�֣���Χ0-86s��ʱ����1s
figure(6);
plot(t_1,data_2);%�����ʱ��-���������������ͼ�� 
xlabel('X-������/s');
ylabel('Y-�����������/A');
title('ʱ��-���������������ͼ�񣨼���½��أ�');
D=data_2;%ɸѡ�����ֵС��-0.012���½���
a1=2.45;%a1=mean(up_rms(1:300));%��һ�ξ�����ƽ��ֵ 
D1=D-a1;%��ǰ300����������
F=[t_1',D1];

% ���ڶ����Ƿ�����0.06������������ʹ���߼�������ȡͬ�е�һ�е�ֵ
column_indices = F(:, 2) > 0.10; % �����߼�����0.14
E2 = F(column_indices, 1); % ��ȡ��Ӧ�ĵ�һ��ֵ
% �߼��������ҳ������д�����ֵ��Ԫ��
inRange1 = (E2 >=10) & (E2 <= 20); % �߼�����
result1 = E2(inRange1); % ����������Ԫ��
% ����������Ƿ���������
if ~all(result1(:) == 0)
    fprintf('�����������������ﵽ0.10�ڣ�%d00�ۣ����ߴ��ڱ��еĿ��ܣ����鵶��״̬\n', result1(1));  % �����һ������
end
%2500���ڴﵽ�������½����Ƽ��
C=diff(data_2);%������������Ĳ��ֵ����ǰֵ���ֵ����
D2=data_2(C<=-0.002);%ɸѡ�����ֵС��-0.012���½���
E=find(C(1:end)<=-0.002)+1;%�ó��½��ز���,��λ��
D3=D2-a1;%��ǰ300����������
F1=[E,D3];
% ���ڶ����Ƿ�����0.06������������ʹ���߼�������ȡͬ�е�һ�е�ֵ
column_indices = F1(:, 2) > 0.075; % �����߼�����
E3 = F1(column_indices, 1); % ��ȡ��Ӧ�ĵ�һ��ֵ
% �߼��������ҳ������в�������ֵ��Ԫ��
inRange1 = (E3 >= 17) & (E3 <= 23); % �߼�����17-23
result2 = E3(inRange1); % ����������Ԫ��
if ~all(result2(:) == 0)
    fprintf('�½��������������ﵽ0.075�ڣ�%d00�ۣ����ߴ��ڱ��п��ܣ����鵶��״̬\n', result2(1));
end
% disp(result1); % ��ʾ��� 


