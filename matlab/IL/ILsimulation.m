clear 
clc 
global beta_ia temp_ia 
data1=readtable('ILEMA.csv');
timemark=312;
data1=table2array(data1(:,1));
I_il=data1.*2.74;%每周流感数据
data2=readtable('ILbeta.csv');
beta_ia=table2array(data2(:,2));
temp_ia=table2array(data2(:,3));
N=12830000;
E0=6;%7 4
I0=5;
R0=0;
S0=N-E0-I0-R0;

y0=[S0,E0,I0,R0,0];
tend=timemark*7;
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
sol=ode45(@seir_IL,[0 tend-1],y0,options);
x1=linspace(0,tend-1,timemark); %按周为一个时间单位
y1=deval(sol,x1);
newcase=y1(5,2:end)-y1(5,1:end-1);
week1 = 1:timemark-1;
week2 = 1:timemark-1;
startYear = 2013;
startWeek = 40;
startDate = datetime(startYear, 1, 1) + calweeks(startWeek - 1);
date1 = startDate + calweeks(week1 - 1);
date2 = startDate + calweeks(week2 - 1);
figure;

plot(date1,newcase,'LineWidth',1.5,'Color','b')
hold on
plot(date2,I_il(1:timemark-1),'r.', 'MarkerSize', 15)
xtickformat('yyyy')
xlabel('year');
ylabel('ILI cases');
legend('Predict curve','Real data')
% save('ILfit.mat','newcase')




