clear 
clc 
global beta_ia temp_ia 
data1=readtable('INEMA.csv');
timemark=312;
data1=table2array(data1(:,1));
I_il=data1.*2.74;%每周流感数据
data2=readtable('INbeta.csv');
beta_ia=table2array(data2(:,2));
temp_ia=table2array(data2(:,3));
% beta=[];
% temp=temp_ia;
% for t=30:2183
% dT=(temp(floor(t+1))-temp(floor(t+1-30)))/30;
% T=temp(floor(t+1))*0.5629+0.1257*temp(floor(t))+0.1178*temp(floor(t-1))+0.0906*temp(floor(t-2))+...
%     0.0512*temp(floor(t-3))+0.0258*temp(floor(t-4))+0.0147*temp(floor(t-5))+0.0113*temp(floor(t-6));
% % beta(t+1)=1.0125*exp(1.132181*(0.00001892*T^3-0.001225*T^2+0.00203298*T+0.219606*dT^2-0.40504*dT+0.124901));
% beta(t+1)=1.0125*exp(1.155591*((dT*(-0.4051327660725596)-(T *0.02586677997818288)^3)));
% end
% plot(beta(30:end))
% hold on
% plot(beta_ia(30:end),'+')
% sqrt(mean((beta_ia(31:end)-beta(31:end)').^2))
N=6630000;
E0=6;
I0=6;
R0=0;
S0=N-E0-I0-R0;

y0=[S0,E0,I0,R0,0];
tend=timemark*7;
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
sol=ode45(@seir_IN,[0 tend-1],y0,options);
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
save('INfit.mat','newcase')


