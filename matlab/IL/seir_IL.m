function dy=seir_IL(t,y)
global beta_ia temp_ia
temp=temp_ia;
if t<31
 beta=beta_ia(floor(t+1))/7;
else
dT=(temp(floor(t+1))-temp(floor(t+1-30)))/30;
T=temp(floor(t+1))*0.5629+0.1257*temp(floor(t))+0.1178*temp(floor(t-1))+0.0906*temp(floor(t-2))+...
    0.0512*temp(floor(t-3))+0.0258*temp(floor(t-4))+0.0147*temp(floor(t-5))+0.0113*temp(floor(t-6));
beta=1.0118*1.0064*exp(0.8088143390*(0.0000057972*T^3-0.0011048007*T^2+0.0050645368*T-0.4386011040*dT+0.1641361835))/7;
end
delta=0.0078/7;
sigma=1/2;
gamma=1/7;
dy=zeros(5,1);
S=y(1);
E=y(2);
I=y(3);
R=y(4);
N=12830000;
dy(1)=-beta*S*I/N+delta*R;
dy(2)=beta*S*I/N-sigma*E;
dy(3)=sigma*E-gamma*I;
dy(4)=gamma*I-delta*R;
dy(5)=sigma*E;
end