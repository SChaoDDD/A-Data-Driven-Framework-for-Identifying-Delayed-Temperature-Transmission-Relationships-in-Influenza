function dy=seir_OH(t,y)
global Tweight DT
T=Tweight(floor(t+1));
dT=DT(floor(t+1));
beta=0.9919*exp(0.9833*(0.0000057972*T^3-0.0011048007*T^2+0.0050645368*T-0.4386011040*dT+0.1641361835))/7;
delta=0.0078/7;
sigma=1/2;
gamma=1/7;
dy=zeros(5,1);
S=y(1);
E=y(2);
I=y(3);
R=y(4);
N=11632846;
dy(1)=-beta*S*I/N+delta*R;
dy(2)=beta*S*I/N-sigma*E;
dy(3)=sigma*E-gamma*I;
dy(4)=gamma*I-delta*R;
dy(5)=sigma*E;
end