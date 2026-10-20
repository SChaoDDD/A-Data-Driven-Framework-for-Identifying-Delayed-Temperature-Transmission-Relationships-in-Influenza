d1=readtable('INdaily_beta_ensemble1.csv');
d2=readtable('allbeta.csv');
beta1=table2array(d1(:,3));
beta2=table2array(d2(:,5));
plot(beta1,'b')
hold on
plot(beta2)