clear;
clc;
close all;

%% =========================================================
% Semi-shared nonlinear temperature model (corrected)
%
% y_i = k_i * F(T_eff, DeltaT30)
%
% F =
% a*T_eff^3
% +b*T_eff^2
% +c*T_eff
% +d*DeltaT30^2
% +e*DeltaT30
% +f
%
% Shared:
% a,b,c,d,e,f
%
% Region specific:
% k_IA,k_IL,k_IN
%
% Target:
% y = ln(beta/beta_bar)
%% =========================================================


%% 1. Temperature weights

weights = [
    0.5629;
    0.1257;
    0.1178;
    0.0906;
    0.0512;
    0.0258;
    0.0147;
    0.0113
];


files = {
    'IAbeta.csv'
    'ILbeta.csv'
    'INbeta.csv'
};


region_names = {
    'IA'
    'IL'
    'IN'
};



%% 2. Read and preprocess

IA = preprocess_temperature(readtable(files{1}),weights);
IL = preprocess_temperature(readtable(files{2}),weights);
IN = preprocess_temperature(readtable(files{3}),weights);

regions = {IA,IL,IN};



%% 3. Construct data

X = [];
y = [];
region_id = [];


for i = 1:3

    data = regions{i};

    beta = data.beta;

    beta_bar = mean(beta);

    yi = beta./beta_bar;


    X = [
        X;
        data.temp_eff, data.deltaT30
    ];


    y = [
        y;
        yi
    ];


    region_id = [
        region_id;
        i*ones(length(yi),1)
    ];

end



%% 4. Semi-shared model

model = @(p,X) semi_model(p,X,region_id);



%% 5. Initial parameters

% a b c d e f kIA kIL kIN

p0 = [
    zeros(6,1);
    1;
    1;
    1
];



options = optimoptions( ...
    'lsqcurvefit', ...
    'Display','iter', ...
    'MaxIterations',3000 ...
);



%% 6. Parameter estimation

p_est = lsqcurvefit( ...
    model,...
    p0,...
    X,...
    y,...
    [],...
    [],...
    options ...
);



%% 7. Output coefficients

fprintf('\n==============================\n');
fprintf('Shared coefficients\n');
fprintf('==============================\n');

fprintf('a = %.10f\n',p_est(1));
fprintf('b = %.10f\n',p_est(2));
fprintf('c = %.10f\n',p_est(3));
fprintf('d = %.10f\n',p_est(4));
fprintf('e = %.10f\n',p_est(5));
fprintf('f = %.10f\n',p_est(6));


fprintf('\n==============================\n');
fprintf('Region coefficients\n');
fprintf('==============================\n');

fprintf('k_IA = %.10f\n',p_est(7));
fprintf('k_IL = %.10f\n',p_est(8));
fprintf('k_IN = %.10f\n',p_est(9));



%% 8. Overall evaluation

ypred_all = semi_model(p_est,X,region_id);


RMSE = sqrt(mean((y-ypred_all).^2));

R2 = 1 - sum((y-ypred_all).^2) / ...
    sum((y-mean(y)).^2);


fprintf('\nOverall RMSE = %.6f\n',RMSE);
fprintf('Overall R2 = %.6f\n',R2);



%% 9. Region prediction plots

for i = 1:3

    data = regions{i};

    idx = region_id == i;


    y_pred = semi_model(p_est,X(idx,:),region_id(idx));


    beta_bar = mean(data.beta);


    beta_pred = beta_bar .* y_pred;



    rmse_i = sqrt(mean((data.beta-beta_pred).^2));


    fprintf('\n%s beta RMSE = %.6f\n', ...
        region_names{i}, rmse_i);



    figure;

    plot(...
        data.Day,...
        data.beta,...
        'LineWidth',2);

    hold on;


    plot(data.Day,...
        beta_pred,...
        'LineWidth',2);


    xlabel('Day');
    ylabel('beta');

    title([
        region_names{i}, ...
        ' Semi-shared nonlinear model'
    ]);

% save('INbetafit.mat','beta_pred')


    grid on;



    
end



%% =========================================================
% Functions
%% =========================================================


function yhat = semi_model(p,X,region_id)


F = p(1).*X(:,1).^3+...
    p(2).*X(:,1).^2 + ...
    p(3).*X(:,1) + ...
    p(5).*X(:,2) + ...
    p(6);



k = zeros(length(region_id),1);


for i = 1:length(region_id)

    if region_id(i)==1

        k(i)=p(7);

    elseif region_id(i)==2

        k(i)=p(8);

    elseif region_id(i)==3

        k(i)=p(9);

    end

end


yhat = exp(k .* F);


end




function df = preprocess_temperature(df,weights)


temp = df.temp;

n = length(temp);


lag = length(weights)-1;


Teff = zeros(n-lag,1);


for t = lag+1:n

    value = 0;


    for k = 0:lag

        value = value + weights(k+1)*temp(t-k);

    end


    Teff(t-lag)=value;

end



deltaT30 = zeros(n-30,1);


for t = 31:n

    deltaT30(t-30) = ...
        (temp(t)-temp(t-30))/30;

end



start_day = 31;


df = df(start_day:end,:);


df.temp_eff = Teff(start_day-lag:end);


df.deltaT30 = deltaT30(start_day-30:end);


end
