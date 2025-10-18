

%% Fitting the sigmoid to the mean of the curve (D1)

clc
D21_mean_2 = D161_meanLc2c{3, 1};
coeffs_2 = D161_meanLc2c{6, 1}
gof_2 = D161_meanLc2c{7, 1}
D21_mean_4 = D161_meanLc4c{3, 1};
coeffs_4 = D161_meanLc4c{6, 1}
gof_4 = D161_meanLc4c{7, 1}

figure;
scatter(1:length(D21_mean_2),D21_mean_2, 'b');
hold on;
scatter(1:length(D21_mean_4),D21_mean_4, 'r');

x2=1:length(D21_mean_2);
x4=1:length(D21_mean_4);
plot(x2, coeffs_2.a./(1 + exp(-coeffs_2.b*(x2 - coeffs_2.c))) + coeffs_2.d, 'b', 'LineWidth', 3);
plot(x4, coeffs_4.a./(1 + exp(-coeffs_4.b*(x4 - coeffs_4.c))) + coeffs_4.d, 'r', 'LineWidth', 3);
xlabel('Trial Number');
ylabel('Learning Performance');
legend('Odor:2 No-Stim', 'Odor:4 Stim', 'Sigmoid Fit for Odor:2', 'Sigmoid Fit for Odor:4', 'Location', 'southeast');
title('Sigmoid Curve Fitting for Learning Curve');

%% Fitting the sigmoid to the mean of the curve (D2)

D2_mean_2 = D74_meanLc2c{3, 1};
D2coeffs_2 = D74_meanLc2c{6, 1};
D2gof_2 = D74_meanLc2c{7, 1};
D2_mean_4 = D74_meanLc4c{3, 1};
D2coeffs_4 = D74_meanLc4c{6, 1};
D2gof_4 = D74_meanLc4c{7, 1};

figure;
scatter(1:length(D2_mean_2), D2_mean_2, 'b');
hold on;
scatter(1:length(D2_mean_4), D2_mean_4, 'r');

x2 = 1:length(D2_mean_2);
x4 = 1:length(D2_mean_4);
plot(x2, D2coeffs_2.a./(1 + exp(-D2coeffs_2.b*(x2 - D2coeffs_2.c))) + D2coeffs_2.d, 'b', 'LineWidth', 3);
plot(x4, D2coeffs_4.a./(1 + exp(-D2coeffs_4.b*(x4 - D2coeffs_4.c))) + D2coeffs_4.d, 'r', 'LineWidth', 3);
xlabel('Trial Number');
ylabel('Learning Performance');
legend('Odor:2 No-Stim', 'Odor:4 Stim', 'Sigmoid Fit for Odor:2', 'Sigmoid Fit for Odor:4', 'Location', 'southeast');
title('Sigmoid Curve Fitting for Learning Curve');


