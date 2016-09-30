data = csvread('../prob/2i_data.csv');

plot(data(:, 1), data(:, 2), 'k-', data(:, 1), data(:, 3), 'k--', data(:, 1), data(:, 4), 'k:');
legend('U','K','H=K+U');
xlabel('time (sec)');
xlim([0.0, 20.0]);
ylim([0.0, 1.1]);