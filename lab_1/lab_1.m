% Define constant parameters
clear; close all; clc;
M1 = 150;
% M2 = 60;
k = 100;
F = 250;
mu = 0.05;
g = 9.81;

M2_values = 100:100:2000;

% Store simulation results
x1 = {};
v1 = {}; 
x2 = {}; 
v2 = {};
legends = {}; % To store legend entries
t = {};

% Loop over each value of M2
for i = 1:length(M2_values)
    M2 = M2_values(i); % Set the current value of M2
    simOut = sim('lab1.slx');
    x1{i} = simOut.x1;
    v1{i} = simOut.v1;
    x2{i} = simOut.x2;
    v2{i} = simOut.v2;
    t{i} = simOut.tout;

    legends{i} = ['M2 = ' num2str(M2)];
end

% Plot results
figure;
hold on;
for i = 1:length(x1)
    plot(t{i}, x1{i}); 
end
hold off;

% Add legend and labels
legend(legends, Location="southeast");
xlabel('Time (s)');
ylabel('Robot displacement (m)');
title('Robot displacement with different value of loading mass');
saveas(gcf,'pos_1_multiple.png', "png")

figure;
hold on;
for i = 1:length(x1)
    plot(t{i}, v1{i}); 
end
hold off;

% Add legend and labels
legend(legends, Location="southeast");
xlabel('Time (s)');
ylabel('Robot velocity (m/s)');
title('Robot velocity with different value of loading mass');
saveas(gcf,'vel_1_multiple.png', "png")

figure;
hold on;
for i = 1:length(x1)
    plot(t{i}, x2{i}); 
end
hold off;

% Add legend and labels
legend(legends, Location="southeast");
xlabel('Time (s)');
ylabel('Load displacement (m)');
title('Load displacement with different value of loading mass');
saveas(gcf,'pos_2_multiple.png', "png")

figure;
hold on;
for i = 1:length(x1)
    plot(t{i}, v2{i}); 
end
hold off;

% Add legend and labels
legend(legends, Location="southeast");
xlabel('Time (s)');
ylabel('Load velocity (m/s)');
title('Load velocity with different value of loading mass');
saveas(gcf,'vel_2_multiple.png', "png")

figure;

subplot(221)
plot(t{6}, x1{6}); 
xlabel('Time (s)');
ylabel('Robot displacement (m)');
title('Simulation of robot displacement');

subplot(222)
plot(t{6}, v1{6});
xlabel('Time (s)');
ylabel('Robot velocity (m/s)');
title('Simulation of robot velocity');

subplot(223)
plot(t{6}, x2{6});
xlabel('Time (s)');
ylabel('Load displacement (m)');
title('Simulation of load displacement');

subplot(224)
plot(t{6},v2{6});
xlabel('Time (s)');
ylabel('Load velocity (m/s)');
title('Simulation of load velocity');
saveas(gcf,'robot_sim.png', "png")

mass_index = 19;
% Calculate final value and tolerance
v1_end = v1{mass_index}(end);
tolerance = 0.01 * abs(v1_end);

% Moving average
n = 1000; %n represent how many rows in the final result
v1_averaged = reshape(v1{mass_index}(1:10000), n, []);
v1_averaged = abs(v1_averaged - v1_end);
v1_averaged = mean(v1_averaged, 2);
% Find the settling time index
t_settle_index = find(v1_averaged <= tolerance, 1, 'first');

% Plot the original data
figure;
plot(t{mass_index}, v1{mass_index});
xlabel('Time (s)');
ylabel('Robot velocity (m/s)');
title('Robot velocity settling time');
hold on;

% Find the time at which settling occurs
t_settle_index = round(10000/n*(1/2+t_settle_index));
t_settle = t{mass_index}(t_settle_index);

% Zoom in plot
zoom_window = 300;
zoom_indices = max(t_settle_index - zoom_window, 1) : min(t_settle_index + zoom_window, length(t{mass_index}));

yline(v1_end, '--r', '1% Tolerance');
xline(t_settle, '--k', 'Settling Time');

inset_pos = [0.55 0.55 0.3 0.3]; 
inset_axes = axes('Position', inset_pos);

plot(inset_axes, t{mass_index}(zoom_indices), v1{mass_index}(zoom_indices));
xlabel(inset_axes, 'Time (s)');
ylabel(inset_axes, 'Load velocity (m/s)');
title(inset_axes, 'Zoomed View');
xlim(inset_axes, [t{mass_index}(zoom_indices(1)), t{mass_index}(zoom_indices(end))]);
% ylim(inset_axes, [2.2, 2.5]);

hold(inset_axes, 'on');
yline(inset_axes, v1_end, '--r');
xline(inset_axes, t_settle, '--k');
hold(inset_axes, 'off');
%%
clear; close all; clc;
M1 = 150;
% M2 = 60;
k = 100;
F = 250;
mu = 0.05;
g = 9.81;

M2_values = 100:100:2000;

% Store simulation results
x1 = {};
v1 = {}; 
x2 = {}; 
v2 = {};
legends = {}; % To store legend entries
t = {};

% Loop over each value of M2
for i = 1:length(M2_values)
    M2 = M2_values(i); % Set the current value of M2
    simOut = sim('lab1.slx');
    x1{i} = simOut.x1;
    v1{i} = simOut.v1;
    x2{i} = simOut.x2;
    v2{i} = simOut.v2;
    t{i} = simOut.tout;
end

t_settle_array = [];
for i = 1:length(v1)

% Find the time at which settling occurs

v1_end = v1{i}(end);
tolerance = 0.01 * abs(v1_end);

% Moving average
n = 500; %n represent how many rows in the final result
v1_averaged = reshape(v1{i}(1:10000), n, []);
v1_averaged = abs(v1_averaged - v1_end);
v1_averaged = mean(v1_averaged, 2);
% Find the settling time index
t_settle_index = find(v1_averaged <= tolerance, 1, 'first');

% Find the settling time index
t_settle_index = round(10000/n*(1/2+t_settle_index));
t_settle_index = find(abs(v1{i} - v1_end) <= tolerance, 1, 'first');
fprintf("%d %.2f\n", M2_values(i), t{i}(t_settle_index));
t_settle_array(i) = t{i}(t_settle_index);
end
figure;
plot(M2_values,t_settle_array, '-o');
title("Settling time vs load mass");
xlabel("load mass (kg)")
ylabel("Robot settling time (s)");
%%
clear; close all; clc;
M1 = 150;
M2 = 60;
k = 100;
F = 250;
mu = 0.05;
g = 9.81;

k_values = 20:20:200;

% Store simulation results
 
x2 = {};
v2 = {};
a2 = {};
legends = {};
t = {};

% Loop over each value of M2
for i = 1:length(k_values)
    k = k_values(i); % Set the current value of M2
    simOut = sim('lab1.slx');
    x2{i} = simOut.x2;
    v2{i} = simOut.v2;
    a2{i} = simOut.a2;
    t{i} = simOut.tout;
    legends{i} = ['k = ' num2str(k)];
end

figure;
hold on;
for i = 1:length(x2)
    plot(t{i}, x2{i}); 
end
xlim([0,30])
legend(legends, Location="southeast");
xlabel('Time (s)');
ylabel('Load displacement (m)');
title('Load displacement with different spring constant value');
saveas(gcf,'x2_2_multiple_spring.png', "png")

figure;
hold on;
for i = 1:length(v2)
    plot(t{i}, v2{i}); 
end
xlim([0,30])
legend(legends, Location="southeast");
xlabel('Time (s)');
ylabel('Load velocity (m/s)');
title('Load velocity with different spring constant value');
saveas(gcf,'x2dot_2_multiple_spring.png', "png")

figure;
hold on;
for i = 1:length(a2)
    plot(t{i}, a2{i}); 
end
hold off;
xlim([0,30])
legend(legends, Location="northeast");
xlabel('Time (s)');
ylabel('Load acceleration (m/s^2)');
title('Load acceleration with different spring constant value');
saveas(gcf,'x2dotdot_pos_2_multiple_spring.png', "png")