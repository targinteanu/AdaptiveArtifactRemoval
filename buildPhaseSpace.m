function [phaseSpace, av] = buildPhaseSpace(filelabel, x, t, window, showPlotBool)
 %{
    showPlotBool=true;
    x=[1,2,4,8,16];
    t=[1,2,3,4,5];
 %}
    x=movmean(x,window);
    deriv = diff(x)./diff(t);
    phaseSpace = [x(1:(end-1));deriv]';
    [k,av] = convhull(phaseSpace, 'Simplify', true);
    if showPlotBool == true
       plot(phaseSpace(:,1),phaseSpace(:,2),".")
       title(filelabel); xlabel('x'); ylabel('dx/dt');
       hold on
       plot(phaseSpace(k,1),phaseSpace(k,2))
    end
    av
% end
%% 
%{
x=[1, 2, 1, 8, 4, 32, 64];
t=1:length(x);

data=[x',t'];

% Extract time intervals and data values
sample_rate=100; %samples per second
time = data(:, 1);
%time = 0:length(data)/sample_rate:1/sample_rate
data_values = data(:, 2); 

% Calculate phase space
phase_space = [data_values(1:end-1), data_values(2:end)];

% Plot phase space
figure;
plot(phase_space(:, 1), phase_space(:, 2), '-o');
xlabel('Data(t)');
ylabel('Data(t+1)');
title('Phase Space Plot');
grid on;
%}

end