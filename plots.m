function plots(t1, y1, t2, y2)
% PLOTS  Compares the solutions of two different umerical solvers 
% (for 2D states).

figure('Name','Numerical methods demo');
tiledlayout(2,2, 'Padding','compact','TileSpacing','compact');

% Phase portrait (EE)
nexttile; hold on; grid on; box on;
plot(y1(:,1), y1(:,2), 'LineWidth', 1.2);
xlabel('x'); ylabel('y'); title('Method 1 — phase portrait'); axis padded;

% Time series (EE)
nexttile; hold on; grid on; box on;
plot(t1, y1(:,1), 'LineWidth', 1.2);
plot(t1, y1(:,2), 'LineWidth', 1.2);
xlabel('t'); ylabel('state'); title('Method 1 — time series');
legend('x','y','Location','northeast'); ylim('padded')

% Phase portrait (IE)
nexttile; hold on; grid on; box on;
plot(y2(:,1), y2(:,2), 'LineWidth', 1.2);
xlabel('x'); ylabel('y'); title('Method 2 — phase portrait'); axis padded;

% Time series (IE)
nexttile; hold on; grid on; box on;
plot(t2, y2(:,1), 'LineWidth', 1.2);
plot(t2, y2(:,2), 'LineWidth', 1.2);
xlabel('t'); ylabel('state'); title('Method 2 — time series');
legend('x','y','Location','northeast'); ylim('padded')

end

