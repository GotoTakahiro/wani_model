l_CFL0 = 0.35;
alpha = 0.01;
t_dec = 9;

t = 0:0.01:13;

for i = 1:length(t)
    if 0 <= t(i) && t(i) < 1
        l_CFL(i) = l_CFL0;
    elseif 1 <= t(i) && t(i) < t_dec
        l_CFL(i) = l_CFL0 - alpha*(t(i)-1);
    elseif t_dec <= t(i) && t(i) <= 13
        l_CFL(i) = l_CFL0 - alpha*(t_dec-1) - alpha*(1-exp(-(t(i)-t_dec)));
    end
end

figure;
plot(t, l_CFL, 'LineWidth', 2);
xlabel('Time [s]', 'FontSize', 25);
ylabel(['$l_{\mathrm{CFL_{ref}}}$', ' \textsf{[m]}'], 'Interpreter', 'latex', 'FontSize', 25);
xlim([0 13]);
ylim([0.25 0.37]);
h_axes = gca;
h_axes.FontSize = 20;
h_axes.XAxis.FontSize = 20;
h_axes.YAxis.FontSize = 20;
exportgraphics(gca, 'CFL_length.pdf', 'ContentType', 'vector');