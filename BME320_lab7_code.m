raw_data = (xlsread('BME320_lab7_sc_data.csv'))';
Qb = 0.1398;
len_t1 = 720;
len_t2 = 600;
len_t3 = 480;

time_t1 = linspace(1,len_t1,len_t1);
time_t2 = linspace(1,len_t2,len_t2);
time_t3 = linspace(1,len_t3,len_t3);
for i = 1:len_t1
    W(1,i) = Qb * (raw_data(1,1) - raw_data(1,i));
end
for i = 1:len_t2
    W(2,i) = Qb * (raw_data(2,1) - raw_data(2,i));
end
for i = 1:len_t3
    W(3,i) = Qb * (raw_data(3,1) - raw_data(3,i));
end
xinterval = 20;
lw = 2;
xq_1 = (0:xinterval:len_t1);
interp_W1 = spline(time_t1,W(1,1:len_t1),xq_1);
xq_2 = (0:xinterval:len_t2);
interp_W2 = spline(time_t2,W(2,1:len_t2),xq_2);
xq_3 = (0:xinterval:len_t3);
interp_W3 = spline(time_t3,W(3,1:len_t3),xq_3);

figure(1)
subplot(3,1,1)
p1 = plot(time_t1,W(1,1:len_t1),xq_1,interp_W1)
p1(1).LineWidth = lw;
p1(2).LineWidth = lw;
title('Dialysate pump flow rate: 130 mL/min'); ...
    xlim([0 800]); ...
    ylim([-0.001 .02]);...
    xlabel('Time (sec)');...
    ylabel('Total Mass Transferred (mol / min)');...
    legend('Data Values', 'Line of Best Fit');
    

subplot(3,1,2)
p2 = plot(time_t2,W(2,1:len_t2),xq_2,interp_W2)
p2(1).LineWidth = lw;
p2(2).LineWidth = lw;
title('Dialysate pump flow rate: 260 mL/min'); ...
    xlim([0 800]); ...
    ylim([-0.001 .02]);...
    xlabel('Time (sec)');...
    ylabel('Total Mass Transferred (mol / min)');...
    legend('Data Values', 'Line of Best Fit');

subplot(3,1,3)
p3 = plot(time_t3,W(3,1:len_t3),xq_3,interp_W3)
p3(1).LineWidth = lw;
p3(2).LineWidth = lw;
title('Dialysate pump flow rate: 390 mL/min'); ...
    xlim([0 800]); ...
    ylim([-0.001 .02]);...
    xlabel('Time (sec)');...
    ylabel('Total Mass Transferred (mol / min)');...
    legend('Data Values', 'Line of Best Fit');

