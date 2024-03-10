% ================================================== %
% Average Radial Distribution Computation %

r_mat_1 = average_calc('GFP3.xlsx','MCH3.xlsx',1,2,1,'1/29/20 GFP No peptide 1','1/29/20 mCherry No peptide 1',1);
r_mat_2 = average_calc('GFP6.xlsx','MCH6.xlsx',3,4,1,'1/29/20 GFP No peptide 2','1/29/20 mCherry No peptide 2',1);
r_mat_3 = average_calc('GFP9.xlsx','MCH9.xlsx',5,6,1,'1/29/20 GFP No peptide 3','1/29/20 mCherry No peptide 3',1);
r_mat_4 = average_calc('GFP13.xlsx','MCH13.xlsx',7,8,1,'2/3/20 GFP No peptide 1','2/3/20 mCherry No peptide 1',1);
r_mat_5 = average_calc('GFP16.xlsx','MCH16.xlsx',9,10,1,'2/3/20 GFP No peptide 2','2/3/20 mCherry No peptide 2',1);
r_mat_6 = average_calc('GFP19.xlsx','MCH19.xlsx',11,12,2,'2/3/20 GFP No peptide 3','2/3/20 mCherry No peptide 3',1);
% Stores values for radial distribution for each data set for later
% incorperation into average radially dependent distribution plot

GFP_ave = zeros(length(r_mat_1(:,1)),1);
MCH_ave = zeros(length(r_mat_1(:,1)),1);

A = cat(3,r_mat_1,r_mat_2,r_mat_3,r_mat_4,r_mat_5,r_mat_6);
standard_err = std(A,[],3,'omitnan')/sqrt(6);

for i = 1:length(r_mat_1(:,1))
    a = r_mat_1(i,2);
    b = r_mat_2(i,2);
    c = r_mat_3(i,2);
    a2 = r_mat_4(i,2);
    b2 = r_mat_5(i,2);
    c2 = r_mat_6(i,2);
    new_g = (a+b+c+a2+b2+c2)/6;

    d = r_mat_1(i,3);
    f = r_mat_2(i,3);
    g = r_mat_3(i,3);
    d2 = r_mat_4(i,3);
    f2 = r_mat_5(i,3);
    g2 = r_mat_6(i,3);
    new_m = (d+f+g+d2+f2+g2)/6;

    GFP_ave(i) = new_g;
    MCH_ave(i) = new_m;
end
% The above for loop calculates average cell fractions at each radius
% considered

figure(13)
scatter(r_mat_1(:,1),GFP_ave,15,'r','filled')
hold on
scatter(r_mat_1(:,1),MCH_ave,15,'b','filled')
errorbar(r_mat_1(:,1),GFP_ave,standard_err(:,2),'k','LineStyle','none')
errorbar(r_mat_1(:,1),MCH_ave,standard_err(:,3),'k','LineStyle','none')
title('Average Radially Dependent Cell Distribution (No Protein)')
xlabel('Radius Considered')
ylabel('Fraction of Cells Encompassed')
legend('GFP', 'mCherry', 'Location', 'northwest')
hold off

function rad_mat = average_calc(gfpf, mchf, a, b, art_num, pro_type_title_1, pro_type_title_2, art_num_m)
GFP_table = readtable(gfpf,'Sheet','DATA','PreserveVariableNames',true);
mCh_table = readtable(mchf,'Sheet','DATA','PreserveVariableNames',true);
% Convert excel data tables to matlab readable table

x = GFP_table.Var1; % Access column one of the excel table (containing x coordinates for GFP cells)
y = GFP_table.Var2; % Access column two of the excel table (containing y coordinates for GFP cells)

for p = 1:art_num
    ave_x = mean(x);
    ave_y = mean(y);
    outlier_index = 0;
    outlier_dist = 0;

    for i = 1:length(x)
        d = sqrt((x(i)-ave_x)^2 + (y(i)-ave_y)^2);
        if d > outlier_dist
            outlier_dist = d;
            outlier_index = i;
        end
    end

    x(outlier_index) = [];
    y(outlier_index) = [];
end

% Average of points on one graph

% Check for distance away from average cell position (approx. center of the
% microscope view) D = sqrt((x-xi)^2 + (y-yi)^2)

x_max = max(x); % Find and store x coordinate for right-most cell
y_max = max(y); % Find and store y coordinate for upper-most cell
x_min = min(x); % Find and store x coordinate for left-most cell
y_min = min(y); % Find and store y coordinate for lower-most cell

diff_x = x_max - x_min; % Compute and store difference between x max and x min
diff_y = y_max - y_min; % Compute and store difference between y max and y min

rad_x = diff_x/2; % Compute and store radius as determined by diff_x
rad_y = diff_y/2; % Compute and store radius as determined by diff_y

x = x - (rad_x + x_min); % Center all x coordinates about y = 0
y = y - (rad_y + y_min); % Center all y coordinates about x = 0

x = x/(max(x)); % Compress all points horizontally to result in a horizontal radius of 1.0
y = y/(max(y)); % Compress all points vertically to result in a vertical radius of 1.0

% Below is the same procedure as above but for the corresponding mCherry
% cells in the same well
xm = mCh_table.Var1;
ym = mCh_table.Var2;

for p = 1:art_num_m
    ave_x = mean(xm);
    ave_y = mean(ym);
    outlier_index = 0;
    outlier_dist = 0;

    for i = 1:length(xm)
        d = sqrt((xm(i)-ave_x)^2 + (ym(i)-ave_y)^2);
        if d > outlier_dist
            outlier_dist = d;
            outlier_index = i;
        end
    end

    xm(outlier_index) = [];
    ym(outlier_index) = [];
end

xm_max = max(xm);
ym_max = max(ym);
xm_min = min(xm);
ym_min = min(ym);

diff_xm = xm_max - xm_min;
diff_ym = ym_max - ym_min;

rad_xm = diff_xm/2;
rad_ym = diff_ym/2;

xm = xm - (rad_xm + xm_min);
ym = ym - (rad_ym + ym_min);

xm = xm/(max(xm));
ym = ym/(max(ym));

% Plot heat map cell distribution for GFP cells
figure(a)
scatplot(x,y)
title(pro_type_title_1);
xlabel('X Position');
ylabel('Y Position');

% Plot heat map cell distribution for the mCherry cells
figure(b)
scatplot(xm,ym)
title(pro_type_title_2);
xlabel('X Position');
ylabel('Y Position');

% Add conditional if to check for other cells around outlier cells
% CellProfiler pipeline ==> check for outliers

% End of Heat Map Coding Portion
% ==================================================== %
% ==================================================== %
% Radially Dependent Cell Concentration Analysis

num_points = 40;
y_densities = zeros(num_points,1);
ym_densities = zeros(num_points,1);
A = zeros(num_points,1);

for j = 1:num_points
    r = j/num_points;
    A(j) = pi * r^2;
    for k = 1:length(x)
        distance = sqrt((x(k))^2 + (y(k))^2);
        if distance <= r
            y_densities(j) = y_densities(j) + 1;
        end
    end
    for k = 1:length(xm)
        distance2 = sqrt((x(k))^2 + (y(k))^2);
        if distance2 <= r
            ym_densities(j) = ym_densities(j) + 1;
        end
    end
end

y_densities(num_points) = length(x);
y_densities = y_densities/length(x);
ym_densities(num_points) = length(xm);
ym_densities = ym_densities/length(xm);
r_mat = sqrt(A/pi);

rad_mat = zeros(length(r_mat),3);
rad_mat(:,1) = r_mat;
rad_mat(:,2) = y_densities;
rad_mat(:,3) = ym_densities;

% Instantaneous Slope = Cell Density (in terms of circular area centered at
% origin)
end