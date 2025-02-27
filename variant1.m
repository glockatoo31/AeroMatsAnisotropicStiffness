clear; clc;

%% Input
E11   = 76e9;         % Pa (76 GPa)
E22   = 5.5e9;        % Pa (5.5 GPa)
G12   = 2.1e9;        % Pa (2.1 GPa)
nu12  = 0.34;         % unitless

theta_deg_single = 65.5;   % degrees

sigma_xx = 1730e6;    % Pa (1730 MPa)
sigma_yy = 683e6;     % Pa (683 MPa)
tau_xy   = 314e6;     % Pa (314 MPa)

% Pack stress vector (3x1)
sigma_x = [sigma_xx; tau_xy; sigma_yy];

%% Compliance matrix in 1-2 coords
S1 = [  1/E11,      -nu12/E11,       0;
       -nu12/E11,    1/E22,          0;
        0,           0,            1/G12 ];

%% Function for Variant 1 transformations
computeVariant1 = @(thDeg) doVariant1(S1, sigma_x, thDeg);

%% Single calculation at chosen angle (65.5 deg)
[exx_single, eyy_single, gxy_single, ...
 e11_single, e22_single, g12_single, ...
 s11_single, s22_single, t12_single] = computeVariant1(theta_deg_single);

disp('=== Variant 1: Single Angle ===');
fprintf('theta = %g deg\n', theta_deg_single);
fprintf('epsilon_xx  = %g\n', exx_single);
fprintf('epsilon_yy  = %g\n', eyy_single);
fprintf('gamma_xy    = %g\n', gxy_single);
fprintf('epsilon_11  = %g\n', e11_single);
fprintf('epsilon_22  = %g\n', e22_single);
fprintf('gamma_12    = %g\n', g12_single);
fprintf('sigma_11    = %g Pa\n', s11_single);
fprintf('sigma_22    = %g Pa\n', s22_single);
fprintf('tau_12      = %g Pa\n\n', t12_single);

%% Sweep theta from 0..90 for plots
thetaVec = 0:1:90;
nAng     = length(thetaVec);

exx_array  = zeros(nAng,1);
eyy_array  = zeros(nAng,1);
gxy_array  = zeros(nAng,1);
e11_array  = zeros(nAng,1);
e22_array  = zeros(nAng,1);
g12_array  = zeros(nAng,1);
s11_array  = zeros(nAng,1);
s22_array  = zeros(nAng,1);
t12_array  = zeros(nAng,1);

for i = 1:nAng
    th = thetaVec(i);
    [exx, eyy, gxy, e11, e22, g12, s11, s22, t12] = computeVariant1(th);
    exx_array(i) = exx;
    eyy_array(i) = eyy;
    gxy_array(i) = gxy;
    e11_array(i) = e11;
    e22_array(i) = e22;
    g12_array(i) = g12;
    s11_array(i) = s11;
    s22_array(i) = s22;
    t12_array(i) = t12;
end

%% Plot: Local Strains
figure('Name','Variant 1: Local Strains');
plot(thetaVec, e11_array, '-o', thetaVec, e22_array, '-s', thetaVec, g12_array, '-^','LineWidth',1.2);
xlabel('$\theta$ (deg)','Interpreter','latex');
ylabel('Local Strain (1-2)','Interpreter','latex');
legend({'$\epsilon_{11}$','$\epsilon_{22}$','$\gamma_{12}$'},'Interpreter','latex','Location','best');
title('Local Strains in 1-2 vs. $\theta$','Interpreter','latex');
hold on;
labelMinMax(thetaVec, e11_array);
labelMinMax(thetaVec, e22_array);
labelMinMax(thetaVec, g12_array);
hold off;

%% Plot: Global Strains
figure('Name','Variant 1: Global Strains');
plot(thetaVec, exx_array, '-o', thetaVec, eyy_array, '-s', thetaVec, gxy_array, '-^','LineWidth',1.2);
xlabel('$\theta$ (deg)','Interpreter','latex');
ylabel('Global Strain (x-y)','Interpreter','latex');
legend({'$\epsilon_{xx}$','$\epsilon_{yy}$','$\gamma_{xy}$'},'Interpreter','latex','Location','best');
title('Global Strains in x-y vs. $\theta$','Interpreter','latex');
hold on;
labelMinMax(thetaVec, exx_array);
labelMinMax(thetaVec, eyy_array);
labelMinMax(thetaVec, gxy_array);
hold off;

%% Plot: Local Stresses
figure('Name','Variant 1: Local Stresses');
plot(thetaVec, s11_array, '-o', thetaVec, s22_array, '-s', thetaVec, t12_array, '-^','LineWidth',1.2);
xlabel('$\theta$ (deg)','Interpreter','latex');
ylabel('Local Stress (Pa)','Interpreter','latex');
legend({'$\sigma_{11}$','$\sigma_{22}$','$\tau_{12}$'},'Interpreter','latex','Location','best');
title('Local Stresses in 1-2 vs. $\theta$','Interpreter','latex');
hold on;
labelMinMax(thetaVec, s11_array);
labelMinMax(thetaVec, s22_array);
labelMinMax(thetaVec, t12_array);
hold off;

%% Min/Max Summary (printed in Command Window)
disp('=== Variant 1: Min/Max over 0..90 deg ===');
printMinMax('epsilon_{11}', thetaVec, e11_array);
printMinMax('epsilon_{22}', thetaVec, e22_array);
printMinMax('gamma_{12}',   thetaVec, g12_array);
printMinMax('epsilon_{xx}', thetaVec, exx_array);
printMinMax('epsilon_{yy}', thetaVec, eyy_array);
printMinMax('gamma_{xy}',   thetaVec, gxy_array);
printMinMax('sigma_{11}',   thetaVec, s11_array);
printMinMax('sigma_{22}',   thetaVec, s22_array);
printMinMax('tau_{12}',     thetaVec, t12_array);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supporting functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [exx, eyy, gxy, e11, e22, g12, s11, s22, t12] = ...
    doVariant1(S1, sigma_x, thetaDeg)
    theta = deg2rad(thetaDeg);
    c = cos(theta);
    s = sin(theta);

    T1 = [ c^2,   s^2,   2*c*s;
           s^2,   c^2,  -2*c*s;
          -c*s,   c*s,   c^2 - s^2 ];
    T2_inv = T1;
    T2     = inv(T1);

    % Transform compliance
    Sx = T2_inv * S1 * T1;

    % Global strains
    eps_x = Sx * sigma_x;
    exx = eps_x(1);
    gxy = eps_x(2);
    eyy = eps_x(3);

    % Local stresses
    sigma_1 = T1 * sigma_x;
    s11 = sigma_1(1);
    t12 = sigma_1(2);
    s22 = sigma_1(3);

    % Local strains
    eps_1 = T2 * eps_x;
    e11 = eps_1(1);
    g12 = eps_1(2);
    e22 = eps_1(3);
end

function printMinMax(labelStr, thetaArr, dataArr)
    [mi, iMin] = min(dataArr);
    [ma, iMax] = max(dataArr);
    fprintf(' %s: min=%g @%g deg, max=%g @%g deg\n',...
        labelStr, mi, thetaArr(iMin), ma, thetaArr(iMax));
end

function labelMinMax(thetaArr, dataArr)
    % Find min and max
    [valMin, idxMin] = min(dataArr);
    [valMax, idxMax] = max(dataArr);
    % Plot markers for min and max, hide from legend
    plot(thetaArr(idxMin), valMin, 'kd', 'MarkerFaceColor','k', ...
        'HandleVisibility','off');
    plot(thetaArr(idxMax), valMax, 'kd', 'MarkerFaceColor','k', ...
        'HandleVisibility','off');
    % Add text labels with a box
    text(thetaArr(idxMin), valMin, sprintf(' Min=%.3g', valMin), ...
        'HorizontalAlignment','left','VerticalAlignment','top', ...
        'BackgroundColor','white','EdgeColor','black','Interpreter','latex');
    text(thetaArr(idxMax), valMax, sprintf(' Max=%.3g', valMax), ...
        'HorizontalAlignment','left','VerticalAlignment','bottom', ...
        'BackgroundColor','white','EdgeColor','black','Interpreter','latex');
end
