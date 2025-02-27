clear; clc;

%% Input
E11 = 76e9;          % Pa (76 GPa)
E22 = 5.5e9;         % Pa (5.5 GPa)
G12 = 2.1e9;         % Pa (2.1 GPa)
nu12 = 0.34;         % unitless

theta_deg_single = 65.5;  % degrees

sigma_xx = 1730e6;   % Pa (1730 MPa)
sigma_yy = 683e6;    % Pa (683 MPa)
tau_xy   = 314e6;    % Pa (314 MPa)

sigma_x = [sigma_xx; tau_xy; sigma_yy];

%% Compliance & stiffness
S1 = [  1/E11,      -nu12/E11,       0;
       -nu12/E11,    1/E22,          0;
        0,           0,            1/G12 ];
C1 = inv(S1);

%% Variant 2 transformations
% [Cx] = T1^-1 * C1 * T2, then eps_x = inv(Cx)*sigma_x
computeVariant2 = @(thDeg) doVariant2(C1, sigma_x, thDeg);

%% Single calculation at chosen angle
[exx_single, eyy_single, gxy_single, ...
 e11_single, e22_single, g12_single, ...
 s11_single, s22_single, t12_single] = computeVariant2(theta_deg_single);

disp('=== Variant 2: Single Angle ===');
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

%% Sweep theta from 0..90
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
    [exx, eyy, gxy, e11, e22, g12, s11, s22, t12] = computeVariant2(th);
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

%% Plots: Strains
figure('Name','Variant 2: Strains');
subplot(2,1,1);
plot(thetaVec,e11_array,'-o',thetaVec,e22_array,'-s',thetaVec,g12_array,'-^','LineWidth',1.2);
xlabel('$\theta$ (deg)','Interpreter','latex');
ylabel('Local Strain (1-2)','Interpreter','latex');
legend({'$\epsilon_{11}$','$\epsilon_{22}$','$\gamma_{12}$'}, ...
       'Interpreter','latex','Location','best');
title('Strains in 1-2 vs. $\theta$','Interpreter','latex');

subplot(2,1,2);
plot(thetaVec,exx_array,'-o',thetaVec,eyy_array,'-s',thetaVec,gxy_array,'-^','LineWidth',1.2);
xlabel('$\theta$ (deg)','Interpreter','latex');
ylabel('Global Strain (x-y)','Interpreter','latex');
legend({'$\epsilon_{xx}$','$\epsilon_{yy}$','$\gamma_{xy}$'}, ...
       'Interpreter','latex','Location','best');
title('Strains in x-y vs. $\theta$','Interpreter','latex');

%% Plots: Local Stresses
figure('Name','Variant 2: Local Stresses');
plot(thetaVec,s11_array,'-o',thetaVec,s22_array,'-s',thetaVec,t12_array,'-^','LineWidth',1.2);
xlabel('$\theta$ (deg)','Interpreter','latex');
ylabel('Local Stress (Pa)','Interpreter','latex');
legend({'$\sigma_{11}$','$\sigma_{22}$','$\tau_{12}$'}, ...
       'Interpreter','latex','Location','best');
title('Stresses in 1-2 vs. $\theta$','Interpreter','latex');

%% Min/Max
disp('=== Variant 2: Min/Max over 0..90 deg ===');
printMinMax('epsilon_{11}',thetaVec,e11_array);
printMinMax('epsilon_{22}',thetaVec,e22_array);
printMinMax('gamma_{12}',  thetaVec,g12_array);
printMinMax('epsilon_{xx}',thetaVec,exx_array);
printMinMax('epsilon_{yy}',thetaVec,eyy_array);
printMinMax('gamma_{xy}',  thetaVec,gxy_array);
printMinMax('sigma_{11}',  thetaVec,s11_array);
printMinMax('sigma_{22}',  thetaVec,s22_array);
printMinMax('tau_{12}',    thetaVec,t12_array);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supporting functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [exx, eyy, gxy, e11, e22, g12, s11, s22, t12] = ...
    doVariant2(C1, sigma_x, thetaDeg)

    theta = deg2rad(thetaDeg);
    c = cos(theta);
    s = sin(theta);

    T1 = [ c^2,   s^2,   2*c*s;
           s^2,   c^2,  -2*c*s;
          -c*s,   c*s,   c^2 - s^2 ];
    T2 = inv(T1);

    % build Cx
    Cx = inv(T1) * C1 * T2;

    % global strains
    eps_x = inv(Cx) * sigma_x;
    exx   = eps_x(1);
    gxy   = eps_x(2);
    eyy   = eps_x(3);

    % local stresses
    sigma_1 = T1 * sigma_x;
    s11     = sigma_1(1);
    t12     = sigma_1(2);
    s22     = sigma_1(3);

    % local strains
    eps_1 = T2 * eps_x;
    e11   = eps_1(1);
    g12   = eps_1(2);
    e22   = eps_1(3);
end

function printMinMax(labelStr,thetaArr,dataArr)
    [mi,iMin] = min(dataArr);
    [ma,iMax] = max(dataArr);
    fprintf(' %s: min=%g @%g deg, max=%g @%g deg\n',...
        labelStr,mi,thetaArr(iMin),ma,thetaArr(iMax));
end
