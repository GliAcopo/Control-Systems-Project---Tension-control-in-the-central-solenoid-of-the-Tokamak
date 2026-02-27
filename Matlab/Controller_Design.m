%% Clear
clc
clear all
close all

%% Start settings
% Dock all figures -> Does not work
set(0, "DefaultFigureWindowStyle", "docked");
set(groot, 'defaultFigureWindowStyle', 'docked');

% Permanently activate "grid on"
set(groot,'DefaultAxesXGrid','on','DefaultAxesYGrid','on','DefaultAxesZGrid','on');

% format long g % printing to maximum "double" precision
format default;

% To save pictures in the latex folder
figFolder = fullfile('..', 'Figures'); % goes up one level, then into Figures

% activate figure saving into files
active = false;

enablesave = false;

% Number of figure
figcounter = 1;
%% Variables and measures
s = tf('s');

% Reference
rho_r = 0;
r = 1/(s^rho_r);

rp = 1/s^3; % parable reference
rr = 1/s^2; % ramp reference
rs = 1/s^1; % step reference
rc = 1;     % Constant reference

% Disturbances


% Measuring units
milli = 1e-3; 
micro = 1e-6;

% Variables
coupling_coefficient = 0.8;
L_P = 1.9 * micro;  % Henry
L_T = 80 * milli;
R_T = 4 * micro;    % OHM
R_P = 1 * milli;

M = coupling_coefficient*sqrt(L_T*L_P);

% System P(s) = (y(s))/(u(s)) = (I_P)/(v_T) = (-M*s)/((R_T+L_T*s)*(R_P+L_P*s)-M^2*s^2)
P = minreal( ...
    (-M*s)/((R_T+L_T*s)*(R_P+L_P*s)-M^2*s^2) ...
    );

%% System elements

k = 1;     % System gain

rho = 3; % Number of poles in the Controller

C  = (k/s^rho);     % Controller
H  = 1;    % Sensor (Filter on measure)
Ff = 0;    % Feed Forward
Fr = 1;    % Filter on reference


L = minreal(C*P*H)
wyr = minreal((P*Ff+P*C*Fr)/(1+P*C*H))

%% Analisys of the process
gainP = dcgain(P)
% The gain of the system is zero, so this means that it is not a
% non-minimum phase system
polesP = pole(P)
% We have no unstable poles
zerosP = zero(P)
% we have one zero in the origin
%% Steady state performance

% Thanks to TVF we get that if we want perfect asymptotic following of:
% Because in the process P we have a derivator, we need to account for 
% rho + 1 integrators
% 1. Parable -> rho = 4
% 2. Ramps -> rho = 3
% 3. Step -> rho = 2

% We verify the TVF, if the system is able to follow the reference with
% limited error, then the dcgain will be different from inf

rho_r = 0;
r = 1/(s^rho_r);
lim = s * r * (1/(1+L));
followcostant = dcgain(lim) == 0

rho_r = 1;
r = 1/(s^rho_r);
lim = s * r * (1/(1+L));
followstep = dcgain(lim) == 0

rho_r = 2;
r = 1/(s^rho_r);
lim = s * r * (1/(1+L));
followramp = dcgain(lim) == 0

rho_r = 3;
r = 1/(s^rho_r);
lim = s * r * (1/(1+L));
followparable = dcgain(lim) == 0

%% Transient performances
%% Analysis STARTING OBSERVATION
close all
clc;
k = 1;     % System gain

C = minreal(k/s^rho);     % Controller
H = 1;     % Sensor
Ff = 1;    % Feed Forward
Fm = 1;    % Filter on Measure
Fr = 1;    % Filter on reference

L = minreal(C*P*H)
wyr = minreal((P*Ff+P*C*Fr)/(1+P*C*H))


stepinfowyr  = stepinfo(wyr)
allmarginwyr = allmargin(wyr)
isstablewyr  = isstable(wyr) 
% Figures
figure(1)
    rlocus(L);
    figcounter = 1;
    figName = sprintf('%d-Rlocus - Loop function of base system - positive gain', figcounter);
    if active; SaveFigure(figName, figFolder, true, true); end
figure(2)
    margin(L);
    figcounter = 2;
    figName = sprintf('%d-Bode - Loop function of base system', figcounter);
    figcounter = figcounter + 1;
    if active; SaveFigure(figName, figFolder, true, true); end
figure(3)
    step(wyr);
    figcounter = 3;
    figName = sprintf('%d-Step - closed Loop function of base system', figcounter);
    figcounter = figcounter + 1;
    if active; SaveFigure(figName, figFolder, true, true); end
 figure(4)
    rlocus(-1*L);
    figcounter = 4;
    figName = sprintf('%d-Rlocus - Loop function of base system - negative gain', figcounter);
    if active; SaveFigure(figName, figFolder, true, true); end

% A first analysis shows us that the system is unstable, with a zero in 
% the orgin and three poles in the origin we will need to handle attract
% the rlocus with multiple zeroes near the origin
% With three poles in the controller we can add up to three zeroes before
% the C is not phisically achievable
%% -----------------------------------------------------------------------
%                                 ANALYSYS (1)
% ------------------------------------------------------------------------
% Load variables from workspace
%load("starting_workspace.mat")
% In this analisys we precisely delete all the poles in the process P
P = zpk( ...
    (-M*s)/...
    (((R_T+L_T*s)*(R_P+L_P*s)-M^2*s^2)) ...
    )
% Zero to cancel the poles of the system
cancellationC = (((R_T+L_T*s)*(R_P+L_P*s)-M^2*s^2));

%% .................................(1.A)..................................
% In this analysis we observe that the system is always marginally stable
% so there is no need for us to take any action other than choose the gain
% But this happens only with the negative rlocus
clc;
close all
% (1.A.I) POSITIVE RLOCUS
k = 1;     % System gain
C = minreal(cancellationC*(k/s^rho)) % Controller
biproperFCN = size(pole(C), 1) >= size(zero(C), 1)

% System components
H = 1;     % Sensor
Ff = 0;    % Feed Forward
Fr = 1;    % Filter on reference
L = minreal(C*P*H);
wyr = minreal((P*Ff+P*C*Fr)/(1+P*C*H));

figure(1)
    rlocus(L);
    axis equal
    if active; SaveFigure("5-Rlocus - Case (1.A.I)", figFolder, true, true); end

% (1.A.II) NEGATIVE RLOCUS
k = -1;     % System gain
C = minreal(cancellationC*(k/s^rho)) % Controller
biproperFCN = size(pole(C), 1) >= size(zero(C), 1)

% System components
H = 1;     % Sensor
Ff = 1;    % Feed Forward
Fm = 1;    % Filter on Measure
Fr = 1;    % Filter on reference
L = minreal(C*P*H);
wyr = minreal((P*Ff+P*C*Fr)/(1+P*C*H));

figure(2)
    rlocus(L);
    axis equal
    if active; SaveFigure("6-Rlocus - Case (1.A.II)", figFolder, true, true); end

stepinfowyr  = stepinfo(wyr)
allmarginwyr = allmargin(wyr)
isstablewyr  = isstable(wyr)

figure(3)
    margin(L);
    grid on
figure(4)
    step(wyr);
    axis equal
    if active; SaveFigure("7-Step Response - Case (1.A.II)", figFolder, true, true); end


%% .................................(1.B)..................................
% We place a zero in order to GUARANTEE stability of the system
clc;
close all

%{
re = -1; % real part of complex zero
im = 1;  % Imaginary part of complex zero
% IT IS BETTER TO NOT USE COMPLEX ZEROES SINCE THEY SLOW THE RESPONSE DOWN
zeroesC = (s^2 + -2*re*s + re^2 + im^2)*... % two complex zeroes
          (s + 1);                          % One fast zero (small real)
%}
zeroesC = (s + 1);
polesC  = 1%(s + 5); % Fat pole to compensate for the zeroes (biproper function needed)

k = -64.371;     % System gain

% Lag network
alpha = 100;
tau = 0.5;
Rr = (1+tau*alpha*s)/(1+tau*s);

C = minreal(Rr * cancellationC * 1/polesC * zeroesC * (k/s^rho)) % Controller
biproperFCN = size(pole(C), 1) >= size(zero(C), 1)                          % szdim2 = size(A,2) Interrogare la lunghezza della seconda dimensione di A.

% System components
H = 1;     % Sensor
Ff = 0;    % Feed Forward
           % Filter on reference
Fr  = 1;%/(1 + s/0.9)
L = minreal(C*P*H);
wyr = minreal((P*Ff+P*C*Fr)/(1+P*C*H));


stepinfowyr  = stepinfo(wyr)
%allmarginwyr = allmargin(wyr)
isstablewyr  = isstable(wyr)


figure(1)
    step(wyr); grid on;
    if active; SaveFigure("8-Step Response - Case (1.B)", figFolder, true, true); end
figure(2)
    margin(L); grid on;
    if active; SaveFigure("9-Bode L - Case (1.B)", figFolder, true, true); end
figure(3)
    rlocus(L); grid on;
    axis([-2.1 0.1 -0.03 0.03]); % Set [xmin xmax ymin ymax]
    if active; SaveFigure("10-Rlocus L - Case (1.B)", figFolder, true, true); end
figure(4)
    NyquistLog(L); grid on;
    % axis([-2.1 0.1 -0.03 0.03]); % Set [xmin xmax ymin ymax]
    if active; SaveFigure("11-Logarithmic Nyquist of L - Case (1.B)", figFolder, true, true); end
%% Robustness analisys of -> Analysis (1)
%% Initial sample
% Save old variables
close all
M_old = M;
P_old = P;
L_old = L;
wyr_old = wyr;

M = M * 1.20;

P = minreal( (-M * s) / (((R_T + L_T * s) * (R_P + L_P * s) -M^2 * s^2)) );
L = minreal(C * P * H);
wyr = minreal( (P * Ff + P * C * Fr)/(1 + P * C * H) );

stepinfowyr_old = stepinfo(wyr_old); % WARNING, RUN PRECEDENT SECTION TO CALCULATE WYR
stepinfowyr  = stepinfo(wyr);
% Differences in percentages
pct_settling = (stepinfowyr.SettlingTime - stepinfowyr_old.SettlingTime) ...
               / stepinfowyr_old.SettlingTime * 100

pct_overshoot = (stepinfowyr.Overshoot - stepinfowyr_old.Overshoot) ...
                / stepinfowyr_old.Overshoot * 100

pct_rise = (stepinfowyr.RiseTime - stepinfowyr_old.RiseTime) ...
           / stepinfowyr_old.RiseTime * 100


figure (1); step(wyr_old); title('Response before change in M');
figure (2); step(wyr); title('Response after change in M');
figure (3); margin(L_old); title('Stability margins before change in M');
figure (4); margin(L); title('Stability margins after change in M');


M = M_old;
P = P_old;
L = L_old;
wyr = wyr_old;
%% ------------------------------------------------------------------------
% Thorough analysis of Robustness of -> Analysis (1)
close all
s = tf('s');
P_old = P; L_old = L; C_old = C; H_old = H; Ff_old = Ff; Fr_old = Fr;
wyr_old = wyr;
M_old = M; R_T_old = R_T; R_P_old = R_P; L_T_old = L_T; L_P_old = L_P;

% Choose what variables to change
changeM = true;
changeR_T = true;
changeR_P = false;
changeL_T = false;
changeL_P = false;
% Figures names

    try
        % Build list of changed variable short names
        changed = {};
        if changeM;   changed{end+1} = 'M';   end
        if changeR_T; changed{end+1} = 'R_T'; end
        if changeR_P; changed{end+1} = 'R_P'; end
        if changeL_T; changed{end+1} = 'L_T'; end
        if changeL_P; changed{end+1} = 'L_P'; end
    
        % Create suffix and filename
        if isempty(changed)
            suffix = 'base';
        else
            suffix = strjoin(changed, '-');   % e.g. 'R_T' or 'M_R_T'
        end
        figName1 = sprintf('1-R(1.B)-Robustness_%s_step_response_of_(1.B)_system', suffix);
        figName2 = sprintf('2-R(1.B)-Robustness_%s_bode_plot_of_(1.B)_system', suffix);
    catch e %e is an MException struct
        fprintf(2,'ERROR at %s (line %d)\n', e.stack(1).name, e.stack(1).line);
        fprintf(2,'Identifier: %s\n', e.identifier);
        fprintf(2,'Message:\n%s\nCause:\n%s\n', e.message, char(join(string(arrayfun(@(c)string(c.message), e.cause, 'UniformOutput', true)), newline)));
        fprintf(2,'Correction:\n%s\n', char(join(string(e.Correction), newline)));
    end

try 
    num_of_colors = 15;
    % set(gca,'colororder',parula(num_of_colors)) % More colors for the graph lines yay!
    ax = gca;
    if isa(ax,'matlab.graphics.axis.Axes')
        set(ax,'ColorOrder',parula(num_of_colors));
    end
    
    stepinfowyr_old = stepinfo(wyr_old) % WARNING, RUN PRECEDENT SECTION TO CALCULATE WYR
    
    [y,t] = step(wyr_old);
    figure(1); clf;
    plot(t, squeeze(y), 'LineWidth', 2);
    hold on; grid on;
    leg = {'nominal'};
    
    for st = 0.50 : 0.25 : 1.50 % st stands for step, but calling it "step" will overshadow the function name
        if st == 1.0; fprintf("jumping nominal\n"); continue; end % Jump nominal case
        if changeM;   M   = M_old   * st; end
        if changeR_T; R_T = R_T_old * st; end
        if changeR_P; R_P = R_P_old * st; end
        if changeL_T; L_T = L_T_old * st; end
        if changeL_P; L_P = L_P_old * st; end
    
        P = minreal( (-M * s) / (((R_T + L_T * s) * (R_P + L_P * s) -M^2 * s^2)) );
        wyr = minreal( (P * Ff + P * C * Fr)/(1 + P * C * H) );

        % Check if some system is unstable
            y = squeeze(step(wyr,t));
            if ~(isstable(wyr))
                fprintf(2,'st=%.2f exceeds ylim (max|y|=%.3g)\n', st, max(abs(y)));
                plot(t,y,'LineWidth', 2, 'Color', "red", "LineStyle", "--");
                ylim([0 1.3]);
            else
                plot(t,y,'LineWidth',1.7);
            end
            
        %plot(t, squeeze(y), 'LineWidth', 1.7);
        leg{end+1} = sprintf('%s %+0.0f%%', suffix, (st-1)*100);
    end
    legend(leg, 'Location','best');
    hold off;
    if active; SaveFigure(figName1, figFolder, true, true); end
    
    
    figure(2); clf;
    bp = bodeplot(L_old); hold on; grid on;
    bp.Responses(1).LineWidth = 1.5;
    % opt = getoptions(bp); opt.ColorOrder = parula(num_of_colors); setoptions(bp, opt);
    leg = {'nominal'};
    for st = 0.50 : 0.25 : 1.50
        if st == 1.0; fprintf("jumping nominal\n"); continue; end % Jump nominal case
        if changeM;   M   = M_old   * st; end
        if changeR_T; R_T = R_T_old * st; end
        if changeR_P; R_P = R_P_old * st; end
        if changeL_T; L_T = L_T_old * st; end
        if changeL_P; L_P = L_P_old * st; end
    
        P = minreal( (-M*s) / (((R_T + L_T*s)*(R_P + L_P*s) - M^2*s^2)) );
        L = minreal(C * P * H);
    
        bp = bodeplot(L);
        bp.Responses(1).LineWidth = 1.3;
        leg{end+1} = sprintf('%s %+0.0f%%', suffix, (st-1)*100);
    end
    legend(leg, 'Location','best');
    hold off;
    if active; SaveFigure(figName2, figFolder, true, true); end
catch e %e is an MException struct
    fprintf(2,'ERROR at %s (line %d)\n', e.stack(1).name, e.stack(1).line);
    fprintf(2,'Identifier: %s\n', e.identifier);
    fprintf(2,'Message:\n%s\nCause:\n%s\n', e.message, char(join(string(arrayfun(@(c)string(c.message), e.cause, 'UniformOutput', true)), newline)));
    fprintf(2,'Correction:\n%s\n', char(join(string(e.Correction), newline)));
    % We do not care about what we catch we only want to restore our
    % variables safely
end

M   = M_old;   R_T = R_T_old;   R_P = R_P_old;   L_T = L_T_old;   L_P = L_P_old;
P   = P_old; L   = L_old; H   = H_old; Ff  = Ff_old; Fr  = Fr_old;
wyr = wyr_old;

%% -----------------------------------------------------------------------
%                                 ANALYSYS (2)
% ------------------------------------------------------------------------
% Load variables from workspace
%load("starting_workspace.mat")
%% .................................(2.A)..................................
% With a positive gain we need three zeroes in order to stabilize the
% system
clc;
close all

rho = 3; % Number of poles in the Controller
k = 1;

zeroesC = (s+1)^3;
polesC  = 1;
C = minreal(1/polesC * zeroesC * (k/s^rho));     % Controller
if ~(size(pole(C), 1) >= size(zero(C), 1))
    fprintf(2, "Warning! C not biproper\n")
end

H = 1;     % Sensor
Ff = 0;    % Feed Forward
Fr = 1;    % Filter on reference

L = minreal(C*P*H)
wyr = minreal((P*Ff+P*C*Fr)/(1+P*C*H))


stepinfowyr  = stepinfo(wyr)
allmarginwyr = allmargin(wyr)
isstablewyr  = isstable(wyr)

figure(1)
    step(wyr); grid on;
    if active; SaveFigure("12-Step Response - Case (2.A)", figFolder, true, true); end
figure(2)
    margin(L); grid on;
    if active; SaveFigure("13-Bode L - Case (2.A)", figFolder, true, true); end
figure(3)
    rlocus(L); grid off;
    %axis([-2.1 0.1 -0.03 0.03]); % Set [xmin xmax ymin ymax]
    if active; SaveFigure("14-Rlocus L - Case (2.A)", figFolder, true, true); end
if active
    figure(4)
        NyquistLog(L); grid on;
        % axis([-2.1 0.1 -0.03 0.03]); % Set [xmin xmax ymin ymax]
        SaveFigure("15-Logarithmic Nyquist of L - Case (2.A)", figFolder, true, true);
    figure(5)
        rlocus(L); grid off;
        axis([-1.7 0.1 -0.04 0.04]); % Set [xmin xmax ymin ymax]
        if active; SaveFigure("14-Rlocus L - Case (2.A) ZOOM", figFolder, true, true); end
end

%% .................................(2.B)..................................
% Negative gain
% With a negative gain we are able to stabilize the system with just 2
% zeroes
clc;
close all

rho = 3; % Number of poles in the Controller
k = -1;

zeroesC = (s+1)^2;
polesC  = 1;
C = minreal(1/polesC * zeroesC * (k/s^rho));     % Controller
if ~(size(pole(C), 1) >= size(zero(C), 1))
    fprintf(2, "Warning! C not biproper\n")
end

H = 1;     % Sensor
Ff = 0;    % Feed Forward
Fr = 1;    % Filter on reference

L = minreal(C*P*H)
wyr = minreal((P*Ff+P*C*Fr)/(1+P*C*H))


stepinfowyr  = stepinfo(wyr)
allmarginwyr = allmargin(wyr)
isstablewyr  = isstable(wyr)

figure(1)
    step(wyr); grid on;
    if active; SaveFigure("12-Step Response - Case (2.B)", figFolder, true, true); end
figure(2)
    margin(L); grid on;
    if active; SaveFigure("13-Bode L - Case (2.B)", figFolder, true, true); end
figure(3)
    rlocus(L); grid off;
    %axis([-2.1 0.1 -0.03 0.03]); % Set [xmin xmax ymin ymax]
    if active; SaveFigure("14-Rlocus L - Case (2.B)", figFolder, true, true); end
if active
    figure(4)
        NyquistLog(L); grid on; hold off;
        % axis([-2.1 0.1 -0.03 0.03]); % Set [xmin xmax ymin ymax]
        SaveFigure("15-Logarithmic Nyquist of L - Case (2.B)", figFolder, true, true);
    figure(5)
        rlocus(L); grid off;
        axis([-3.5 0.5 -2 2]); % Set [xmin xmax ymin ymax]
        if active; SaveFigure("14-Rlocus L - Case (2.B) ZOOM", figFolder, true, true); end
end

%% .................................(2.B.I)................................
% Negative gain
% With a negative gain we are able to stabilize the system with just 2
% zeroes
close all
section = "2.B.I";
namestep = sprintf("12-Step Response - Case (%s)", section);
namebode = sprintf("13-Bode L - Case (%s)", section);
namerlocus = sprintf("14-Rlocus L - Case (%s)", section);
namenyquist = sprintf("15-Logarithmic Nyquist of L - Case (%s)", section);
namerlocuszoom = sprintf("14-Rlocus L - Case (%s) ZOOM", section);
clc;

rho = 3; % Number of poles in the Controller
k = -108; %% 77, 87 0.012, 106 0.01

zeroesC = (s+1)^2;
polesC  = 1;
C = minreal(1/polesC * zeroesC * (k/s^rho));     % Controller
if ~(size(pole(C), 1) >= size(zero(C), 1))
    fprintf(2, "Warning! C not biproper\n")
end

H = 1;     % Sensor
Ff = 0;    % Feed Forward
Fr = 1;    % Filter on reference

L = minreal(C*P*H)
wyr = minreal((P*Ff+P*C*Fr)/(1+P*C*H))


stepinfowyr  = stepinfo(wyr)
isstablewyr  = isstable(wyr)

figure(1)
    step(wyr); grid on;
    if active; SaveFigure(namestep, figFolder, true, true); end
figure(2)
    margin(L); grid on;
    if active; SaveFigure(namebode, figFolder, true, true); end
figure(3)
    rlocus(L); grid off;
    %axis([-2.1 0.1 -0.03 0.03]); % Set [xmin xmax ymin ymax]
    if active; SaveFigure(namerlocus, figFolder, true, true); end
if active
    figure(4)
        NyquistLog(L); grid on; hold off;
        % axis([-2.1 0.1 -0.03 0.03]); % Set [xmin xmax ymin ymax]
        SaveFigure(namenyquist, figFolder, true, true);
    figure(5)
        rlocus(L); grid off;
        axis([-3.5 0.5 -2 2]); % Set [xmin xmax ymin ymax]
        if active; SaveFigure(namerlocuszoom, figFolder, true, true); end
end
%% Disturbance tolerance
%load("ANALYSYS_2.B.I.mat")
%%
clc;
close all
%{
W_er   = 1/(1+P*C*H)           % Sensitivity function
Wyd1   = P/(1+P*C*H)
Wyd2   = -(P*C*H)/(1+P*C*H)    % Complementary sensitivity function

% More commonly known as:
L      = minreal(P*C*H);
S      = minreal(W_er);
T      = minreal(L/(1+L)); % abs(Wyd2)

% The other way around is:
W_er   = S;        % if Fr = 1
Wyd1   = minreal(P*S);
Wyd2   = -T;
Wymd2  = S;
%}

Wer   = minreal((Fr - H*P*Ff)/(1 + P*C*H));   % general
Wed1  = minreal(-(H*P)/(1 + P*C*H));
Wed2  = minreal(-H/(1 + P*C*H));

% Output paths ------------------------------------------------------------
Wyr   = minreal((P*Ff + P*C*Fr)/(1 + P*C*H));   % r -> y
Wyd1  = minreal(P/(1 + P*C*H));                 % d1 -> y
Wyd2  = minreal(-(P*C*H)/(1 + P*C*H));            % d2 -> y   (sensor disturbance)

dt = 1e-4;
t0 = 0;
t = t0 : dt : 0.5;

% d_2 ---------------------------------------------------------------------
A_2 = 1;
f = 60; % Hz
phi_2 = 0;
d2 = A_2*sin(2*pi*f*t + phi_2);

% d_1 ---------------------------------------------------------------------
% d1 as impulse approximation (one-sample pulse with area A1)
A1 = 0.05;

d1 = zeros(size(t));
idx = find(t >= t0, 1, 'first'); % Finds the first sample index whose time is at or after t0
d1(idx) = A1/dt;   % approximates A1*delta(t-t0)

% error lsim --------------------------------------------------------------
r = ones(size(t));
% Individual output contributions
y_r  = lsim(Wyr,  r,  t);
y_d1 = lsim(Wyd1, d1, t);
y_d2 = lsim(Wyd2, d2, t);

% Plot output response to sensor disturbance d2 only ----------------------
figure(6);
    plot(t, d2, 'LineWidth', 1.0); hold on;
    plot(t, y_d2, 'LineWidth', 1.2);
    grid on;
    legend('d_2(t)', 'y_{d2}(t)');
    xlabel('Time [s]');
    ylabel('Amplitude');
    named2 = sprintf('(%s) - Output response to sensor ripple disturbance d_2', section);
    title(named2);
    if active; SaveFigure(named2, figFolder, true, true); end

% Total output
y_total = y_r + y_d1 + y_d2;

figure(7);
    plot(t, y_total, 'LineWidth', 1.3); grid on;
    xlabel('Time [s]');
    ylabel('y(t)');
    title('Total output y(t) under r, d_1, d_2');
    nameytotal = sprintf('(%s) - Total output y(t) under r, d_1, d_2', section);
    if active; SaveFigure(nameytotal, figFolder, true, true); end

if ~isstable(wyr)
    warning('Closed-loop is unstable! Output/disturbance analysis is not valid.');
end
%%
save("Disturbance_tolerance_2.B.I.mat");


%% .................................(2.B.II)...............................
%load("Disturbance_tolerance_2.B.I.mat");
%%
% We add a notch filter in order to attenuate the response to the d2 error
close all
section = "2.B.II";
namestep = sprintf("12-Step Response - Case (%s)", section);
namebode = sprintf("13-Bode L - Case (%s)", section);
namerlocus = sprintf("14-Rlocus L - Case (%s)", section);
namenyquist = sprintf("15-Logarithmic Nyquist of L - Case (%s)", section);
namerlocuszoom = sprintf("16-Rlocus L - Case (%s) ZOOM", section);

rho = 3; % Number of poles in the Controller
k = -30;

zeroesC = (s+1)^2;
polesC  = 1;
C = minreal(1/polesC * zeroesC * (k/s^rho));     % Controller
if ~(size(pole(C), 1) >= size(zero(C), 1))
    fprintf(2, "Warning! C not proper\n")
end

% BUTTER
w1 = 2*pi*(f-5); % 50Hz
w2 = 2*pi*(f+5); % 60Hz
n = 2;
[b,a] = butter(n, [w1, w2], 'stop', 's')
H = tf(b, a)

Ff = 0;    % Feed Forward
Fr = 1;    % Filter on reference

L = minreal(C*P*H)
wyr = minreal((P*Ff+P*C*Fr)/(1+P*C*H))

stepinfowyr  = stepinfo(wyr)
isstablewyr  = isstable(wyr)

figure(1)
    step(wyr); grid on;
    if active; SaveFigure(namestep, figFolder, true, true); end
figure(2)
    margin(L); grid on;
    if active; SaveFigure(namebode, figFolder, true, true); end
figure(3)
    rlocus(L); grid off;
    %axis([-2.1 0.1 -0.03 0.03]); % Set [xmin xmax ymin ymax]
    if active; SaveFigure(namerlocus, figFolder, true, true); end
if active
    
        NyquistLog(L); grid on; hold off;
        % axis([-2.1 0.1 -0.03 0.03]); % Set [xmin xmax ymin ymax]
        SaveFigure(namenyquist, figFolder, true, true);
    figure(5)
        rlocus(L); grid off;
        axis([-3.5 0.5 -2 2]); % Set [xmin xmax ymin ymax]
        if active; SaveFigure(namerlocuszoom, figFolder, true, true); end
end

%% .................................(2.B.III)..............................
% We attenuate the error by decreasing the system gain crossover frequency
% well below 55Hz
close all
butteractive = true;
FfP1 = false;
FfP2 = true;
FfC = false;

section = "2.B.III-gain mod";
clc;
% Moving the crossover frequency away:
rho = 3; % Number of poles in the Controller
k = -1;

zeroesC = (s+2)^2;
polesC  = 1;
C = minreal(1/polesC * zeroesC * (k/s^rho));     % Controller
if ~(size(pole(C), 1) >= size(zero(C), 1))
    fprintf(2, "Warning! C not proper\n")
end

if butteractive
    w1 = 2*pi*(f-5); % 50Hz
    w2 = 2*pi*(f+5); % 60Hz
    n = 2;
    [b,a] = butter(n, [w1, w2], 'stop', 's')
    H = tf(b, a)
    section = section + "-butter";
else
    H = 1
end

if FfP1  
    % Approximate plant inverse
    Ff = 1/P % raw inverse
    
    % Make it proper + roll-off (essential if P has relative degree ≥1)
    wn = 50 * Wcp; % cutoff ~50× crossover (tune 20-100)
    zeta = 0.7;
    Ff = minreal(Ff * (wn^2 / (s^2 + 2*zeta*wn*s + wn^2)));  % 2nd-order low-pass
    if ~(size(pole(Ff), 1) >= size(zero(Ff), 1))
        fprintf(2, "Warning! Ff not proper\n")
    end
    section = section + "-FfP1";
elseif FfP2
    
    fc_ff = 500;                    % Hz (circa 500–1200)
    wn    = 2*pi*fc_ff;
    zeta  = 0.707;
    % Ff = 1/P + Low pass to make proper
    Ff    = minreal( (1/P) * (wn^2 / (s^2 + 2*zeta*wn*s + wn^2)) );
    if ~(size(pole(Ff), 1) >= size(zero(Ff), 1))
        fprintf(2, "Warning! Ff not proper\n")
    end
    section = section + "-FfP2";
elseif FfC
    Ff = C               % doubles the forward path → faster but more overshoot
    section = section + "-FfC";
else
    Ff = 0
end

Fr = 1;    % Filter on reference


L = minreal(C*P*H)
wyr = minreal((P*Ff+P*C*Fr)/(1+P*C*H));

stepinfowyr  = stepinfo(wyr)
[Gm,Pm,Wcg,Wcp] = margin(L);
Wcp = Wcp
WcpHz = Wcp/(2*pi)
isstablewyr  = isstable(wyr)

namestep = sprintf("12-Step Response - Case (%s)", section);
namebode = sprintf("13-Bode L - Case (%s)", section);
namerlocus = sprintf("14-Rlocus L - Case (%s)", section);
namenyquist = sprintf("15-Logarithmic Nyquist of L - Case (%s)", section);
namerlocuszoom = sprintf("16-Rlocus L - Case (%s) ZOOM", section);

figure(1)
    step(wyr); grid on;
    if active; SaveFigure(namestep, figFolder, true, true); end
figure(2)
    margin(L); grid on;
    if active; SaveFigure(namebode, figFolder, true, true); end
figure(3)
    rlocus(L); grid off;
    %axis([-2.1 0.1 -0.03 0.03]); % Set [xmin xmax ymin ymax]
    if active; SaveFigure(namerlocus, figFolder, true, true); end
if active
    
        NyquistLog(L); hold off;
        % axis([-2.1 0.1 -0.03 0.03]); % Set [xmin xmax ymin ymax]
        SaveFigure(namenyquist, figFolder, true, true);
    figure(5)
        rlocus(L); grid off;
        axis([-3.5 0.5 -2 2]); % Set [xmin xmax ymin ymax]
        if active; SaveFigure(namerlocuszoom, figFolder, true, true); end
end

if enablesave; save("2.B.III.mat"); end
%% DISTURBANCE TOLERANCE ANALYSYS 2.B.III
% Output paths ------------------------------------------------------------
close all
Wyr   = minreal((P*Ff + P*C*Fr)/(1 + P*C*H));   % r -> y
Wyd1  = minreal(P/(1 + P*C*H));                 % d1 -> y
Wyd2  = minreal(-(P*C*H)/(1 + P*C*H));            % d2 -> y   (sensor disturbance)

dt = 1e-4;
t0 = 0;
t = t0 : dt : 0.5;

% d_2 ---------------------------------------------------------------------
A_2 = 1;
f = 55; % Hz
phi_2 = 0;
d2 = A_2*sin(2*pi*f*t + phi_2);

% d_1 ---------------------------------------------------------------------
% d1 as impulse approximation (one-sample pulse with area A1)
A1 = 0.05;

d1 = zeros(size(t));
idx = find(t >= t0, 1, 'first'); % Finds the first sample index whose time is at or after t0
d1(idx) = A1/dt;   % approximates A1*delta(t-t0)

% error lsim --------------------------------------------------------------
r = ones(size(t));
% Individual output contributions
y_r  = lsim(Wyr,  r,  t);
y_d1 = lsim(Wyd1, d1, t);
y_d2 = lsim(Wyd2, d2, t);

% Plot output response to sensor disturbance d2 only ----------------------
figure(6);
    plot(t, d2, 'LineWidth', 1.0); hold on;
    plot(t, y_d2, 'LineWidth', 1.2);
    grid on;
    legend('d_2(t)', 'y_{d2}(t)');
    xlabel('Time [s]');
    ylabel('Amplitude');
    named2 = sprintf('(%s) - Output response to sensor ripple disturbance d_2', section);
    title(named2);
    if active; SaveFigure(named2, figFolder, true, true); end
hold off;
% Total output
y_total = y_r + y_d1 + y_d2;

figure(7);
    plot(t, y_total, 'LineWidth', 1.3); grid on;
    xlabel('Time [s]');
    ylabel('y(t)');
    title('Total output y(t) under r, d_1, d_2');
    nameytotal = sprintf('(%s) - Total output y(t) under r, d_1, d_2', section);
    if active; SaveFigure(nameytotal, figFolder, true, true); end
hold off;

w0 = 2*pi*f;
magWyd2 = abs(squeeze(freqresp(Wyd2, w0)));
fprintf('|Wyd2(j2pi%d)| = %.6f (%.2f dB)\n', f, magWyd2, 20*log10(magWyd2));

if ~isstable(wyr)
    warning('Closed-loop is unstable! Output/disturbance analysis is not valid.');
end

if enablesave; save("DISTURBANCE TOLERANCE ANALYSYS 2.B.III.mat"); end
%% ------------------------------------------------------------------------
% Thorough analysis of Robustness of -> Analysis (2.B.III)
close all
section = "2.B.III";
s = tf('s');
P_old = P; L_old = L; C_old = C; H_old = H; Ff_old = Ff; Fr_old = Fr;
wyr_old = wyr;
M_old = M; R_T_old = R_T; R_P_old = R_P; L_T_old = L_T; L_P_old = L_P;

% Choose what variables to change
changeM = true;
changeR_T = true;
changeR_P = false;
changeL_T = false;
changeL_P = false;
% Figures names

    try
        % Build list of changed variable short names
        changed = {};
        if changeM;   changed{end+1} = 'M';   end
        if changeR_T; changed{end+1} = 'R_T'; end
        if changeR_P; changed{end+1} = 'R_P'; end
        if changeL_T; changed{end+1} = 'L_T'; end
        if changeL_P; changed{end+1} = 'L_P'; end
    
        % Create suffix and filename
        if isempty(changed)
            suffix = 'base';
        else
            suffix = strjoin(changed, '-');   % e.g. 'R_T' or 'M_R_T'
        end
        figName1 = sprintf('1-R(%s)-Robustness_%s_step_response_of_(%s)_system', section, suffix, section);
        figName2 = sprintf('2-R(%s)-Robustness_%s_bode_plot_of_(%)_system', section, suffix, section);
    catch e %e is an MException struct
        fprintf(2,'ERROR at %s (line %d)\n', e.stack(1).name, e.stack(1).line);
        fprintf(2,'Identifier: %s\n', e.identifier);
        fprintf(2,'Message:\n%s\nCause:\n%s\n', e.message, char(join(string(arrayfun(@(c)string(c.message), e.cause, 'UniformOutput', true)), newline)));
        fprintf(2,'Correction:\n%s\n', char(join(string(e.Correction), newline)));
    end

try 
    num_of_colors = 15;
    % set(gca,'colororder',parula(num_of_colors)) % More colors for the graph lines yay!
    ax = gca;
    if isa(ax,'matlab.graphics.axis.Axes')
        set(ax,'ColorOrder',parula(num_of_colors));
    end
    
    stepinfowyr_old = stepinfo(wyr_old) % WARNING, RUN PRECEDENT SECTION TO CALCULATE WYR
    
    [y,t] = step(wyr_old);
    figure(1); clf;
    plot(t, squeeze(y), 'LineWidth', 2);
    hold on; grid on;
    leg = {'nominal'};
    
    for st = 0.50 : 0.25 : 1.50 % st stands for step, but calling it "step" will overshadow the function name
        if st == 1.0; fprintf("jumping nominal\n"); continue; end % Jump nominal case
        if changeM;   M   = M_old   * st; end
        if changeR_T; R_T = R_T_old * st; end
        if changeR_P; R_P = R_P_old * st; end
        if changeL_T; L_T = L_T_old * st; end
        if changeL_P; L_P = L_P_old * st; end
    
        P = minreal( (-M * s) / (((R_T + L_T * s) * (R_P + L_P * s) -M^2 * s^2)) );
        wyr = minreal( (P * Ff + P * C * Fr)/(1 + P * C * H) );

        % Check if some system is unstable
            YMAX = 2; 
            y = squeeze(step(wyr,t));
            if ~(isstable(wyr))
                fprintf(2,'st=%.2f exceeds ylim (max|y|=%.3g)\n', st, max(abs(y)));
                plot(t,y,'LineWidth', 2, 'Color', "red", "LineStyle", "--");
                ylim([-YMAX YMAX]);
            else
                plot(t,y,'LineWidth',1.7);
            end
            
        %plot(t, squeeze(y), 'LineWidth', 1.7);
        leg{end+1} = sprintf('%s %+0.0f%%', suffix, (st-1)*100);
    end
    legend(leg, 'Location','best');
    hold off;
    if active; SaveFigure(figName1, figFolder, true, true); end
    
    
    figure(2); clf;
    bp = bodeplot(L_old); hold on; grid on;
    bp.Responses(1).LineWidth = 1.5;
    % opt = getoptions(bp); opt.ColorOrder = parula(num_of_colors); setoptions(bp, opt);
    leg = {'nominal'};
    for st = 0.50 : 0.25 : 1.50
        if st == 1.0; fprintf("jumping nominal\n"); continue; end % Jump nominal case
        if changeM;   M   = M_old   * st; end
        if changeR_T; R_T = R_T_old * st; end
        if changeR_P; R_P = R_P_old * st; end
        if changeL_T; L_T = L_T_old * st; end
        if changeL_P; L_P = L_P_old * st; end
    
        P = minreal( (-M*s) / (((R_T + L_T*s)*(R_P + L_P*s) - M^2*s^2)) );
        L = minreal(C * P * H);
    
        bp = bodeplot(L);
        bp.Responses(1).LineWidth = 1.3;
        leg{end+1} = sprintf('%s %+0.0f%%', suffix, (st-1)*100);
    end
    legend(leg, 'Location','best');
    hold off;
    if active; SaveFigure(figName2, figFolder, true, true); end
catch e %e is an MException struct
    fprintf(2,'ERROR at %s (line %d)\n', e.stack(1).name, e.stack(1).line);
    fprintf(2,'Identifier: %s\n', e.identifier);
    fprintf(2,'Message:\n%s\nCause:\n%s\n', e.message, char(join(string(arrayfun(@(c)string(c.message), e.cause, 'UniformOutput', true)), newline)));
    fprintf(2,'Correction:\n%s\n', char(join(string(e.Correction), newline)));
    % We do not care about what we catch we only want to restore our
    % variables safely
end

M   = M_old;   R_T = R_T_old;   R_P = R_P_old;   L_T = L_T_old;   L_P = L_P_old;
P   = P_old; L   = L_old; H   = H_old; Ff  = Ff_old; Fr  = Fr_old;
wyr = wyr_old;


%% Verify system performaces
wyrstepinfo = stepinfo(wyr)
stepinfo(wyr, "SettlingTimeThreshold", 0.05)

%% Max delay
[Gm,Pm,Wcg,Wcp] = margin(L)
wgc = 4.63;
Njw = evalfr(H, 1j*wgc)      % N(jwgc)
phiN = angle(Njw)            % fase in radianti
tau_eq = -phiN / wgc         % ritardo equivalente [s]

phiN_deg = phiN*180/pi
tau_ms = tau_eq*1e3
%% Simulink data
clc

A_2 = 1;
f = 55; % Hz
phi_2 = 0;
% d2 = A_2*sin(2*pi*f*t + phi_2);

P_tf = tf(P)
[numP, denP] = tfdata(P_tf, 'v');  % 'v' gives plain vectors
C_tf = tf(C)
[numC, denC] = tfdata(C_tf, 'v');  % 'v' gives plain vectors
H_tf = tf(H)
[numH, denH] = tfdata(H_tf, 'v');  % 'v' gives plain vectors
Ff_tf = tf(Ff)
[numFf, denFf] = tfdata(Ff_tf, 'v');  % 'v' gives plain vectors
Fr_tf = tf(Fr)
[numFr, denFr] = tfdata(Fr_tf, 'v');  % 'v' gives plain vectors

% Saturation on controller
V_MAX = 10; % 

wyr = minreal((P*Ff+P*C*Fr)/(1+P*C*H));
if ~isstable(wyr)
    fprintf("WARNING: SYSTEM UNSTABLE, NOT SAVING")
    step(wyr)
else
    save("simulink_data", "numP", "denP", "numC", "denC", "numH", "denH", "numFf", "denFf", "numFr", "denFr", "A_2", "f", "phi_2", "V_MAX")
    fprintf("saved!\n")
end

%% Different inputs tests
% Impulse
figure(1), impulse(wyr); grid on

% ramp
figure(2), step(wyr); grid on

% step(wyr / s);  % Ramp response
t = 0:1e-3:10;  % the time vector
input = t;   % assuming ramp with slope=1
figure(3)
y = lsim(wyr, input, t);
plot(y, t)
grid on
title("Impulse response")



