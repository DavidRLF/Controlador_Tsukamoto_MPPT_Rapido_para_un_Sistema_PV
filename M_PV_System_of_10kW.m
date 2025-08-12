clc; clear all; close all;

%% ------------------------------------------------------------------------
%  GLOBAL PLOT COSMETICS 
% -------------------------------------------------------------------------
set(groot,'defaultFigureColor','w');          
set(groot,'defaultAxesFontSize',12);          
set(groot,'defaultLineLineWidth',1.4);       

%% ------------------------------------------------------------------------
%  PV ARRAY: Canadian Solar CS6X-280P
% -------------------------------------------------------------------------
No_Cadenas_Paralelo = 2;   % Number of parallel strings
No_Modulos_PorCadena = 18; % Modules per series string
Irr_in=0;                  % Initial irradiance placeholder

%% ------------------------------------------------------------------------
%  BOOST CONVERTER (DC-DC) SIZING & OPERATING POINT
%  Given operating point at ~25°C: Pmp = 10.07 kW, Vmp = 640.8 V
%  vo : PV-side operating voltage (approx Vmp, used as input to boost)
%  vdc: desired DC bus voltage (boost output)
%  Fs : switching frequency
%  Po : PV input power at MPP
%  The following computes duty, load seen by converter, L, output ripple, etc.
% -------------------------------------------------------------------------
vo=640.8;       % PV-side operating voltage ~Vmp [V]
vdc=1000;       % Desired DC bus voltage [V]
Fs=5e3;         % Switching frequency [Hz]
Po=10.07e3;     % PV input power at MPP [W]
io=Po/vo;       % Input current at MPP [A]
iL=io;          % Inductor current (steady-state approx) [A]
d=(vdc-vo)/vdc; % Ideal boost duty [-]
Ro = vo/io;     % MPP equivalent resistance [Ohm]
Iload=io*(1-d); % DC bus current [A]
Rload=vdc/(Iload); % DC bus load [Ohm]
Lmin =(d*(1-d)^2*Rload)/(2*Fs);  % Minimum inductance [H]
L=68.2729*Lmin                  % Selected inductance [H]
ILmax=vo/((1-d)^2*Rload)+(vo*d)/(2*L*Fs) % Inductor current (max) [A]
ILmin=vo/((1-d)^2*Rload)-(vo*d)/(2*L*Fs) % Inductor current (min) [A]
Cload=470e-6;                   % DC bus capacitor [F]
dVdc=(d*vdc)/(Rload*Cload*Fs)   % DC bus ripple [V]
Co=1e-6;                        % Input capacitor [F]
b=16;                           % DPWM resolution [bits]
Dmax=(2^b-1);                   % Max DPWM count
Kdpwm = 1/Dmax;                 % DPWM gain [duty/count]
Kh1=1/io;                       % Sensor gain  

% Small-signal model: X=[iL; vdc], U=[d], Y=[iL]
A=[0 (d-1)/L
   1*(1-d)/Cload -1/(Rload*Cload)];  % State matrix
B=[vdc/L
   -iL/(Cload)];                     % Input matrix
B=B*Kdpwm;                           % Include DPWM gain
C=[1 0];                             % Output iL
D=[0];                               % Direct term
X = {'iL' 'vdc'};    % State names
U = {'d'};           % Input name
Y = {'iL'};          % Output name
sys_cont = ss(A,B,C,D,'statename',X,'inputname',U,'outputname',Y)
tf(sys_cont)

% GsiLd: TF duty->iL (continuous)
GsiLd=ans(1,1)

% --- FIGURE 1: Open-loop current response to a duty step ---
figure(1)
step(((1-d)*Dmax)/2*GsiLd);     % Open-loop response of iL to a duty step
grid on; box on;
title('Open-Loop Response G_{i_L,d} (Duty \rightarrow Inductor Current)');
xlabel('Time [s]');
ylabel('Inductor current i_L [A]');
legend(['Step (scaled, DPWM=',num2str(Dmax),' counts, d=',num2str(d,'%.3f'),')'], ...
       'Location','best');
hold on;

%% ------------------------------------------------------------------------
%  MAMDANI FUZZY CONTROLLER (placeholder/timing alignment)
% -------------------------------------------------------------------------
Ts1=1/Fs;   % Sampling time [s]
Ts2 = Ts1;

%% ------------------------------------------------------------------------
%  TSUKAMOTO FUZZY CONTROLLER (CoDif) FOR BOOST
%  - Two-rule Tsukamoto with monotone triangular consequents
%  - Builds control curve y(x) by weighted average of monotone outputs
% -------------------------------------------------------------------------
paso = 0.001;
Univ_discurso=0.03;     % Duty variation half-range [p.u.]
x=-1:paso:1;            % Error universe Ep [p.u.]
y=-1:paso:1;            % Consequent universe Vm (?d) [p.u.]
a= 1;                   % Sigmoid slope parameter

% Input MFs (Ep): N and P (sigmoids)
mA1 = sigmf (x,[-a 0]); % N
mA2 = sigmf (x,[a 0]);  % P


fig3 = figure(3); clf(fig3);
set(fig3,'Units','pixels','Position',[100 100 1200 700]); % Larger canvas


cN = [0 94 184]/255;   % blue
cP = [220 50 32]/255;  % red

subplot(3,1,1);
plot(x,mA1,'-','LineWidth',2.2,'Color',cN,'DisplayName','N (Negative)'); hold on;
plot(x,mA2,'-','LineWidth',2.2,'Color',cP,'DisplayName','P (Positive)');
grid on; box on;
title('Input Membership Functions (Error Ep)');
xlabel('Error Ep [p.u.]'); ylabel('Membership \mu(Ep) [-]');
xlim([min(x) max(x)]); ylim([0 1.02]);
legend('Location','northwest');


if exist('xline','file') && exist('yline','file')
    yline(0.5,':','\mu = 0.5','LabelHorizontalAlignment','left');
    xline(0,':','Ep = 0','LabelVerticalAlignment','bottom');
else
    plot([min(x) max(x)],[0.5 0.5],'k:'); text(min(x)+0.02,0.53,'\mu = 0.5');
    plot([0 0],[0 1],'k:');           text(0.02,0.95,'Ep = 0');
end

% Output MFs (Vm=?d): N and P (monotone triangles)
mB1 = trimf (y,[-1 1 1 1]); % N (open to the right)
mB2 = trimf (y,[-1 -1 1]);  % P (open to the left)
subplot(3,1,2);
plot(y,mB1,'-','LineWidth',2.2,'Color',cN,'DisplayName','N (monotone)'); hold on;
plot(y,mB2,'-','LineWidth',2.2,'Color',cP,'DisplayName','P (monotone)');
grid on; box on;
title('Output Membership Functions (Command \Delta d)');
xlabel('Duty variation \Delta d [p.u.]'); ylabel('Membership \mu(Vm) [-]');
xlim([min(y) max(y)]); ylim([0 1.02]);
legend('Location','northwest');

% Tsukamoto inference curve (two-rule weighted average)
imax = length(x);
for i=1:imax
    w = [mA1(i) mA2(i)];        % Firing strengths
    z1 = -Univ_discurso*(mA1(i)-0.5); % Rule 1 crisp output [p.u.]
    z2 =  Univ_discurso*(mA2(i)-0.5); % Rule 2 crisp output [p.u.]
    y(i) = (w(1)*z1 + w(2)*z2)/(w(1)+w(2)); % Weighted average
end

subplot(3,1,3);
plot(x,y,'-','LineWidth',2.2,'DisplayName','Control curve'); hold on;
Nmk = max(1,round(numel(x)/40));
plot(x(1:Nmk:end), y(1:Nmk:end), 'o', 'MarkerSize',3, 'HandleVisibility','off');
grid on; box on;
title('Control Curve (Tsukamoto)');
xlabel('Error Ep [p.u.]'); ylabel('Command \Delta d [p.u.]');
xlim([min(x) max(x)]);
ylim([-Univ_discurso/2, Univ_discurso/2]);
yticks(linspace(-Univ_discurso/2,Univ_discurso/2,7));
if exist('xline','file') && exist('yline','file')
    yline(0,':'); xline(0,':');
else
    plot([min(x) max(x)],[0 0],'k:');
    plot([0 0],[-Univ_discurso/2, Univ_discurso/2],'k:');
end

ttlStr = 'Tsukamoto Fuzzy Controller: MFs and Control Curve';
if exist('sgtitle','file')
    sgtitle(ttlStr);
elseif exist('suptitle','file')
    suptitle(ttlStr);
else
    annotation('textbox',[0 0.96 1 0.04], 'String', ttlStr, ...
        'EdgeColor','none', 'HorizontalAlignment','center', ...
        'FontWeight','bold', 'FontSize', 14);
end


Ts3=Ts1; % Sampling time [s]

%% ------------------------------------------------------------------------
%  STATE-FEEDBACK CONTROL (Discrete) WITH INTEGRAL AUGMENTATION
% -------------------------------------------------------------------------
[Ad, Bd, Cd, Dd] = c2dm(A, B, C, D, Ts1,'ZOH');  % Discrete model
X = {'iL' 'vdc'};                                 % State names
U = {'d'};                                        % Input name
Y = {'iL'};                                       % Output name
sys_disc = ss(Ad,Bd,Cd,Dd,Ts1,'statename',X,'inputname',U,'outputname',Y)
tf(sys_disc);
GziLd=ans(1,1)               % Discrete TF duty->iL

% --- FIGURE 4: Open-loop discrete response (duty->iL) ---
figure(4);
step(((1-d)*Dmax)/2*GziLd);  % Open-loop discrete step
grid on; box on;
title('Open-Loop Discrete Response G_{i_L,d}');
xlabel('Time [s]');
ylabel('Inductor current i_L [A]');
legend(['Discrete step (ZOH, T_s=',num2str(Ts1,'%.6g'),' s)'], 'Location','best');
hold on;

poly(eig(sys_disc))          % Char. polynomial (open-loop discrete)
roots([ans])                 % Open-loop poles (discrete)

% Augmented system with integral action
Aa = [Ad Bd;zeros(1,length(Ad)) 0];
Ba = [zeros(length(Bd),1);1];
Ca = [Cd 0];
Da = [Dd 0];
sys_a =ss(Aa,Ba,Ca,0,Ts1)

% --- FIGURE 4 (cont.): Augmented system step (scaled) ---
step((1-d)*Dmax/2*tf(sys_a))
grid on; box on;
title('Augmented System Step (with Integral State)');
xlabel('Time [s]'); ylabel('Output (scaled units)');
legend('Augmented open-loop step','Location','best');
hold off;

% Controllability, characteristic polynomial
Coo = ctrb(sys_a)
rank(Coo) % Should equal size of Aa
E_Ca=poly(eig(sys_a))
roots(E_Ca)

% Desired poles (selected set)
Ed = [0.975;0.99295;0.99];
Ps = poly(Ed)

% State-feedback gain (acker) and recovery of K1/K2
k = acker(Aa,Ba,Ed)
Km = [Ad-eye(length(Ad)) Bd;
      Cd*Ad              Cd*Bd];
In = [zeros(1,length(k)-1) 1];
K = (k+In)/Km;
K1 = K(1:end-1)
K2 = K(end)

% Closed-loop augmented system
Af= [Ad-Bd*K1         Bd*K2;
     -Cd*Ad+Cd*Bd*K1  1-Cd*Bd*K2];
Bf = [zeros(length(Bd),1);1];
Cf = Ca;
Cf = [Ca;zeros(0,length(Ca))] 
slc= ss(Af,Bf,Cf,0,Ts1)

eig(Af)   % Closed-loop poles

% --- FIGURE 5: Closed-loop step response (state feedback + integral) ---
figure(5)
step(io*slc);
grid on; box on;
title('Closed-Loop Step: State-Feedback with Integral Action');
xlabel('Time [s]');
ylabel('Inductor current i_L [A]');
legend(['Closed-loop step (scaled by i_o=',num2str(io,'%.2f'),' A)'], ...
       'Location','best');

%% ------------------------------------------------------------------------
%  GRID-TIED INVERTER MODELING (dq control loops, PWM delay, filters)
% -------------------------------------------------------------------------
Fc=2e3;               % Inverter PWM carrier frequency [Hz]
fs=60;                % Grid frequency [Hz]
ws=2*pi*fs;           % Grid angular frequency [rad/s]
vdc=1000;             % DC bus voltage [V]
S=10e3*3/2;           % Apparent power base [VA]
VL=(381.051/2);       % Line-to-line base voltage [Vrms]
Vf=VL/sqrt(3);        % Phase base voltage [Vrms]
vdc_ref=vdc;          % DC bus reference [V]
Lbase=vdc^2/(ws*S);   % Base inductance [H]

% L filter and internal R (selected values)
L_pu=0.3;
Linv=L_pu*Lbase;
Linv=0.0133;                        % Inductance [H]
Rinv=(Linv*377/fs/2);               % Internal resistance estimate [Ohm]
Rinv=0.0417;                        % Final R [Ohm]
Cinv=470e-6;                        % Inverter input capacitor [F]
Ts4=1/(4*Fc);                       % Inner-loop sampling [s]
Td=Ts4/2;                           % PWM update delay [s]
Ts_inv=1/(100*Fc);                  % PLL sampling [s]
Vp=1;                               % PWM carrier peak [V]

% PI design (id/iq loops, DC bus, PLL)
Porcentaje1=1;  Mp1=Porcentaje1/100;
Zeta1=sqrt(log(Mp1)^2/(log(Mp1)^2+pi^2));
tset1=10e-3;    Porcentaje=2; E1=Porcentaje/100;
wn1=-log(E1)/(Zeta1*tset1);
Kp2=2*Linv*Vp*Zeta1*wn1-Rinv*Vp;
Ki2=Linv*Vp*wn1^2;

Tlpf=0.02; Alpha=2;
Kp3=(Cinv)/(Alpha*Tlpf);
Ki3=(Cinv)/(Alpha^3*Tlpf^2);

Porcentaje=1;  Mp2=Porcentaje/100;
Zeta2=sqrt(log(Mp2)^2/(log(Mp2)^2+pi^2));
n=2; tset2=tset1*n; Porcentaje=4; E2=Porcentaje/100;
wn2=-log(E2)/(Zeta2*tset2);
Kp4=2*Cinv*Zeta2*wn2;
Ki4=Cinv*wn2^2;

Zeta3=0.707; wn3=(120*pi)/2;
Kp5=(2*Zeta3*wn3);
Ki5=wn3^2;

% Plants (continuous with PWM delay) and discrete (ZOH)
Gidedp_s=(tf([1],[Linv Rinv],'InputDelay',Td));
Giqedp_s=(tf([1],[Linv Rinv],'InputDelay',Td));
Gidedp_z1=c2d(Gidedp_s,Ts4,'ZOH');
Giqedp_z1=c2d(Giqedp_s,Ts4,'ZOH');

m=0.5; % Fractional delay factor

% Discrete plant coefficients for Simulink blocks (id)
ad2=0; ad1=1-exp(-(Rinv*Ts4*m)/Linv);
ad0=(exp(-(Rinv*Ts4*m)/Linv)-exp(-(Rinv*Ts4)/Linv));
bd2=Rinv; bd1=-Rinv*exp(-(Rinv*Ts4)/Linv); bd0=0;
Gidedp_z2=zpk(tf([ad2 ad1 ad0],[bd2 bd1 bd0],Ts4))
% Should match:
Gidedp_z1

% Discrete plant coefficients for Simulink blocks (iq)
aq2=0; aq1=1-exp(-(Rinv*Ts4*m)/Linv);
aq0=(exp(-(Rinv*Ts4*m)/Linv)-exp(-(Rinv*Ts4)/Linv));
bq2=Rinv; bq1=-Rinv*exp(-(Rinv*Ts4)/Linv); bq0=0;
Giqedp_z2=zpk(tf([aq2 aq1 aq0],[bq2 bq1 bq0],Ts4))
% Should match:
Giqedp_z1

% LPF for id feedback (continuous and discrete)
Glpf_s=1/Tlpf*(tf([1],[1 1/Tlpf]));
Glpf_z1=c2d(Glpf_s,Ts4,'ZOH');

% Discrete LPF coefficients for Simulink
alpf1=0; alpf0=1-exp(-Ts4/Tlpf);
blpf1=1; blpf0=-exp(-Ts4/Tlpf);
Glpf_z2=zpk(tf([alpf1 alpf0],[blpf1 blpf0],Ts4))
% Should match:
Glpf_z2

% DC bus plant (continuous and discrete)
Gvdcid_s=tf([1],[Cinv 0]);
Gvdcid_z1=c2d(Gvdcid_s,Ts4,'ZOH');

% Discrete DC bus coefficients for Simulink
avdc1=0; avdc0=Ts4/Cinv;
bvdc1=1; bvdc0=-1;
Gvdcid_z2=zpk(tf([avdc1 avdc0],[bvdc1 bvdc0],Ts4))
% Should match:
Gvdcid_z1

% PWM gain and delay (continuous and discrete)
Kpwm_s=1/Vp; Kpwm_z=1/Vp;
Gret_s=tf([1],[Td 1]); Gret_z1=c2d(Gret_s,Ts4,'ZOH');

% Discrete PWM delay coefficients for Simulink
aRet1=0; aRet0=1-1*exp(-Ts4/Td);
bRet1=1; bRet0=-exp(-Ts4/Td);
Gret_z2=zpk(tf([aRet1 aRet0],[bRet1 bRet0],Ts4))
% Should match:
Gret_z2

% PLL plant (continuous and discrete)
GPLL_s=tf([0 1],[1 0]);
GPLL_z1=c2d(GPLL_s,Ts_inv,'ZOH');

% Discrete PLL coefficients for Simulink
aPLL1=0; aPLL0=Ts_inv;
bPLL1=1; bPLL0=-1;
GPLL_z2=zpk(tf([aPLL1 aPLL0],[bPLL1 bPLL0],Ts_inv))
% Should match:
GPLL_z1

% Continuous-to-discrete PI gains
KI2=Ki2*Ts4; KP2=Kp2-KI2/2;
KI3=Ki3*Ts4; KP3=Kp3-KI3/2;
KI4=Ki4*Ts4; KP4=Kp4-KI4/2;
KI5=Ki5*Ts_inv; KP5=Kp5*Ts_inv-KI5/2;

%% ------------------------------------------------------------------------
%  TRANSFORMER (50 kVA exemplar — small building proxy)
% -------------------------------------------------------------------------
Ptrafo   = 100e3;     % Nominal power [W]
Ftrafo   = fs;        % Operating frequency [Hz]
VLLsec   = 33000;     % Line-to-line voltage (grid side) [Vrms]
VLLprim  = VL;        % Line-to-line voltage (inverter side) [Vrms]

%% ------------------------------------------------------------------------
%  GRID PARAMETERS
% -------------------------------------------------------------------------
VLngrid = 33000/sqrt(3); % Line-to-neutral grid voltage [Vrms]
Fred    = fs;            % Grid frequency [Hz]
In = 1e6/(1.73*(Vf));   % Nominal current (short-circuit calc) [A]
Icc_calc = In/(5/100);  % Short-circuit current estimate (5% impedance) [A]
