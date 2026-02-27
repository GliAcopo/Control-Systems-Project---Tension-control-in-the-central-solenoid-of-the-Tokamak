# Control Systems Project

Author: Jacopo Rizzuto  
Date: February 12, 2026  
Presented: February 25, 2026

---

## The problem

In the Tokamak:

- Right after the plasma cools, the *plasma current* $I_p$ collapses (1-10 ms).
- Because $I_p$ is dropping fast, the voltage raises fast **($V \approx - L \frac{d I_p}{d t}$)**.

That strong electric field is one of the main reasons **runaway electrons** can be generated/avalanche during CQ.

- Magnetic control **is not fast enough to prevent RE creation** due to:
  - Coils of the magnets having huge inductance.
  - Power supplies having phisical limits on voltage/current erogation.

---

## The solution

Control the position of the RE current by:

1. Control the current in position ($I$ sent to the position coils)
2. Safely ramp down the $I_{RE}$ current by waiting for it to be absorbed by the central current magnet.

Like shown in the graph.

Figure link: [Progetto controlli-3 equivalent](./Figures/1.jpg)

The system: controls **the electric field that drives toroidal current** by acting on the **central solenoid (also called transformer)**.

- If $I_p$ is **below** the reference $r$: it commands the solenoid to produce **more loop voltage** (toroidal $E$), so the plasma current *tends* to **increase**.
- If $I_p$ is **above** $r$: it commands **less loop voltage**, so the current tends to **decrease**.

Control procedure:

- We compare the **plasma current reference** $r$ with the measured **plasma current** $I_p$.
- The **Controller** will output the correction signal $u_1$.
- That correction is **added to a preprogrammed feedforward waveform** $I_{T\,prep}$.
- The AL-T will handle controlling of the solenoid in order to get the current to $I_{Tref}$.

---

## The model

The model is as follows, the relationship between central transformer and the plasma can be seen as two circuits coupled with a solenoid.

Figure link: [Tokamak diagram](./Figures/TOKAMAK_DIAGRAM.svg)

---

## The model - Concatenated flux

(References in original slides: link1, link2, pdf page 14)

*$M$ by definition is how much concatenated flux is produced in $P$ by the current flowing in $T$.*  
*By convention we choose $M$ as positive (the plasma/transformer coupling is considered positive).*  
*We also assume that $M_{T,P}=M_{P,T}$, that is not always the case and we assume that in the tokamak it can happen.*

Being these two circuits coupled by a "trasformer-like" link we can calculate the **concatenated flux** between them:

$$
\lambda = N\Phi
$$

$$
N := \text{number of loops of the bobine}
$$

$$
\Phi := \text{magnetic flux that each loop percieves}
$$

Fortunately we can ignore $\Phi$ and write the formula in terms of total inductance of the circuits:

- $i_t, i_p$: plasma and transformer current
- $L_T, L_P$: "autoinductance" of circuit $n$, $L_n = \tfrac{\partial \lambda_n}{\partial i_n}$
- $M_{T,P}$: mutual inductance of the circuits, $M_{T,P}=\tfrac{\partial\lambda_1}{\partial i_2}=\tfrac{\partial\lambda_2}{\partial i_1}$
- In general $M = k\sqrt{L_TL_P}$ with coupling coefficient $k\in[0,1]$

$$
\lambda_T = L_T i_T + M i_P
$$

$$
\lambda_P = M i_T + L_P i_P
$$

$$
\begin{bmatrix}
\lambda_T\\
\lambda_P
\end{bmatrix}
=
\begin{bmatrix}
L_T & M\\
M & L_P
\end{bmatrix}
\begin{bmatrix}
i_T\\
i_P
\end{bmatrix}
$$

---

## The model - From concatenated flux to voltages

We know from Kirchhoff:

$$
v = Ri + \frac{d\lambda}{dt}
$$

So:

$$
v_T = R_T i_T + \frac{d}{dt}(L_T i_T + M i_P)
$$

$$
v_P = R_P i_P + \frac{d}{dt}(M i_T + L_P i_P)
$$

Note: $v_p$ is a purely induced f.e.m. by the central transformer circuit $T$.

The nonlinear system (matricial form):

$$
\begin{bmatrix}v_T\\v_P\end{bmatrix}
=
\begin{bmatrix}R_T(K^\circ) & R_P(K^\circ)\end{bmatrix}
\begin{bmatrix}i_T\\i_P\end{bmatrix}
+
\begin{bmatrix}
\frac{d}{dt}L_T(i_T) & M(x,y,z)\\
M(x,y,z) & \frac{d}{dt}L_P(i_P)
\end{bmatrix}
\begin{bmatrix}\frac{d}{dt}i_T\\\frac{d}{dt}i_P\end{bmatrix}
$$

$$
L(i):=\text{inductance dependent on current }i
$$

$$
M(x,y,z,L):=\text{mutual inductance dependent on position and circuit inductance}
$$

$$
R(K^\circ):=\text{Resistance dependent on plasma temperature}
$$

---

## Linearization - Assumptions

We place the equilibrium point in a plateau state $\overline{x}_0$ of the plasma:

1. We consider resistence $R$ to be constant. (constant temperature $K^\circ = 150$ millions)
2. We consider inductance $L$ dependent on $i$, so we choose the equilibrium point as:

$$
\begin{bmatrix}I_T\\I_P\end{bmatrix}=\overline{0}
$$

3. We consider $M$ to be constant (plasma position does not change).

From nonlinear to "linear" system:

$$
\begin{bmatrix}v_T\\v_P\end{bmatrix}
=
\begin{bmatrix}R_T & R_P\end{bmatrix}
\begin{bmatrix}i_T\\i_P\end{bmatrix}
+
\begin{bmatrix}L_T & M\\M & L_P\end{bmatrix}
\begin{bmatrix}\frac{d}{dt}i_T\\\frac{d}{dt}i_P\end{bmatrix}
$$

---

## Linearization - Obtaining $\dot{x}$

$$
\begin{bmatrix}L_T & M\\M & L_P\end{bmatrix}
\begin{bmatrix}\dot{i}_T\\\dot{i}_P\end{bmatrix}
=
\begin{bmatrix}v_T\\v_P\end{bmatrix}
-
\begin{bmatrix}R_T & R_P\end{bmatrix}
\begin{bmatrix}i_T\\i_P\end{bmatrix}
$$

$$
\frac{d}{dt}\begin{bmatrix}i_T\\i_P\end{bmatrix}
=
\begin{bmatrix}L_T & M\\M & L_P\end{bmatrix}^{-1}
\left(
\begin{bmatrix}v_T\\v_P\end{bmatrix}
-
\begin{bmatrix}R_T I_T\\R_P I_P\end{bmatrix}
\right)
$$

So:

$$
\dot{x} =
\begin{bmatrix}L_T & M\\M & L_P\end{bmatrix}^{-1}
\left(
\begin{bmatrix}v_T\\v_P\end{bmatrix}
-
\begin{bmatrix}R_T I_T\\R_P I_P\end{bmatrix}
\right)
$$

---

## Linearization - Variables

- State: $x(t)=\begin{bmatrix}I_T\\I_P\end{bmatrix}$
- Reference: $r(t)=I_{ref}(t)$
- Output: $y(t)=I_P(t)$
- Control: $u(t)=v_T(t)$

Since the system is SISO, we can control only one of the two tensions, $v_T$ and not $v_P$. Since $v_P$ is induced and dependent on $v_T$, we treat plasma-side tension as exogenous disturbance $e_{v_p}$:

$$
\dot{x} =
\begin{bmatrix}L_T & M\\M & L_P\end{bmatrix}^{-1}
\left(
\begin{bmatrix}v_T\\0\end{bmatrix}
-
\begin{bmatrix}R_T I_T\\R_P I_P + e_{v_p}\end{bmatrix}
\right)
$$

---

## Linearization - State space representation

Define

$$
K = \begin{bmatrix}L_T & M\\M & L_P\end{bmatrix}^{-1}
= \frac{1}{L_TL_P-M^2}
\begin{bmatrix}L_P & -M\\-M & L_T\end{bmatrix}
$$

with

$$
K_\Delta := L_TL_P - M^2 > 0
$$

Then

$$
\dot{x} = A x + B u + E e_{v_p}, \qquad y = Cx + D
$$

with

$$
A=\frac{1}{K_\Delta}
\begin{bmatrix}
-L_P R_T & MR_T\\
MR_P & -L_TR_P
\end{bmatrix},
\quad
B=\frac{1}{K_\Delta}\begin{bmatrix}L_P\\-M\end{bmatrix},
\quad
E=\frac{1}{K_\Delta}\begin{bmatrix}M\\-L_T\end{bmatrix}
$$

$$
C=\begin{bmatrix}0 & 1\end{bmatrix},\quad D=\begin{bmatrix}0\end{bmatrix},\quad y=I_P
$$

---

## Process $P(s)$ expression

Assuming exogenous disturbances $e=0$:

$$
P(s)=\frac{y(s)}{u(s)}=\frac{I_P}{v_T}
$$

From

$$
v_T = R_T I_T + L_T\dot{I}_T + M\dot{I}_P
$$

and, with $v_P=0$,

$$
R_P I_P + L_P\dot{I}_P = -M\dot{I}_T
$$

Laplace-domain derivation gives:

$$
P(s)=\frac{I_P}{v_T}=
\frac{-Ms}{(R_T+L_T s)(R_P+L_P s)-M^2 s^2}
$$

---

## Constant values assumptions

- Coupling coefficient $k=0.8$.
- Inductance $L_P\approx1.9\,\mu H$; $L_T\approx80\,mH$.
- Resistance $R_P\approx4\,\mu\Omega$; $R_T\approx1\,m\Omega$.

---

## System requirements

1. $|S|\% < 20\%$
2. $T_a < 0.001\,s$
3. Required by slides: $e_\infty=0$ for $r=\frac{1}{s^2}$
4. $d_1$: impulsive disturbance on actuator (T-coil):

$$
u(t)=u_c(t)+d_1(t),\quad d_1(t)=A_1\delta(t-t_0),\quad \int d_1(t)dt=A_1
$$

5. $d_2$: disturbance ripple in sensor reading:

$$
y_m(t)=y(t)+d_2(t),\quad d_2(t)=A_2\sin(2\pi f t+\phi_2),\quad f\in\{50,60\}\,Hz
$$

6. Robustness $R_p\pm30\%$ (optimally $\pm50\%$) for changes in $R_T$ and plasma displacement $M$.

---

## Controller Design Process

### Block algebra to find $W_{yr}$

$$
y = P\left(d_1 + rF_f + C(rF_r - H(y+d_2))\right)
$$

$$
y(1+PCH)=Pd_1 + PF_f r + PCF_r r - PCHd_2
$$

With $d_1=d_2=0$:

$$
W_{yr}=\frac{y}{r}=\frac{PF_f+PCF_r}{1+PCH}
$$

### Final Value Theorem

Our system does not have any integrators. Because process has a zero in $s=0$, overall loop type is $\rho_c-1$.

Condition derived in slides:

$$
\rho_c-\rho_r\ge1 \Rightarrow \rho_c\ge\rho_r+1
$$

---

## Transient state analysis

### ANALYSYS (1)

First naive attempt: cancel all poles of process.

```matlab
P = zpk((-M*s)/(((R_T+L_T*s)*(R_P+L_P*s)-M^2*s^2)))
cancellationC = (((R_T+L_T*s)*(R_P+L_P*s)-M^2*s^2));
```

(1.A) No added zero: system at stability limit with self-sustained oscillations.

Figure links:
- [Rlocus - Case (1.A.I)](./Figures/Rlocus_-_Case_(1.A.I).pdf)
- [Rlocus - Case (1.A.II)](./Figures/Rlocus_-_Case_(1.A.II).pdf)
- [Step Response - Case (1.A.II)](./Figures/Step_Response_-_Case_(1.A.II).pdf)

(1.B) Add one zero to guarantee stability.

```matlab
zeroesC = (s + 1);
polesC = 1;
k = -64.371;
alpha = 100;
tau = 0.5;
Rr = (1+tau*alpha*s)/(1+tau*s);
C = minreal(Rr * cancellationC * 1/polesC * zeroesC * (k/s^rho))
```

Performances (reported):

- RiseTime: 2.4878
- TransientTime: 4.4700
- SettlingTime: 4.4700
- Overshoot: 1.6271

Robustness note from slides: changing $R_T$ gives very small changes (order $10^{-4}$ in Bode), changing $M$ gives much larger changes (order $10^2$).

### ANALYSYS (2)

More robust approach: do **not** cancel all poles.

(2.A) Positive gain, even with three zeroes system cannot be made stable.

```matlab
k = 1;
zeroesC = (s+1)^3;
polesC = 1;
C = minreal(1/polesC * zeroesC * (k/s^rho));
H = 1; Ff = 0; Fr = 1;
```

(2.B) Negative gain, two zeroes are enough.

```matlab
k = -1;
zeroesC = (s+1)^2;
polesC = 1;
C = minreal(1/polesC * zeroesC * (k/s^rho));
H = 1; Ff = 0; Fr = 1;
```

(2.B.I) Tuning only gain gives very good performances ($k=-108$):

- RiseTime: 0.0038
- TransientTime: 0.0062
- SettlingTime: 0.0062
- Overshoot: 0.5021

System considered stable (all poles in left half-plane).

---

## Disturbance handling

From

$$
y = \frac{Pd_1 + PF_f r + PCF_r r - PCHd_2}{1+PCH}
$$

and

$$
W_{y,r}=\frac{PF_f+PCF_r}{1+PCH},\quad
W_{y,d_1}=\frac{P}{1+PCH},\quad
W_{y,d_2}=-\frac{PCH}{1+PCH}
$$

The disturbance from $d_2$ had strong effect in intermediate designs.

Used in simulation:

```matlab
Wer = minreal((Fr - H*P*Ff)/(1 + P*C*H));
Wed1 = minreal(-(H*P)/(1 + P*C*H));
Wed2 = minreal(-H/(1 + P*C*H));

dt = 1e-4;
t = 0:dt:0.5;
A_2 = 1; f = 60; phi_2 = 0;
d2 = A_2*sin(2*pi*f*t + phi_2);
```

Impulse approximation:

```matlab
A1 = 0.05;
d1 = zeros(size(t));
idx = find(t >= t0, 1, 'first');
d1(idx) = A1/dt;
```

---

## ANALYSYS (2.B.II) - notch filter

A notch filter was added:

```matlab
w1 = 2*pi*(f-5);
w2 = 2*pi*(f+5);
n = 2;
[b,a] = butter(n, [w1, w2], 'stop', 's');
H = tf(b, a)
```

Performance example reported:

- RiseTime: 0.0144
- SettlingTime: 0.0508
- Overshoot: 2.1594

But crossover frequency was in the disturbed 50-60 Hz region, so design was remade.

Figure links:
- [Step Response - Case (2.B.II)](./Figures/12-Step_Response_-_Case_(2.B.II).pdf)
- [Bode L - Case (2.B.II)](./Figures/13-Bode_L_-_Case_(2.B.II).pdf)
- [Rlocus L - Case (2.B.II)](./Figures/14-Rlocus_L_-_Case_(2.B.II).pdf)
- [Nyquist - Case (2.B.II)](./Figures/15-Logarithmic_Nyquist_of_L_-_Case_(2.B.II).pdf)
- [Rlocus Zoom - Case (2.B.II)](./Figures/16-Rlocus_L_-_Case_(2.B.II)_ZOOM.pdf)

---

## ANALYSYS (2.B.III)

Move crossover frequency much lower than disturbance band.

For $k=-1$:

- RiseTime: 0.2290
- SettlingTime: 2.8346
- Overshoot: 41.1302
- $WcpHz = 0.7364$

This reduced disturbance overlap but badly degraded performance.

### Add aggressive feedforward

```matlab
fc_ff = 1000;
wn = 2*pi*fc_ff;
zeta = 0.707;
Ff = minreal((1/P) * (wn^2 / (s^2 + 2*zeta*wn*s + wn^2)));
```

Performance with aggressive $F_f$:

- RiseTime: 6.8262e-04
- SettlingTime: 0.0019
- Overshoot: 4.4729

Note from original: loop became much slower than $F_f$ branch (about $149000$ times).

With notch active, performance dictated by $F_f$ remained essentially unchanged in slides.

Figure links:
- [Step Response - Case (2.B.III-gain_mod)](./Figures/12-Step_Response_-_Case_(2.B.III-gain_mod).pdf)
- [Step Response - Case (2.B.III-gain_mod-FfP2)](./Figures/12-Step_Response_-_Case_(2.B.III-gain_mod-FfP2).pdf)
- [Step Response - Case (2.B.III-gain_mod-butter-FfP2)](./Figures/12-Step_Response_-_Case_(2.B.III-gain_mod-butter-FfP2).pdf)
- [Bode L - Case (2.B.III-gain_mod)](./Figures/13-Bode_L_-_Case_(2.B.III-gain_mod).pdf)
- [Bode L - Case (2.B.III-gain_mod-FfP2)](./Figures/13-Bode_L_-_Case_(2.B.III-gain_mod-FfP2).pdf)
- [Bode L - Case (2.B.III-gain_mod-butter-FfP2)](./Figures/13-Bode_L_-_Case_(2.B.III-gain_mod-butter-FfP2).pdf)
- [Rlocus L - Case (2.B.III-gain_mod)](./Figures/14-Rlocus_L_-_Case_(2.B.III-gain_mod).pdf)
- [Rlocus L - Case (2.B.III-gain_mod-FfP2)](./Figures/14-Rlocus_L_-_Case_(2.B.III-gain_mod-FfP2).pdf)
- [Rlocus L - Case (2.B.III-gain_mod-butter-FfP2)](./Figures/14-Rlocus_L_-_Case_(2.B.III-gain_mod-butter-FfP2).pdf)
- [Nyquist - Case (2.B.III-gain_mod)](./Figures/15-Logarithmic_Nyquist_of_L_-_Case_(2.B.III-gain_mod).pdf)
- [Nyquist - Case (2.B.III-gain_mod-FfP2)](./Figures/15-Logarithmic_Nyquist_of_L_-_Case_(2.B.III-gain_mod-FfP2).pdf)
- [Nyquist - Case (2.B.III-gain_mod-butter-FfP2)](./Figures/15-Logarithmic_Nyquist_of_L_-_Case_(2.B.III-gain_mod-butter-FfP2).pdf)
- [Rlocus Zoom - Case (2.B.III-gain_mod)](./Figures/16-Rlocus_L_-_Case_(2.B.III-gain_mod)_ZOOM.pdf)
- [Rlocus Zoom - Case (2.B.III-gain_mod-FfP2)](./Figures/16-Rlocus_L_-_Case_(2.B.III-gain_mod-FfP2)_ZOOM.pdf)
- [Rlocus Zoom - Case (2.B.III-gain_mod-butter-FfP2)](./Figures/16-Rlocus_L_-_Case_(2.B.III-gain_mod-butter-FfP2)_ZOOM.pdf)
- [Output response to d2 - (2.B.III-gain_mod)](./Figures/(2.B.III-gain_mod)_-_Output_response_to_sensor_ripple_disturbance_d_2.pdf)
- [Total output y(t) - (2.B.III-gain_mod)](./Figures/(2.B.III-gain_mod)_-_Total_output_y(t)_under_r,_d_1,_d_2.pdf)
- [Output response to d2 - (2.B.III-gain_mod-FfP2)](./Figures/(2.B.III-gain_mod-FfP2)_-_Output_response_to_sensor_ripple_disturbance_d_2.pdf)
- [Total output y(t) - (2.B.III-gain_mod-FfP2)](./Figures/(2.B.III-gain_mod-FfP2)_-_Total_output_y(t)_under_r,_d_1,_d_2.pdf)
- [Output response to d2 - (2.B.III-gain_mod-butter-FfP2)](./Figures/(2.B.III-gain_mod-butter-FfP2)_-_Output_response_to_sensor_ripple_disturbance_d_2.pdf)
- [Total output y(t) - (2.B.III-gain_mod-butter-FfP2)](./Figures/(2.B.III-gain_mod-butter-FfP2)_-_Total_output_y(t)_under_r,_d_1,_d_2.pdf)

---

## System stability

The system is confirmed stable in final analysis branch (2.B.III).

Figure links:
- [NyquistLogFf](./Figures/NyquistLogFf.pdf)
- [Margin FfP2](./Figures/Margin_FfP2.pdf)

---

## Parameter robustness analysis (2.B.III)

Observed:

- Stable for about $\pm50\%$ on $R_T$.
- Stops around $30\%$ robustness regarding $M$.

Figure links:
- [Robustness step response (R_T)](./Figures/1-R(2.B.III)-Robustness_R_T_step_response_of_(2.B.III)_system.pdf)
- [Robustness step response (M, R_T)](./Figures/1-R(2.B.III)-Robustness_M-R_T_step_response_of_(2.B.III)_system.pdf)

---

## Maximum acceptable delay

Loop delay margin:

$$
T_{max} \approx \frac{PM\,[rad]}{\omega_{gc}\,[rad/s]}
$$

With reported values:

$$
PM=42.9^\circ,\quad \omega_{gc}=4.63\,rad/s\Rightarrow T_{max}\approx0.162\,s
$$

Delay from notch at $\omega^*=\omega_{gc}$:

$$
T_{notch,\omega^*}=-\frac{\angle N(j\omega^*)}{\omega^*}
$$

with reported:

$$
\angle N(j\omega^*)=-0.0029\Rightarrow T_{notch,\omega^*}=0.0006297\,s
$$

---

## System limits

1. **Saturation and limits on real actuators**
   - voltage saturation: $|v_T|\le V_{max}$
   - slew-rate: $|\dot I_T|\le\dot I_{T,max}$
2. **Plasma movement changes circuit coupling**
   - $M$ depends on reciprocal position of circuits.
3. **Disturbance on reference would not be very filtered**
   - because of aggressive feedforward and relatively slow loop, fragility to reference disturbance.

---

## System simulation

Since the system is linear, the dynamics are represented entirely by block $P$.

---

## Appunti - The problem / Tokamak phases

In a tokamak disruption there are usually two main fast phases:

1. **Thermal Quench (TQ)**
   - Plasma temperature collapses very quickly.
   - Happens in sub-ms to a few ms.
   - Energy is lost mainly by transport + radiation.
2. **Current Quench (CQ)**
   - Right after plasma cools, plasma current $I_p$ collapses (drops fast).
   - Happens in a few ms to tens of ms.
   - Because $I_p$ drops fast, voltage raises fast **($V\approx-L\frac{dI_p}{dt}$)**.
   - That strong electric field is one of the main reasons runaway electrons can be generated/avalanche during CQ.

---

## Notes

This file is a direct Markdown conversion of `Progetto controlli.tex`, preserving wording and equations from the slides and only adapting formatting/syntax.

Please refer to the original `.tex` file for the original formatting and to the slides for the original figures.

## All Figure File Links

- [(2.B.I)_-_Output_response_to_sensor_ripple_disturbance_d_2.pdf](<./Figures/(2.B.I)_-_Output_response_to_sensor_ripple_disturbance_d_2.pdf>)
- [(2.B.I)_-_Total_output_y(t)_under_r,_d_1,_d_2.pdf](<./Figures/(2.B.I)_-_Total_output_y(t)_under_r,_d_1,_d_2.pdf>)
- [(2.B.II)_-_Output_response_to_sensor_ripple_disturbance_d_2.pdf](<./Figures/(2.B.II)_-_Output_response_to_sensor_ripple_disturbance_d_2.pdf>)
- [(2.B.II)_-_Total_output_y(t)_under_r,_d_1,_d_2.pdf](<./Figures/(2.B.II)_-_Total_output_y(t)_under_r,_d_1,_d_2.pdf>)
- [(2.B.III-gain_mod)_-_Output_response_to_sensor_ripple_disturbance_d_2.pdf](<./Figures/(2.B.III-gain_mod)_-_Output_response_to_sensor_ripple_disturbance_d_2.pdf>)
- [(2.B.III-gain_mod)_-_Total_output_y(t)_under_r,_d_1,_d_2.pdf](<./Figures/(2.B.III-gain_mod)_-_Total_output_y(t)_under_r,_d_1,_d_2.pdf>)
- [(2.B.III-gain_mod-FfP2)_-_Output_response_to_sensor_ripple_disturbance_d_2.pdf](<./Figures/(2.B.III-gain_mod-FfP2)_-_Output_response_to_sensor_ripple_disturbance_d_2.pdf>)
- [(2.B.III-gain_mod-FfP2)_-_Total_output_y(t)_under_r,_d_1,_d_2.pdf](<./Figures/(2.B.III-gain_mod-FfP2)_-_Total_output_y(t)_under_r,_d_1,_d_2.pdf>)
- [(2.B.III-gain_mod-butter-FfP2)_-_Output_response_to_sensor_ripple_disturbance_d_2.pdf](<./Figures/(2.B.III-gain_mod-butter-FfP2)_-_Output_response_to_sensor_ripple_disturbance_d_2.pdf>)
- [(2.B.III-gain_mod-butter-FfP2)_-_Total_output_y(t)_under_r,_d_1,_d_2.pdf](<./Figures/(2.B.III-gain_mod-butter-FfP2)_-_Total_output_y(t)_under_r,_d_1,_d_2.pdf>)
- [1-R(1.B)-Robustness_M-R_T_step_response_of_(1.B)_system.pdf](<./Figures/1-R(1.B)-Robustness_M-R_T_step_response_of_(1.B)_system.pdf>)
- [1-R(1.B)-Robustness_R_T_step_response_of_(1.B)_system.pdf](<./Figures/1-R(1.B)-Robustness_R_T_step_response_of_(1.B)_system.pdf>)
- [1-R(1.B)-ZOOM-Robustness_R_T_step_response_of_(1.B)_system.pdf](<./Figures/1-R(1.B)-ZOOM-Robustness_R_T_step_response_of_(1.B)_system.pdf>)
- [1-R(2.B.III)-Robustness_M-R_T_step_response_of_(2.B.III)_system.pdf](<./Figures/1-R(2.B.III)-Robustness_M-R_T_step_response_of_(2.B.III)_system.pdf>)
- [1-R(2.B.III)-Robustness_R_T_step_response_of_(2.B.III)_system.pdf](<./Figures/1-R(2.B.III)-Robustness_R_T_step_response_of_(2.B.III)_system.pdf>)
- [1.jpg](<./Figures/1.jpg>)
- [10-Rlocus_L_-_Case_(1.B).pdf](<./Figures/10-Rlocus_L_-_Case_(1.B).pdf>)
- [11-Logarithmic_Nyquist_of_L_-_Case_(1.B).pdf](<./Figures/11-Logarithmic_Nyquist_of_L_-_Case_(1.B).pdf>)
- [12-Step_Response_-_Case_(2.A).pdf](<./Figures/12-Step_Response_-_Case_(2.A).pdf>)
- [12-Step_Response_-_Case_(2.B).pdf](<./Figures/12-Step_Response_-_Case_(2.B).pdf>)
- [12-Step_Response_-_Case_(2.B.I).pdf](<./Figures/12-Step_Response_-_Case_(2.B.I).pdf>)
- [12-Step_Response_-_Case_(2.B.II).pdf](<./Figures/12-Step_Response_-_Case_(2.B.II).pdf>)
- [12-Step_Response_-_Case_(2.B.III-gain_mod).pdf](<./Figures/12-Step_Response_-_Case_(2.B.III-gain_mod).pdf>)
- [12-Step_Response_-_Case_(2.B.III-gain_mod-FfP2).pdf](<./Figures/12-Step_Response_-_Case_(2.B.III-gain_mod-FfP2).pdf>)
- [12-Step_Response_-_Case_(2.B.III-gain_mod-butter-FfP2).pdf](<./Figures/12-Step_Response_-_Case_(2.B.III-gain_mod-butter-FfP2).pdf>)
- [13-Bode_L_-_Case_(2.A).pdf](<./Figures/13-Bode_L_-_Case_(2.A).pdf>)
- [13-Bode_L_-_Case_(2.B).pdf](<./Figures/13-Bode_L_-_Case_(2.B).pdf>)
- [13-Bode_L_-_Case_(2.B.I).pdf](<./Figures/13-Bode_L_-_Case_(2.B.I).pdf>)
- [13-Bode_L_-_Case_(2.B.II).pdf](<./Figures/13-Bode_L_-_Case_(2.B.II).pdf>)
- [13-Bode_L_-_Case_(2.B.III-gain_mod).pdf](<./Figures/13-Bode_L_-_Case_(2.B.III-gain_mod).pdf>)
- [13-Bode_L_-_Case_(2.B.III-gain_mod-FfP2).pdf](<./Figures/13-Bode_L_-_Case_(2.B.III-gain_mod-FfP2).pdf>)
- [13-Bode_L_-_Case_(2.B.III-gain_mod-butter-FfP2).pdf](<./Figures/13-Bode_L_-_Case_(2.B.III-gain_mod-butter-FfP2).pdf>)
- [14-Rlocus_L_-_Case_(2.A).pdf](<./Figures/14-Rlocus_L_-_Case_(2.A).pdf>)
- [14-Rlocus_L_-_Case_(2.A)_ZOOM.pdf](<./Figures/14-Rlocus_L_-_Case_(2.A)_ZOOM.pdf>)
- [14-Rlocus_L_-_Case_(2.B).pdf](<./Figures/14-Rlocus_L_-_Case_(2.B).pdf>)
- [14-Rlocus_L_-_Case_(2.B)_ZOOM.pdf](<./Figures/14-Rlocus_L_-_Case_(2.B)_ZOOM.pdf>)
- [14-Rlocus_L_-_Case_(2.B.I).pdf](<./Figures/14-Rlocus_L_-_Case_(2.B.I).pdf>)
- [14-Rlocus_L_-_Case_(2.B.I)_ZOOM.pdf](<./Figures/14-Rlocus_L_-_Case_(2.B.I)_ZOOM.pdf>)
- [14-Rlocus_L_-_Case_(2.B.II).pdf](<./Figures/14-Rlocus_L_-_Case_(2.B.II).pdf>)
- [14-Rlocus_L_-_Case_(2.B.III-gain_mod).pdf](<./Figures/14-Rlocus_L_-_Case_(2.B.III-gain_mod).pdf>)
- [14-Rlocus_L_-_Case_(2.B.III-gain_mod-FfP2).pdf](<./Figures/14-Rlocus_L_-_Case_(2.B.III-gain_mod-FfP2).pdf>)
- [14-Rlocus_L_-_Case_(2.B.III-gain_mod-butter-FfP2).pdf](<./Figures/14-Rlocus_L_-_Case_(2.B.III-gain_mod-butter-FfP2).pdf>)
- [15-Logarithmic_Nyquist_of_L_-_Case_(2.A).pdf](<./Figures/15-Logarithmic_Nyquist_of_L_-_Case_(2.A).pdf>)
- [15-Logarithmic_Nyquist_of_L_-_Case_(2.B).pdf](<./Figures/15-Logarithmic_Nyquist_of_L_-_Case_(2.B).pdf>)
- [15-Logarithmic_Nyquist_of_L_-_Case_(2.B.I).pdf](<./Figures/15-Logarithmic_Nyquist_of_L_-_Case_(2.B.I).pdf>)
- [15-Logarithmic_Nyquist_of_L_-_Case_(2.B.II).pdf](<./Figures/15-Logarithmic_Nyquist_of_L_-_Case_(2.B.II).pdf>)
- [15-Logarithmic_Nyquist_of_L_-_Case_(2.B.III-gain_mod).pdf](<./Figures/15-Logarithmic_Nyquist_of_L_-_Case_(2.B.III-gain_mod).pdf>)
- [15-Logarithmic_Nyquist_of_L_-_Case_(2.B.III-gain_mod-FfP2).pdf](<./Figures/15-Logarithmic_Nyquist_of_L_-_Case_(2.B.III-gain_mod-FfP2).pdf>)
- [15-Logarithmic_Nyquist_of_L_-_Case_(2.B.III-gain_mod-butter-FfP2).pdf](<./Figures/15-Logarithmic_Nyquist_of_L_-_Case_(2.B.III-gain_mod-butter-FfP2).pdf>)
- [16-Rlocus_L_-_Case_(2.B.II)_ZOOM.pdf](<./Figures/16-Rlocus_L_-_Case_(2.B.II)_ZOOM.pdf>)
- [16-Rlocus_L_-_Case_(2.B.III-gain_mod)_ZOOM.pdf](<./Figures/16-Rlocus_L_-_Case_(2.B.III-gain_mod)_ZOOM.pdf>)
- [16-Rlocus_L_-_Case_(2.B.III-gain_mod-FfP2)_ZOOM.pdf](<./Figures/16-Rlocus_L_-_Case_(2.B.III-gain_mod-FfP2)_ZOOM.pdf>)
- [16-Rlocus_L_-_Case_(2.B.III-gain_mod-butter-FfP2)_ZOOM.pdf](<./Figures/16-Rlocus_L_-_Case_(2.B.III-gain_mod-butter-FfP2)_ZOOM.pdf>)
- [2-R(1.B)-Robustness_M-R_T_bode_plot_of_(1.B)_system.pdf](<./Figures/2-R(1.B)-Robustness_M-R_T_bode_plot_of_(1.B)_system.pdf>)
- [2-R(1.B)-Robustness_R_T_bode_plot_of_(1.B)_system.pdf](<./Figures/2-R(1.B)-Robustness_R_T_bode_plot_of_(1.B)_system.pdf>)
- [2-R(2.B.III)-Robustness_M-R_T_bode_plot_of_(.pdf](<./Figures/2-R(2.B.III)-Robustness_M-R_T_bode_plot_of_(.pdf>)
- [2-R(2.B.III)-Robustness_R_T_bode_plot_of_(.pdf](<./Figures/2-R(2.B.III)-Robustness_R_T_bode_plot_of_(.pdf>)
- [8-Step_Response_-_Case_(1.B).pdf](<./Figures/8-Step_Response_-_Case_(1.B).pdf>)
- [9-Bode_L_-_Case_(1.B).pdf](<./Figures/9-Bode_L_-_Case_(1.B).pdf>)
- [Bode_-_Loop_function_of_base_system.pdf](<./Figures/Bode_-_Loop_function_of_base_system.pdf>)
- [Margin_FfP2.pdf](<./Figures/Margin_FfP2.pdf>)
- [NyquistLogFf.pdf](<./Figures/NyquistLogFf.pdf>)
- [Rlocus_-_Case_(1.A.I).pdf](<./Figures/Rlocus_-_Case_(1.A.I).pdf>)
- [Rlocus_-_Case_(1.A.II).pdf](<./Figures/Rlocus_-_Case_(1.A.II).pdf>)
- [Rlocus_-_Loop_function_of_base_system_-_negative_gain.pdf](<./Figures/Rlocus_-_Loop_function_of_base_system_-_negative_gain.pdf>)
- [Rlocus_-_Loop_function_of_base_system_-_positive_gain.pdf](<./Figures/Rlocus_-_Loop_function_of_base_system_-_positive_gain.pdf>)
- [Science_TCV-purple-plasma-visible-light-cam.pdf](<./Figures/Science_TCV-purple-plasma-visible-light-cam.pdf>)
- [Step_-_closed_Loop_function_of_base_system.pdf](<./Figures/Step_-_closed_Loop_function_of_base_system.pdf>)
- [Step_Response_-_Case_(1.A.II).pdf](<./Figures/Step_Response_-_Case_(1.A.II).pdf>)
- [TOKAMAK_DIAGRAM.svg](<./Figures/TOKAMAK_DIAGRAM.svg>)
