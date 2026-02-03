**Status as of January 30, 2026**

The project followed a structured experimental workflow beginning with Mach–Zehnder modulator characterization, then extending to evaluation of the TFLN modulation platform, and culminating in link level assessment of IM/DD performance and signal integrity.

### Phase 1 — Modulator Characterization

The Mach–Zehnder modulator (MZM) was characterized to establish key device parameters including insertion loss (IL) and half-wave voltage (Vπ). These measurements provided a baseline understanding of modulation efficiency and biasing requirements for high-speed operation.

The test setup plan was as shown below:

<img src="/images/yooo.jpg" width="400">


The normalized Transfer Function (TF) is shown below

<img src="/images/mzmr1.jpg" width="400">

**Insertion Loss Calculation**

$$
IL(\text{dB}) = P_{\text{ref}}(\text{dBm}) - P_{\text{max}}(\text{dBm})
$$

$$
= 9.63 - 2.39
$$

$$
= 7.24\ \text{dB}
$$

**Half-Wave Voltage Calculation**

$$
V_{\pi} = \left| V_{p,\text{max}} - V_{p,\text{min}} \right|
$$

$$
= \left| 2 - 7.2 \right|
$$

$$
= 5.2\ \text{V}
$$

**Extinction Ratio Calculation**

$$
P_1 = 2.39\ \text{dBm}
$$

$$
= 10^{2.39/10}
$$

$$
= 1.733\ \text{mW}
$$

$$
P_0 = -9.83\ \text{dBm}
$$

$$
= 10^{-9.83/10}
$$

$$
= 0.104\ \text{mW}
$$

$$
\mathrm{ER(dB)} = 10 \log_{10}\left(\frac{P_1}{P_0}\right)
$$

$$
= 10 \log_{10}\left(\frac{1.733}{0.104}\right)
$$

$$
= 12.22\ \text{dB}
$$



### Phase 2 — BPG Electrical Eye Diagram

A Bit Pattern Generator (BPG), synchronized to an external clock, was connected to a Bit Error Analyzer (BEA), which interfaced with a Digital Communication Analyzer (DCA) to capture electrical eye diagrams and evaluate the signal-to-noise ratio (SNR) at data rates up to 70 Gbps.

The test set up is shown below

<img src="/images/BPGGG.jpg" width="400">

The resulting electrical eye diagrams for data rates of 20, 32, 50, 60, and 70 Gbps are shown below. Since On-Off Keying (OOK) is employed, the baud rate is equal to the bit rate. In contrast, for PAM-4 modulation, the bit rate would be twice the baud rate.

**20Gbps**
**SNR 17.37**

<img src="/images/20G.jpg" width="400">

**32Gbps**
**SNR 16.40**

<img src="/images/32G.jpg" width="400">

**50Gbps**
**SNR 13.86**

<img src="/images/50G.jpg" width="400">

**60Gbps**
**SNR 10.13**

<img src="/images/60G.jpg" width="400">

**70Gbps**
**SNR 7.23**

<img src="/images/70G.jpg" width="400">


The SNR decreases with increasing data rate as higher bandwidth is required to preserve fast signal transitions, increasing susceptibility to noise and intersymbol interference (ISI). A significant SNR degradation occurs beyond 40 Gbps due to the 40 GHz frequency response limit of the modulator.


### Phase 3 — Eye Diagram and Signal Integrity Evaluation

The BPG was directly interfaced with the DCA to generate eye diagrams and evaluate baseline signal integrity. These measurements allowed for observation of timing jitter, noise, and amplitude variations independent of optical modulation effects.

### Phase 4 — Optical Link Testing with MZM

The MZM was integrated into the transmission setup with a laser source, bias voltage control, BPG input, and DCA receiver. Eye diagrams were captured after transmission over 2 km and 5 km fiber spans to evaluate jitter accumulation and signal degradation with distance.

### Next Phase — TFLN PIC Evaluation (In progress)

Upcoming work will focus on testing TFLN-based modulator chips integrated on photonic integrated circuits operating in the O-band. Special emphasis will be placed on waveguide crossings to quantify insertion loss, crosstalk, and reflection effects at data rates exceeding 20 Gbps.
