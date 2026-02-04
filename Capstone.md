**Status as of January 30, 2026**

The project followed a structured experimental workflow beginning with Mach–Zehnder modulator characterization, then extending to evaluation of the TFLN modulation platform, and culminating in link level assessment of IM/DD performance and signal integrity.

### Phase 1 — Modulator Characterization

First, let us understand how a Mach–Zehnder modulator (MZM) operates. An MZM works by splitting an incoming optical signal into two interferometer arms. A relative phase shift is then introduced between the two paths. This phase shift can either be applied entirely to one arm as φ, or symmetrically as ±φ/2 across both arms. A visual illustration of the Mach–Zehnder interferometer used in this setup is shown below.

<img src="/images/mzmI.jpg" width="800">

In our experiment, the MZM splits the optical signal and applies a +φ/2 phase shift to one arm and a -φ/2 phase shift to the other. By controlling the relative phase between the two optical waves, constructive or destructive interference occurs when the signals are recombined. This interference determines the resulting optical intensity, which can be detected by a photodiode and interpreted as an amplitude-modulated signal.

The phase shift is induced through the **Pockels electro-optic effect**. When a voltage is applied across the electro-optic material, an electric field is generated that causes a distortion of the crystal lattice. This distortion alters the polarization response of the material, resulting in a change in the refractive index. The change in refractive index directly modifies the optical phase accumulated along the modulator arms.

The optical phase accumulated over a waveguide of length \( L \) is given by

$$
\varphi = \frac{2\pi n L}{\lambda}
$$

A change in the refractive index therefore produces a corresponding phase shift

$$
\Delta \varphi = \frac{2\pi \Delta n L}{\lambda}
$$

This relationship shows that electrically induced changes in the refractive index lead to controlled phase modulation, enabling the MZM to convert an electrical signal into an intensity-modulated optical signal.


The Mach–Zehnder modulator (MZM) was characterized to establish key device parameters including insertion loss (IL) and half-wave voltage (Vπ). These measurements provided a baseline understanding of modulation efficiency and biasing requirements for high-speed operation.

The test setup plan was as shown below:

<img src="/images/yooo.jpg" width="400">


The normalized Transfer Function (TF) is shown below

<img src="/images/mzmr1.jpg" width="600">

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

<img src="/images/20G.jpg" width="600">

**32Gbps**
**SNR 16.40**

<img src="/images/32G.jpg" width="600">

**50Gbps**
**SNR 13.86**

<img src="/images/50G.jpg" width="600">

**60Gbps**
**SNR 10.13**

<img src="/images/60G.jpg" width="600">

**70Gbps**
**SNR 7.23**

<img src="/images/70G.jpg" width="600">


The SNR decreases with increasing data rate as higher bandwidth is required to preserve fast signal transitions, increasing susceptibility to noise and intersymbol interference (ISI). A significant SNR degradation occurs beyond 40 Gbps due to the 40 GHz frequency response limit of the modulator.


### Phase 3 — IM/DD MZM

For this experiment, an IM/DD MZM setup was evaluated at three data rates over fiber lengths of 0 km and 5 km.

The test setup is shown below

<img src="/images/imddmzm.jpg" width="600">

The resulting optical eye diagrams for the 0km fiber optic capbles are shown below for 
10, 21, 28 Gbps respectively

<u><b>0km</b></u>

**10Gbps**

<img src="/images/0_10G.jpg" width="600">

**21Gbps**

<img src="/images/0_21G.jpg" width="600">

**28Gbps**

<img src="/images/0_28G.jpg" width="600">


Although the eye diagrams remain reasonably open, they exhibit noticeably degraded quality compared to the previously measured electrical eye diagrams. Increased timing jitter is observed, resulting in thicker eye traces and reduced clarity.

<u><b>5km</b></u>

**10Gbps**

<img src="/images/5_10G.jpg" width="600">

**21Gbps**

<img src="/images/5_21G.jpg" width="600">

**28Gbps**

<img src="/images/5_28G.jpg" width="600">

The measured eye diagrams no longer resemble ideal eye openings, which is expected since this experiment was performed in the C-band. At these wavelengths, optical signals are more susceptible to **chromatic dispersion** over kilometer-scale fiber lengths. For this reason, the O-band is commonly used in industry for IM/DD systems, as it operates near the fiber’s zero-dispersion wavelength. The C-band was used in this experiment due to its near-zero optical attenuation in standard single-mode fiber (SMF).

Chromatic dispersion occurs when different wavelengths within a light signal propagate at slightly different velocities through an optical fiber, causing the signal to spread out in time and become distorted.

In this system, chromatic dispersion arises from two main mechanisms. First, the optical source is not perfectly monochromatic. Although the laser is centered at 1550 nm, small wavelength variations exist around the carrier. Over long fiber lengths, slower wavelength components from earlier symbols can overlap with faster components from later symbols, resulting in ISI.

Second, chromatic dispersion is introduced by amplitude modulation itself. The optical carrier frequency is given by

$$
f_0 = \frac{c}{\lambda}
$$

$$
= \frac{3.0 \times 10^8}{1550 \times 10^{-9}}
$$

$$
\approx 193\ \text{THz}
$$

Let the modulated optical electric field be expressed as

$$
E(t) = A_0 \left[ 1 + m \cos(\Omega t) \right] \cos(\omega_0 t)
$$

Expanding the expression yields

$$
E(t) = A_0 \cos(\omega_0 t) + m A_0 \cos(\Omega t)\cos(\omega_0 t)
$$

Using the identity

$$
\cos(a)\cos(b) = \frac{1}{2}\left[\cos(a+b) + \cos(a-b)\right]
$$

the final form of the electric field becomes

$$
E(t) =
A_0 \cos(\omega_0 t)
+
\frac{m A_0}{2}\cos\big((\omega_0+\Omega)t\big)
+
\frac{m A_0}{2}\cos\big((\omega_0-\Omega)t\big)
$$

This shows that amplitude modulation generates upper and lower sidebands at different frequencies (and therefore different wavelengths), which propagate at different speeds in the fiber and lead to chromatic dispersion.

One may wonder why the O-band does not significantly suffer from chromatic dispersion. Chromatic dispersion is the combined effect of two mechanisms: material dispersion and waveguide dispersion.

At shorter wavelengths, material dispersion is negative while waveguide dispersion is close to zero. As the wavelength increases, material dispersion increases, whereas waveguide dispersion decreases. In the O-band, these two effects approximately cancel each other, resulting in a chromatic dispersion that is near zero.

As the wavelength continues to increase into the C-band, material dispersion grows more rapidly than waveguide dispersion decreases. Consequently, the cancellation no longer occurs, and the resulting chromatic dispersion becomes non-zero. This behavior is illustrated in the figure below.

<img src="/images/disp.jpg" width="600">


### Next Phase - Upcoming work: TFLN Waveguide Crossing Chip & Full IM/DD TX System

The MZM was integrated into the transmission setup with a laser source, bias voltage control, BPG input, and DCA receiver. Eye diagrams were captured after transmission over 2 km and 5 km fiber spans to evaluate jitter accumulation and signal degradation with distance.

