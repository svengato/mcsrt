# mcsrt
#### Monte Carlo Simulation of Radiation Trapping<br>in a Dense Optically Pumped Alkali Vapor

[Optical pumping](https://en.wikipedia.org/wiki/Optical_pumping) of dense vapors is
limited by radiation trapping, in which the scattered (and generally unpolarized) light
gets reabsorbed before leaving the vapor.

This Monte Carlo simulation of the process is adapted from my thesis, which involved
[spin exchange](https://en.wikipedia.org/wiki/Spin-exchange_interaction) optical pumping.
The main insight was to detune the laser frequency away from the center of the atomic
resonance line, whose shape is a [Voigt profile](https://en.wikipedia.org/wiki/Voigt_profile).
Only atoms whose velocity Doppler shifts them into resonance are likely to absorb the
detuned photon, and will scatter it with a frequency distribution centered at the original
detuning. This reduces the likelihood of reabsorption, achieving a higher polarization in
a denser vapor and facilitating spin exchange.

In the simulation, a circularly polarized laser beam enters a cylindrical cell containing
an alkali vapor. Each incident photon scatters one or more times, and eventually escapes
from the cell. Its total spin angular momentum transfer &Delta;Jz is +&frac13; for the
initial scatter, but may be negative after multiple scatters. The equilibrium polarization
is where &lt;&Delta;Jz&gt; = 0.

#### Input

You may specify the alkali vapor (K, Rb, Cs), its polarization and temperature, and the
frequency detuning of the laser. The vapor concentration and resulting optical depth of
the cell increase nonlinearly with temperature.

#### Controls

You may run the simulation at different speeds: **Run** for continuously, **One Photon**
to pause after each incident photon escapes, and **One Scatter** to pause after each
scatter. **Pause** pauses the simulation, and **Reset** clears the statistics (as does
changing the input values).

#### Output

For now, it displays the trajectory, &Delta;Jz, and the number of scatters for the current
photon, and the mean value and probability density of these quantities across all photons.

#### (Non-)License

Copyright &copy; 1988-2018 svengato

All rights reserved.
