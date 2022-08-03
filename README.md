# accretion-turbulence-ps
Turbulence power spectrum due to accretion of CGM gas by galactic disk

This is a code which tries to reproduce the work outlined in: ```Gas accretion onto galaxies and Kelvin-Helmholtz turbulence (ArXiv: 2207.08685)```
This predicts a power law in velocity power spectrum of with slope of around -1 at scales of around 8 kpc. This can serve as an indirect probe of accretion of CGM onto galactic disk. However, another work ```Turbulent universal galactic Kolmogorov velocity cascade over 6 decades (ArXiv: 2204.13760)``` reports an observed power spectrum slope close to Kolmogorov in a large range of scales. It would pehaps be interesting to study these two work together more closely.

The first order ordinary differential equation solved in this code that gives the rate of energy cascade from larger to smaller eddies due to turbulence is,
![image](https://user-images.githubusercontent.com/39578361/182538488-6ea415a8-2d4c-4720-9593-79a445645585.png)
