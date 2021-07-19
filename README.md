# Intensity-Mapping

> Author: Nialh McCallum <br/> 
> Year: 2021 <br/>
> Email: [nialh.mccallum@googlemail.com](mailto:nialh.mccallum@googlemail.com)
---

**DOI**: [TBC](http://dx.doi.org/TBC)   
**arXiv**: [TBC](https://arxiv.org/abs/TBC)


This is the code repository for the paper Spin-based removal of instrumental systematics in 21cm intensity mapping surveys. All routines are implemented in Python (>=3.5) and we provide scripts for reproducing the results in the paper.


## Contents

1. [Data](#Data)
1. [Requirements](#Requirements)
1. [Scripts](#Scripts)
1. [Acknowledgements](#Acknowledgements)
1. [License](#License)

### Data

We make use of the "fast" simulations of the <a href="http://intensitymapping.physics.ox.ac.uk/CRIME.html" target_="blanck">CRIME</a> software.

The simulated input data is provided at <a href="http://intensitymapping.physics.ox.ac.uk/Simulations/fast1/" target_="blanck">here</a>.


### Requirements

The scripts included have the following prerequisites making extensive use of:
* <a href="https://github.com/isab3lla/gmca4im" target_="blank">gmca4im</a> (recommend >=1.12.9)
* <a href="https://www.python.org/" target_="blank">Python</a> (require >=3.5)
* <a href="https://matplotlib.org/" target_="blank">Matplotlib</a> (recommend >=3.1.1)
* <a href="http://www.numpy.org/" target_="blank">NumPy</a> (recommend >=1.16.4)
* <a href="https://www.scipy.org/" target_="blank">SciPy</a> (recommend >=1.3.0)
* <a href="https://github.com/healpy/" target_="blank">healpy</a> (recommend >=1.12.9)


### Scripts

1. [Gain Mismatch TOD](./scripts/tbc.py)

This script performs the TOD simulations of the Gain Mismatch systematic.

2. [Beam Squint TOD](./scripts/tbc.py)

This script performs the TOD simulations of the Beam Squint systematic.

3. [Foreground Removal](./scripts/tbc.py)

This script performs uses  TOD simulations of the Gain Mismatch systematic.

### Acknowledgements

If you use any of this material in your research, please cite the <a href="http://dx.doi.org/TBC" target_="blanck">companion paper</a>.


### License

GNU General Public License v3.0, read the [LICENSE](LICENSE) for more information.
