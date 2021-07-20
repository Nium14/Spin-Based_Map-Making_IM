# Spin-based removal of instrumental systematics in 21cm intensity mapping surveys

> Author: Nialh McCallum <br/> 
> Year: 2021 <br/>
> Email: [nialh.mccallum@googlemail.com](mailto:nialh.mccallum@googlemail.com)
---

**DOI**: [TBC](http://dx.doi.org/TBC)   
**arXiv**: [TBC](https://arxiv.org/abs/TBC)


This is the code repository for the paper Spin-based removal of instrumental systematics in 21cm intensity mapping surveys. All routines are implemented in Python (>=3.5) and we provide scripts for reproducing the results in the paper.


## Contents

1. [Data](#Data)
2. [Requirements](#Requirements)
3. [Scripts](#Scripts)
4. [Acknowledgements](#Acknowledgements)
5. [License](#License)

### Data

We make use of the "fast" simulations of the <a href="http://intensitymapping.physics.ox.ac.uk/CRIME.html" target_="blanck">CRIME</a> software.

The simulated input data is provided by them <a href="http://intensitymapping.physics.ox.ac.uk/Simulations/fast1/" target_="blanck">here</a>.


### Requirements

The scripts included have the following prerequisites making extensive use of:
* <a href="https://github.com/isab3lla/gmca4im" target_="blank">gmca4im</a>
* <a href="https://pypi.org/project/pyephem/" target_="blank">pyephem</a>
* <a href="https://www.python.org/" target_="blank">Python</a> (required >=3.5)
* <a href="https://matplotlib.org/" target_="blank">Matplotlib</a> (recommended >=3.1.1)
* <a href="http://www.numpy.org/" target_="blank">NumPy</a> (recommended >=1.16.4)
* <a href="https://github.com/healpy/" target_="blank">healpy</a> (recommended >=1.12.9)

### Scripts

1. [Scanning Strategy](./ScanDataGeneration/ScanStrategyGen.py)

This script performs the scanning strategy simulations to generate the pointing data.

2. [Gain Mismatch TOD](./TODSimulations/GainMismatchSimulation.py)

This script performs the TOD simulations of the Gain Mismatch systematic.

3. [Beam Squint TOD](./TODSimulations/BeamSquintSimulation.py)

This script performs the TOD simulations of the Beam Squint systematic.

4. [Foreground Removal](./ForegroundCleaning/GainFGCleaning.py)

This script performs foreground removal on the output maps from the Gain Inluding TOD simulations.

5. [Foreground Removal](./ForegroundCleaning/BeamSquintFGCleaning.py)

This script performs foreground removal on the output maps from the Beam Squint Inluding TOD simulations.

6. [Derivative Fields](./InputandDerivativePlots/derivativefieldsvsfreq.py)

This script generates plots of the derivative fields as a function of frequency.

7. [Input Fields](./InputandDerivativePlots/PlotInputSignals.py)

This script generates plots of the input fields as a function of frequency.


### Acknowledgements

If you use any of this material in your research, please cite the <a href="https://arxiv.org/abs/2107.08058" target_="blanck">companion paper</a>.


### License

GNU General Public License v3.0, read the [LICENSE](LICENSE) for more information.
