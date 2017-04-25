# ConnectivityExclusionSocialNetworkPNAS
Data and code for 
Schmälzle, R., Brook O’Donnell, M., Garcia, J.O., Cascio, C.N., Bayer, J., Bassett, D.S., Vettel, J. & Falk, E.B. (2017). Brain connectivity dynamics during social interaction reflect social network structure. Proceedings of the National Academy of Sciences

The paper can be accessed here: [Link]

***

#### This repository includes:
#### Analysis Code used to conduct main analyses presented in the paper:
* The [notebook to reproduce the analysis of the main effect of exclusion](https://github.com/nomcomm/ConnectivityExclusionSocialNetworkPNAS/blob/master/notebooks/02_main_analysis/Schmaelzle_ConnectivitySociaExclusion.ipynb)
* The [notebook to reproduce the analysis of the association between social network density and brain connectivity](https://github.com/nomcomm/ConnectivityExclusionSocialNetworkPNAS/blob/master/notebooks/02_main_analysis/Schmaelzle_ConnectivityDensity.ipynb)

#### Data:
* [Extracted data matrices for the a-priori-networks](https://github.com/nomcomm/ConnectivityExclusionSocialNetworkPNAS/blob/master/data/clean_ts)
* [Extracted data matrices for the Power-264-node parcellation](https://github.com/nomcomm/ConnectivityExclusionSocialNetworkPNAS/blob/master/data/clean_ts_264)
* [Datasheet with the social network density metrics](https://github.com/nomcomm/ConnectivityExclusionSocialNetworkPNAS/blob/master/data/datasheets/pID_social_networks.csv)

***

#### The following packages are used and we thank their creators
* [Project Jupyter](https://github.com/jupyter)
* [nilearn](https://github.com/nilearn)
* [seaborn](http://seaborn.pydata.org/)
* [numpy](http://www.numpy.org/)
* [scipy](http://www.scipy.org/)
* [matplotlib](http://matplotlib.org/)
* [pandas](http://pandas.pydata.org/)
(not essential: bct !pip install bctpy & mne - if you need to install on the fly: enter !pip install bctpy or so in a code cell)
