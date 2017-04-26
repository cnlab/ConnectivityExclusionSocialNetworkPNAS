# ConnectivityExclusionSocialNetworkPNAS
This repository includes data and code for 
*Schmälzle, R., Brook O’Donnell, M., Garcia, J.O., Cascio, C.N., Bayer, J., Bassett, D.S., Vettel, J. & Falk, E.B. (2017). Brain connectivity dynamics during social interaction reflect social network structure. Proceedings of the National Academy of Sciences.* [![DOI](http://www.ralfschmaelzle.net/wp-content/plugins/papercite/img/pdf.png)](https://doi.org/10.1073/pnas.1511477112) 
[External Link]


***

#### Analysis Code
* The [notebook to reproduce the analysis of the main effect of exclusion](https://github.com/nomcomm/ConnectivityExclusionSocialNetworkPNAS/blob/master/notebooks/02_main_analysis/Schmaelzle_ConnectivitySociaExclusion.ipynb)
* The [notebook to reproduce the analysis of the association between social network density and brain connectivity](https://github.com/nomcomm/ConnectivityExclusionSocialNetworkPNAS/blob/master/notebooks/02_main_analysis/Schmaelzle_ConnectivityDensity.ipynb)

#### Data
* [Extracted data matrices for the a-priori-networks](https://github.com/nomcomm/ConnectivityExclusionSocialNetworkPNAS/blob/master/data/cconnectivity_matrices)
* [Extracted data matrices for the Power-264-node parcellation](https://github.com/nomcomm/ConnectivityExclusionSocialNetworkPNAS/blob/master/data/cconnectivity_matrices_264)
* [Datasheet with the social network density metrics](https://github.com/nomcomm/ConnectivityExclusionSocialNetworkPNAS/blob/master/data/datasheets/pID_social_networks.csv)


#### Dependencies
* Python

[Anaconda](http://continuum.io/downloads) should provide you with most of what you need.

The following packages are used and we feel very indebted to their creators:
* [Project Jupyter](https://github.com/jupyter) 
* [nilearn](https://github.com/nilearn) e.g. `pip install nilearn`
* [seaborn](http://seaborn.pydata.org/)
* [numpy](http://www.numpy.org/)
* [scipy](http://www.scipy.org/)
* [matplotlib](http://matplotlib.org/)
* [pandas](http://pandas.pydata.org/)
* [bctpy](https://github.com/aestrivex/bctpy)
* [mne](https://github.com/mne-tools)

Note: if you run into errors indicating you miss a package, either enter "pip install package" in a terminal or - if in the notebook - insert a cell and write "!pip install package"

***
2017 | Ralf Schmaelzle | Matthew Brook O'Donnell 