# #NPDO
get_ipython().magic(u'matplotlib inline')
from scipy.io import matlab
from scipy import stats
import os, re, nilearn, nibabel
import numpy as np
import sys
import scipy
import matplotlib.pyplot as plt
import seaborn as sns
from nilearn import image, plotting, input_data
from nilearn.connectome import ConnectivityMeasure
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
import sklearn

def plot_sns_func(adj_f, adj_uf, adj_diff):
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18,6))

    cmap = sns.diverging_palette(220, 10, as_cmap=True)
    mask = np.zeros_like(mean_connect_fair, dtype=np.bool)
    mask[np.triu_indices_from(mask)] = True

    #ax1 = fig.add_subplot(131)   
    sns.heatmap(adj_f, mask=mask, cmap=cmap, vmax=.3, vmin = -.3, annot = False,
                    square=True, xticklabels=networks_labels, yticklabels=networks_labels, 
                    linewidths=.5, cbar_kws={"shrink": .3}, ax=ax1)#, cbar=i == 0)       
    ax1.set_title('Connectivity Fair Clean')

    #ax2 = fig.add_subplot(132)   
    sns.heatmap(adj_uf, mask=mask, cmap=cmap, vmax=.3, vmin = -.3, annot = False,
                    square=True, xticklabels=networks_labels, yticklabels=networks_labels, 
                    linewidths=.5, cbar_kws={"shrink": .3}, ax=ax2)#, cbar=i == 0)       
    ax2.set_title('Connectivity Unfair Clean')

    #ax3 = fig.add_subplot(133)   
    sns.heatmap(adj_diff, mask=mask, cmap=cmap, vmax=.3, vmin = -.3, annot = False,
                    square=True, xticklabels=networks_labels, yticklabels=networks_labels, 
                    linewidths=.5, cbar_kws={"shrink": .3}, ax=ax3)#, cbar=i == 0)       
    ax3.set_title('Connectivity GM Unfair-Fair Clean Diff')
    plt.show()




def extract_confounds(path, nifti_file, plot=False):
    '''
        take a NIFTI file from a run and
        return the confounds from movement parameters
        and anatomical segmentation
        
        params:
           * path - path to directory containing nifti file  
           * nifti_file - usually the final stage of preprocessing, e.g. ad_ ... .nii
                          function assumes same folder has the realignment parameters, e.g.
                          rp_ad_.txt
           * plot (optional) - if True then plot slices of the anatomical masks
           
        returns:
            * hv_confounds - high variance confounds
            * matter_confounds - white, gray matter and csf means removed
    '''

    # get the subj id from the path
    subj_id = re.findall('/([^/]+)/BOLD',path)[0]
    
    func_filename = os.path.join(path, nifti_file)  #preprocessed data
                                         
    hv_confounds = nilearn.image.high_variance_confounds(func_filename) #high variance confounds
        
    cut_coords = (0, 0, 0)

    anatomy_dir = re.sub('BOLD/.*','anatomy', path)
    
    t1img = nibabel.load(os.path.join(anatomy_dir, 'regt1spgr.nii'))
    
    # Get Confounds: White Matter + CSF/Ventricles 
    
    matter_confounds = None   

    for c in range(2,4):  #I changed this to 2:4 to get only white matter and CSF - no gray matter...
        m_filename = 'wc'+ str(c) + 'regt1spgr.nii'
        m_filepath = os.path.join(anatomy_dir, m_filename)

        anat_img = nibabel.load(m_filepath);
        m_hand_mask = anat_img.get_data();
        
        m_hand_mask = np.where(m_hand_mask > 0.7, 1, 0)

        m_anat_img = nibabel.Nifti1Image(m_hand_mask, anat_img.affine) 
        
        if plot:
            plotting.plot_anat(m_anat_img, cut_coords=cut_coords,title=('{} {}'.format(subj_id, m_filename)));
        
    
        m_masker = input_data.NiftiLabelsMasker(labels_img = m_anat_img, 
                                                standardize=True, 
                                                memory='nilearn_cache', verbose=0)
        
        ts = m_masker.fit_transform(func_filename)
        
        if matter_confounds is None:
                        matter_confounds = ts

        else:
                matter_confounds = np.append(matter_confounds, ts, axis=1)

    
    return hv_confounds, matter_confounds



def clean_time_series(path, nifti_file, networks_coords, hv_confounds=None, matter_confounds=None, fd_thresh = None, plot=None):  
    
    # get the subj id from the path
    subj_id = re.findall('/([^/]+)/BOLD',path)[0]
    
    func_filename = os.path.join(path, nifti_file)  #preprocessed data

    conf_filename = os.path.join(path, 'rp_adrun_04.txt')  #motion confounds
        
    if matter_confounds is None or hv_confounds is None:
        hv_confounds, matter_confounds = extract_confounds(path, nifti_file)
    
    masker = input_data.NiftiSpheresMasker(networks_coords, 
                                           radius=8, 
                                           allow_overlap=True, 
                                           detrend=True, 
                                           standardize=True, 
                                           low_pass=0.12, high_pass=0.06, 
                                           t_r=2,  
                                           memory='nilearn_cache', memory_level=1); #verbose=2 by default nothing should be printed
        
    cleaned_time_series = masker.fit_transform(func_filename,confounds=[conf_filename, matter_confounds, hv_confounds])  

    if fd_thresh>0:     
        FD = compute_fd(conf_filename)
        bad_fd_vols  = np.where(FD > fd_thresh, 1,0)
    else:
        bad_fd_vols = []
   # if plot:
   #     cm = np.corrcoef(clean_time_series.T)
   #     f, ((ax1, ax2)) = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(10,5))
   #     ax1.plot(cleaned_time_series.T);
   #     ax1.set_title('Subj: {} Clean time series'.format(subj_id)) 
        
    return cleaned_time_series,  FD, bad_fd_vols
    
    
def clean_time_series_FIR(path, nifti_file, networks_coords, hv_confounds=None, matter_confounds=None, fd_thresh = None, plot=None, fir_regressors= None):  
    
    # get the subj id from the path
    subj_id = re.findall('/([^/]+)/BOLD',path)[0]
    
    func_filename = os.path.join(path, nifti_file)  #preprocessed data

    conf_filename = os.path.join(path, 'rp_adrun_04.txt')  #motion confounds
        
    if matter_confounds is None or hv_confounds is None:
        hv_confounds, matter_confounds = extract_confounds(path, nifti_file)
    
    masker = input_data.NiftiSpheresMasker(networks_coords, 
                                           radius=8, 
                                           allow_overlap=True, 
                                           detrend=True, 
                                           standardize=True, 
                                           low_pass=0.12, high_pass=0.06, 
                                           t_r=2,  
                                           memory='nilearn_cache', memory_level=1); #verbose=2 by default nothing should be printed
        
    cleaned_time_series = masker.fit_transform(func_filename,confounds=[fir_regressors, conf_filename, matter_confounds, hv_confounds])  
    print('hvconfounds')
    print(hv_confounds.shape)
    print('fir_regressors')
    print(fir_regressors.shape)

    if fd_thresh>0:     
        FD = compute_fd(conf_filename)
        bad_fd_vols  = np.where(FD > fd_thresh, 1,0)
    else:
        bad_fd_vols = []
   # if plot:
   #     cm = np.corrcoef(clean_time_series.T)
   #     f, ((ax1, ax2)) = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(10,5))
   #     ax1.plot(cleaned_time_series.T);
   #     ax1.set_title('Subj: {} Clean time series'.format(subj_id)) 
        
    return cleaned_time_series,  FD, bad_fd_vols    


def get_TR_ranges(path, subject):
   # Get the Onsets for Fair and Unfair Runs
    behav = matlab.loadmat(os.path.join(path ,subject,'behav','{}.mat'.format(subject)))  
    conds = behav['trials'][0]
  
    #print(subjs[n])
    if len(conds) > 5:     
        fair_trials = conds[4][0];           fair_onset = fair_trials[0]['onset']
        unfair_trials = conds[5][0];         unfair_onset = unfair_trials[0]['onset']
              
        fair_extent = (fair_trials[0]['onset'], fair_trials[-1]['onset'] + fair_trials[-1]['rt'])
        unfair_extent = (unfair_trials[0]['onset'], unfair_trials[7]['onset']+unfair_trials[7]['rt'], unfair_trials[8]['onset'], unfair_trials[-1]['onset'] + unfair_trials[-1]['rt'])
        
        fair_extent_TR = int(fair_extent[0]/2), int(fair_extent[-1]/2)
        unfair_extent_TR = int(unfair_extent[0]/2), int(unfair_extent[-1]/2)      
    else:
        print('Error - subject' + subject)
    return fair_extent_TR, unfair_extent_TR


def compute_fd(conf_filename):
    # This function is based on Powers 2012 as well as Bramila tools and 
    # Poldracks fMRI QA script. Important note: It assumes the input is a 
    # SPM-realignment parameter file. FSL uses a different ordering and 
    # thus cannot be used blindly... Beware!

    motpars = pd.read_csv(conf_filename, header=None, delimiter=r"\s+")
    motpars = motpars.as_matrix()

    # compute absolute displacement
    dmotpars=np.zeros(motpars.shape)
    
    dmotpars[1:,:]=np.abs(motpars[1:,:] - motpars[:-1,:])
    
    # convert rotation to displacement on a 50 mm sphere
    # mcflirt returns rotation in radians
    # from Jonathan Power:
    # The conversion is simple - you just want the length of an arc that a rotational
    # displacement causes at some radius. Circumference is pi*diameter, and we used a 5
    # 0 mm radius. Multiply that circumference by (degrees/360) or (radians/2*pi) to get the 
    # length of the arc produced by a rotation.  
    headradius=50
    disp=dmotpars.copy()
    disp[:,3:6]=np.pi*headradius*2*(disp[:,3:6]/(2*np.pi))
    
    FD=np.sum(disp,1)
      
    return FD

def r_to_z(r):
    # fisher transform
    z=0.5*np.log((1.0+r)/(1.0-r))
    z[np.where(np.isinf(z))]=0
    z[np.where(np.isnan(z))]=0

    return z

def z_to_r(z):
    # inverse transform
    return (np.exp(2.0*z) - 1)/(np.exp(2.0*z) + 1)    
 
def intersect(a, b):
    """ return the intersection of two lists """
    return list(set(a) & set(b))

   
def extract_network_means(adj):
    # helper function for extraction of within-network-metric 
    
    salX = np.array([0,3])
    dmnX = np.array([3,10])

    sal = adj[ salX[0]:salX[1] , salX[0]:salX[1] ]
    dmn = adj[ dmnX[0]:dmnX[1] , dmnX[0]:dmnX[1] ]
    betw = adj[ (dmnX[0]):dmnX[1] , salX[0]:salX[1]  ]

    sal = sal[np.tril_indices_from(sal, k = -1)] 
    dmn = dmn[np.tril_indices_from(dmn, k = -1)] 
       
    sal_m = np.mean(sal)
    dmn_m = np.mean(dmn) 
    betwm = np.mean(betw) 
    res = np.zeros(3)
    res[0] = sal_m
    res[1] = dmn_m
    res[2] = betwm

    return res
    
    
def extract_network_means_3acc(adj):
    # helper function for extraction of within-network-metric 
    
    salX = np.array([0,5])
    dmnX = np.array([5,12])

    sal = adj[ salX[0]:salX[1] , salX[0]:salX[1] ]
    dmn = adj[ dmnX[0]:dmnX[1] , dmnX[0]:dmnX[1] ]
    betw = adj[ (dmnX[0]):dmnX[1] , salX[0]:salX[1]  ]

    sal = sal[np.tril_indices_from(sal, k = -1)] 
    dmn = dmn[np.tril_indices_from(dmn, k = -1)] 
       
    sal_m = np.mean(sal)
    dmn_m = np.mean(dmn) 
    betwm = np.mean(betw) 
    res = np.zeros(3)
    res[0] = sal_m
    res[1] = dmn_m
    res[2] = betwm

    return res    
    
    
def extract_network_means_brainstem(adj):

    salX = np.array([0,7])
    dmnX = np.array([7,14])

    sal = adj[ salX[0]:salX[1] , salX[0]:salX[1] ]
    dmn = adj[ dmnX[0]:dmnX[1] , dmnX[0]:dmnX[1] ]
    betw = adj[ (dmnX[0]):dmnX[1] , salX[0]:salX[1]  ]

    sal = sal[np.tril_indices_from(sal, k = -1)] 
    dmn = dmn[np.tril_indices_from(dmn, k = -1)] 
       
    sal_m = np.mean(sal)
    dmn_m = np.mean(dmn) 
    betwm = np.mean(betw) 
    res = np.zeros(3)
    res[0] = sal_m
    res[1] = dmn_m
    res[2] = betwm

    return res   
    
    
    
def extract_network_means_brainstemsms(adj):

    salX = np.array([0,9])
    dmnX = np.array([9,16])

    sal = adj[ salX[0]:salX[1] , salX[0]:salX[1] ]
    dmn = adj[ dmnX[0]:dmnX[1] , dmnX[0]:dmnX[1] ]
    betw = adj[ (dmnX[0]):dmnX[1] , salX[0]:salX[1]  ]

    sal = sal[np.tril_indices_from(sal, k = -1)] 
    dmn = dmn[np.tril_indices_from(dmn, k = -1)] 
       
    sal_m = np.mean(sal)
    dmn_m = np.mean(dmn) 
    betwm = np.mean(betw) 
    res = np.zeros(3)
    res[0] = sal_m
    res[1] = dmn_m
    res[2] = betwm

    return res  
 
    
def extract_network_means_264(adj, networks_systems_vector):
    # helper function for extraction 

    # assess strenght within
    num_networks = len(np.unique(networks_systems_vector))
    mean_result =[[0]*2 for i in range(num_networks + 1)]
    
    for curr_network in range(num_networks):
        #print('curr_network: ' + str(curr_network))
        curr_network_index = np.unique(networks_systems_vector)[curr_network]
        #print('curr_network_index: ' + str(curr_network_index))
        relevant_indices = np.where(networks_systems_vector == curr_network_index )[0]
        #print('relevant_indices: ' + str(relevant_indices))
        count = 0
        sum = 0
        curr_network_length = len(relevant_indices)  
        for x in range(curr_network_length):
            x_index = relevant_indices[x]
            #print('#####')            
            #if (curr_network == 0):            
               # print('x_index: ' + str(x_index))
        
            for y in range(x+1, curr_network_length,1):  
                y_index = relevant_indices[y]
                #if (curr_network == 0):   
                #    if (x < 2):
                #        print('y_index: ' + str(y_index))
                curr_value = adj[x_index,y_index]
                sum = sum + curr_value
                count = count +1
        mean_result[curr_network][0] = curr_network_index
        print('sum: ' + str(sum))
        print('count: ' + str(count))
        
        mean_result[curr_network][1] = sum/count  
        print('mean_result: ' + str(mean_result[curr_network][1] ))
        print('####next network####')
    # assess strenght between
    count = 0
    sum = 0
    for curr_x_node in range(len(networks_systems_vector)):
        curr_x_node_network = networks_systems_vector[curr_x_node]
        for curr_y_node in range(len(networks_systems_vector)):
            curr_y_node_network = networks_systems_vector[curr_y_node]
            if (curr_x_node_network != curr_y_node_network):
               # print('between networks')
               # print('network x' + str(curr_x_node_network))    
               # print('network y' + str(curr_y_node_network))    
                count = count +1
                curr_value = adj[curr_x_node, curr_y_node]
                sum = sum + curr_value
    mean_result[num_networks ][0] = 99
    mean_result[num_networks ][1] = sum/count    
                 
    return mean_result


def plot_data_fd(timeseries, FD, FD_thresh, subj):
    # takes input timeseries [nsamples, nnodes] and FD parameters as input
    # have to be same length!
    
    fig = plt.figure()
    fig.suptitle(subj, fontsize=12)

    ax1 = plt.subplot(2, 1, 1)
    ax1.plot(FD)
    ax1.axhline(y=FD_thresh, color = 'red')           
    ax1.set_title('Framewise Displacement')
    ax1.set_ylabel('FD')
    ax1.set_ylim(0,2)
    ax1.set_xticklabels([])

    ax2 = plt.subplot(2, 1, 2, sharex = ax1)
    cmap = plt.cm.jet
    ax2.imshow(timeseries.T, interpolation="none", cmap = 'jet', aspect='auto', clim=(-3, 3))
    
    #plt.colorbar()
    
    ax2.grid(b=False) 
    ax2.set_ylim(0,len(timeseries.T))
    ax2.set_title('BOLD Values')
    ax2.set_ylabel('Node #')
    ax2.set_xlabel('Time [TR]')
    plt.show()


def create_fir(path_to_logfile, n_scans, tr, fir_delays ):  
        """
        Finite Impulse Response model applied to remove residual task-effects from button-pressing...
        """
        import numpy as np
        import os
        import nibabel as nib
        import pandas as pd
        import matplotlib.pyplot as plt
        from os import mkdir, path
        from nilearn import plotting
        from nistats.glm import FirstLevelGLM
        from nistats.design_matrix import make_design_matrix, plot_design_matrix
        from nistats import datasets
        from nistats.hemodynamic_models import compute_regressor
        import re
        
        
        infile = path_to_logfile 
        
        important2 = []
        keep_phrases = ["Keypress"]
        
        with open(infile) as f:
            f = f.readlines()
        
        for line in f:
            for phrase in keep_phrases:
                if phrase in line:
                    #ls = line.split('\t')[0]
                    important2.append(line.split('\t')[0])
                    break
        
        print(important2)
        del phrase, line, infile, keep_phrases, f,  
        
        
        #tr=2
        results = map(float, important2)
        button_presses = (np.asarray(results))
        del results, important2
         
        ### Analysis parameters #######################################
        #n_scans = 189
        frame_times = np.arange(n_scans) * tr
        conditions = ['button'] * len(button_presses)
        onsets =     button_presses.tolist()
        duration = 1. * np.ones(len(conditions))
        paradigm = pd.DataFrame({'name': conditions, 'onset': onsets, 'duration': duration})
        
        design_matrix = make_design_matrix(frame_times, paradigm, hrf_model = 'FIR', fir_delays=fir_delays) 
        
        plot_design_matrix(design_matrix)
        plt.show()
        
        dma = np.asarray(design_matrix)
        print(dma)
        regressor_fir = dma[:,0:4]
        print('###')
        print(regressor_fir)
        
        return regressor_fir
