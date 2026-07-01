
# Introduction to TULIPAN
TULIPAN (TUbulin LIgand Pocket ANalysis) is a molecular dynamics analysis protocol for the characterization of long-lived ligand binding modes on multisite proteins. It has been developed as part of the publication "Ligand Binding Free Energy Landscapes at the Tubulin Colchicine Site from Coarse-Grained Metadynamics". 
# Usage
## Part 1: Extraction of binding events with TULIPAN
Before running the analysis, the user should source some important file and set some variables to their appropriate value.
The TULIPAN protocol is structured as a bash script and requires both GROMACS tools and some external functions defined in a separate bash file. Therefore, the first step is to define the appropriate path for

    source "/path/to/gromacs/bin/GMXRC"
    source "/path/to/functions/functions_evanalysis_v3.sh"

The variable definition section is reported at the beginning of the script and can be accessed with a text editor. In particular, variables related to the "Workspace definition" and to the "Description of the biosystem" are highly system dependent.

    ### Variable definition        
    ## Workspace definition
    # Name of the output directory
    workdir="pose_analysis_cgmd"
    # Limit on parallel tasks    
    N=4    
    ## Description of the biosystem   
    # Ligand name    
    ligname=COMB    
    #ligand composition (how many beads)    
    lig_comp=11
    ###

 - `workdir` defines the name of the directory where all the analysis file will be stored;
 - `N` defines the maximum number of parallel tasks allowed to run. This depends on the architecture of the CPU. A low value (2-4) should be fine for most CPUs. For systems with multiple cores, a higher value (12-24) can speed up the analysis;
 - `ligname`: name of the ligand in the `conf.gro` file ( e.g. COMB for combretastatin-A4);
 - `lig_comp`: number of CG-beads associated to each ligand molecule;

 
The first step of the analysis can be performed with:

    bash tulipan.sh -f listfile.txt -m min.mdp -c conf.gro -s topol.tpr -p top_onelig.top -l prot_onelig.tpr

The following input options are required
 -  -f  [<.txt>] (`listfile.txt`) : list of trajectories to be analyzed (with either relative or absolute path);
 - -m [<.mdp>] (`min.mdp`): minimal GROMACS mdp file to be used by the protocol to generate .tpr files;
 - -c [<.gro>] (`conf.gro`): structure file in the gro format with the same content and topology as the trajectories;
 - -s [<.tpr>] (`topol.tpr`): GROMACS executable generated from conf.gro. This is necessary for the GROMACS tools that require information about the masses of the atoms contained in the trajectories in order to perform the analysis;
 - -p [<.top>] (`top_onelig.top`): topology file with the parameters necessary for a protein-single ligand molecule system;
 - -l [<.tpr>] (`prot_onelig.tpr`): GROMACS executable obtained from a protein-single ligand molecule system.
 
 ## Part 2: Classification of binding events
 This second analysis script will build upon the results of the first step to complete the classification of the binding modes and will try to estimate the residence time of the various ligand poses.
 Before running the analysis, please make sure that all the relevant paths and variables are set correctly.
 In particular, the following files should be defined (with the absolute path):


    #Files for fitting procedure with python
    blankfit="/path/to/functions/temp_fit_prob.py"
    blankfitdens="/path/to/functions/temp_fit_dens.py"
    temphisto="/path/to/functions/temp_multihisto_v2.py"
   For visualization with VMD of the final poses, the following templates should be provided:

    #Files for VMD visualization
    refbondpdb="path/to/prot_onelig.pdb"
    blank="/path/to/templates/blank_v3.tcl"
    blankpdb="/path/to/templates//blankpdb_v3.tcl"
   Where:

 - `refbondpd` is a PDB file of the protein with only ligand molecule and explicit CONECT records to be used for visualization of bonds between CG beads.
 - `blank` template `.tcl` file for loading of the top poses and associated concatenated trajectory
 - `blankpdb` template `.tcl` file for loading just the PDB of the top poses (recommended for visualization, as it occupies less RAM)

Other variables that should be checked before running the program are reported at the beginning of the script:

    ###Variable definition    
    ##Workspace definition    
    workdir="pose_analysis_cgmd"    
    N=4    
    ##Biosystem    
    ligname=COMB
    
`method` and `thr` arrays should also be set exactly as in the first script
    
 At the bottom of the definition section, some important aspects of the VMD visualization are specified

 
    extlimit="10"
    vmdcolors=("0"  "1"  "2"  "3"  "4"  "7"  "9"  "10"  "11"  "12")
In particular, `extlimit` regulates how many poses will be included in the final  `.tcl` files for VMD visualization, while the array `vmdcolors` specifies which color should be attributed to each ranked pose (therefore must have the same number of elements as `extlimit`). Please consult VMD's documentation for further info on the color sequence.

To run the analysis,  execute:

    bash tulipan_clust_assign.sh
    
###Python 3 Requirements
The fitting procedure for the residence time estimation requires some basic Python libraries. In particular the modules `numpy`, `matplotlib` and `scipy` will be used.
    
# Technical notes
- Some functions of the script rely on AWK for processing of the text results from GROMACS tools. In particular GNU-AWK is required, since some functionalities rely on the ARGIND keyword;
- The script has been tested with trajectories having a framerate of 1 frame/ns. Different framerates should be applicable by changing the `dt` variable in the `tulipan_clust_assign.sh` (e.g `dt=2` for 1 frame every 2 ns of MD), although this has note been extensively tested;
- During the analysis, the script attempts to generate `.tpr` files for the cluster centroids. Please ensure that the reference files provided via the `-m` and `-p` flags can successfully generate a `.tpr` file without errors. 
You can verify your files beforehand by running a test command such as:

```bash
gmx grompp -f mdp.mdp -p topol_onelig.top -c prot_onelig.gro -o test.tpr
```

# Analysis of the results
## Part 1 
As the analysis progresses, the results will be stored according to this general architecture

    pose_analysis_cgmd/
    └── rep1/
        └── COMB872/
            └── ev_141/
                └── subevents_gromos_0.4/
                    └── sub_1/
First, the protocol creates *n* `rep*/` folders corresponding to the number of trajectories described in `listfile.txt` and performs a ligand-protein contact analysis.
In each replica directory, *k* ligand folders will be generated. For every ligand-trajectory combination, a series of subfolders denoting raw binding events with duration higher than `ll_thr` (typically 100 ns) will be created. For instance, the set of events reported in the table below for COMB872 / replica-1 will yield only two long-lived event folders (`rep1/COMB872/ev_2` and `rep1/COMB872/ev_4` )
| Event | Duration |
|--|--|
|1  | 80 ns |
|2|110 ns|
|3|52 ns|
|4|210 ns|

A `gmx cluster` analysis will be applied to all the raw long-lived events to select highly localized poses. `method` and `cutoff` parameters for `gmx cluster` are defined by the bash arrays defined at the beginning of the script. By the default, only `gromos / 0.4 nm` combination is considered, but the user can also screen other method-cutoff combinations in parallel and compare the results, For instance: 

    method=("gromos" "single-linkage")
    thr=("0.4" "0.2")
Every long-lived event will therefore generate one (or more) `subevents_${method[i]}_${thr[i]} ` folders, yielding several "subevent" directories describing highly localized ligand poses

A series of smoothing and refinement procedures is then applied in order to extract for each subevent (lasting itself more than `sd_thr=100` ns) a representative .pdb file (centroid of subevent.xtc).
All these "accepted" specific events are reported in: 

    $ cat pose_analysis_cgmd/acpt_sp_events_gromos_0.4.txt 
    #SP-serial #Replica LigID Event SubEvent CentroidTime[ps] Start[ns] End[ns] Duration[ns] #Coverage_pc #Avg-RMSD #devst
    1 rep1 COMB872 37 1 2.16e+06 2136 2343 207 0.9760 0.166745 0.085830
    2 rep1 COMB872 37 2 2.62e+06 2356 2959 603 0.9868 0.146382 0.093356
    3 rep1 COMB872 37 3 3.03e+06 2966 3076 110 0.9550 0.209906 0.149684
    ...   
    9 rep1 COMB872 141 1 8.926e+06 8442 9112 670 0.9821 0.153190 0.127478
    10 rep1 COMB872 141 2 9.823e+06 9754 9996 242 0.9383 0.277457 0.144962
    ...
    2501 rep35 COMB875 1 15 1.5821e+07 12236 16019 3783 0.9997 0.132315 0.075171
    2502 rep35 COMB875 1 16 1.6828e+07 16028 20050 4023 1.0000 0.133651 0.062431

This list contains information on the duration of the specific binding event,  the relative population coverage of the centroid and what is the average RMSD observed during the binding event.


The collection of representative poses is subjected to a new cluster analysis with a large cutoff to distinguish topographically distinct binding regions, denoted as "event domains". This allows to partition all the subevents in lists that are saved in 

     pose_analysis_cgmd/
	    └── frames_cntrd_${method[i]}_${thr[i]}
This list is a text file that can be inspected by the user:

    $ head -n 3 pose_analysis_cgmd/frames_cntrd_gromos_0.4/events_cluster_domain_1.txt 
    #CD #SP-Serial #Rep #LigID #Ev #Sub #CentroidTime[ps] #Start[ns] #End[ns] #Duration[ns]
    1 9 rep1 COMB872 141 1 8.926e+06 8442 9112 670 0.9821 0.153190 0.127478
    1 10 rep1 COMB872 141 2 9.823e+06 9754 9996 242 0.9383 0.277457 0.144962

Further classification of these poses will be performed by the second analysis script `tulipan_clust_assign.sh`.
The final step of this first part of the analysis is the computation of statistics on the binding probability.
In particular, the average probability of observing one ligand bound to the protein (regardless of its location) is reported in:

    $ cat pose_analysis_cgmd/summary_onelig.txt 
    #RepID #Pbound-Avg #Stdev
    rep1 0.871506 0.0885833
    rep2 0.883193 0.0979155
    ...
    rep35 0.829734 0.130685
    - - -
    Overall 0.881750 0.033640
## Part 2
The second part of the script identifies the most representative ligand poses within each binding pocket ("event domain"). All the results from this attribution process are save in a folder called `attrib_${method[i]}_${thr[i]}` . For each event domain, the list of its trajectories is copied as `unassigned.txt` and the first element is taken as the reference. The average RMSD distance of every other element of the list list from the reference is computed and if it is below the acceptance criterion `acc_thr=0.4` nm, the event is highlighted. All highlighted elements are classified as belonging to the same supracluster of the reference pose and are removed from the unassigned list. The procedure is iterated until all elements have been attributed to a supracluster.  
In the end, a table describing all supraclusters belonging to the same event domain is generated:

    $cat pose_analysis_cgmd/attrib_gromos_0.4/event_domain_1/table_supracl_dur_domain_1.txt
    #Domain #SupraclusterID #Events #AvgDuration #StDev #TopEventID #TopEventDuration
    1 1 24 161.88 48.88 2 248
    1 2 32 437.41 405.02 208 1476
    1 3 14 245.57 132.85 216 529
    1 4 1 187.0 0 28 187
    1 5 2 152.50 41.72 125 182
    ...
    1 13 1 250.0 0 189 250
    
For the supraclusters with more than `min_ev_restime=10` distinct events attributed, an histogram of the events' duration is created using Python. The histogram is then fitted according to a monoexponential decay $p(t)=A \cdot e^{-\frac{t}{\tau}}$ to estimate the residence time. Several criteria can be followed when building an histogram. The current implementation creates three different histograms according to the following rules:

 - Sturger's criterion
 - Doane's criterion
 - Freedman-Diaconis's criterion (FD) 
 
 Three different fitting attempts are therefore made for each representative pose and the results are stored in

    pose_analysis_cgmd/  
     └── attrib_gromos_0.4/  
          └── event_domain_1/  
              └── cntrd_1/  
                  └── fitresults.txt
For each method, the following values are reported

 - $A$
 - $\epsilon_{A}$
 - $\tau$
 - $\epsilon_{\tau}$
 - RMSE
 - Estimator from the fitted probability distribution

Finally, the script parses the three results and decides which should be reported in the final list, giving precedence to the results from the FD approach and falling-back to the other if FD failed (fitted parameters returned the fallback value of -1.0 as an error)

The final ranking of all the poses observed from the equilibrium simulations is reported in `pose_analysis_cgmd/attrib_gromos_0.4/table_domains_cvg_rstd_cum_RT.txt
`. This table describes:

 -  `Domain`: Where is the pose located;
 - `Supracluster`: which supracluster;
 - `Events` how many events are associated to the ligand pose;
 - `Duration`: sum of all the binding event durations;
 - `Cvg`:  relative coverage;
 - `Cvg-Cum` cumulative coverage;
 - `AvgSk` average duration (skipping top and bottom 5% of the data);
 - `DevSk` standard deviation over the trimmed average duration;
 - `Rule`: rule used for the histogram generation;
 - `Tau` residence time $\tau$;
 -  `ErrTau` uncertainty over the residence time $\epsilon_{\tau}$.

For example:

    $ head -n 5 table_domains_cvg_rstd_cum_RT.txt 
        #Domain #Supracluster #Events #Duration #Cvg #Cvg-Cum #AvgSk #DevSk #Rule #Tau #ErrTau
        2 5 214 150497 0.1177 0.1177 402.4 371.7 FD 175.2 12.3
        3 4 70 82626 0.0646 0.1823 978.3 1002.2 FD 553.0 105.7
        1 1 186 70081 0.0548 0.2371 320.4 203.0 FD 175.0 14.6
        1 2 105 36568 0.0286 0.2657 301.7 199.7 FD 157.1 9.5
        ...

From this final dataset, a selection of poses is extracted for visualization in VMD, depending on the value of `$extlimit`.
In particular, two `.tcl` files are generated, allowing to visualize just the PDBs of the top-poses or the PDBs with a concatenated trajectory of all the events belonging to the pose (beware of RAM requirements).

    pose_analysis_cgmd/attrib_gromos_0.4/load_extraction_selection_pdb.tcl  
    pose_analysis_cgmd/attrib_gromos_0.4/load_extraction_selection.tcl
To visualize the results, it is sufficient to run:

    vmd -e pose_analysis_cgmd/attrib_gromos_0.4/load_extraction_selection_pdb.tcl
    

 # Example files syntax
 ## Part 1
 ### listfile.txt
    /path/to/fitted_trajectory_1.xtc
    /path/to/fitted_trajectory_2.xtc
    /path/to/fitted_trajectory_3.xtc
    /path/to/fitted_trajectory_4.xtc
    /path/to/fitted_trajectory_5.xtc 
 ###  top_onelig.top
 
    #include "/path/to/itp/martini.itp"
    #include "/path/to/itp/protein.itp"
    #include "/path/to/itp/ligand.itp"
        
    [ system ]
    protein-1 ligand
    
    [ molecules ]
    PROT 1
    LIG 1
  ## Part 2
  ### blank_v3.tcl
  This files provide a template for the creation of a VMD-compatible `tcl` files that will load all the 10 most important binding modes and will set the colour of the ligand according to its ranking. The string `TEMP` will be dynamically updated by the script to include absolute path of `pose.pdb` `traj.xtc` and `color_ID`.
  A series of VMD representation is added to the molecule file in VMD. Please modify this template to reflect the structure and sequence numbering of the protein of interest.

     # Load a series of PDB files, each with its own XTC trajectory in VMD
    
    # Define list of PDB and XTC file pairs
    set files {
        TEMP
    }
    
    # Loop through each pair and load them in VMD
    foreach pair $files {
        set pdb_file [lindex $pair 0]  ;# Get PDB filename
        set xtc_file [lindex $pair 1]  ;# Get XTC filename
    	set pose_col [lindex $pair 2]  ;# Get VMD-Colour
        # Load the PDB file and get the molecule ID
        set molID [mol new $pdb_file type pdb waitfor all]
        # Etc [...]
            mol selection "resname LIGNAME"
            mol representation Licorice 0.3 30
            mol color ColorID $pose_col
            mol material AOShiny
            mol addrep top
        # Etc [...]
        # Load the corresponding XTC trajectory
        mol addfile $xtc_file type xtc waitfor all molid $molID
    
        puts "Loaded $pdb_file with trajectory $xtc_file"
    }
    
    puts "All structures and trajectories loaded successfully!"



> Written with [StackEdit](https://stackedit.io/).
