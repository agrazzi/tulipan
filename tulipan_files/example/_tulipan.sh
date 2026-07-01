#!/bin/bash
##
# Step 1 of TULIPAN analysis protocol
##
source "/usr/local/gromacs/bin/GMXRC"
source "/path/to/tulipan_files/functions/functions_tulipan.sh"

Help(){
	echo -e "-c file.gro\n-s topol.tpr\n-f trajlist.txt\n-l prot_onelig.tpr\n-p top_prot_onelig.top\n-m mdp.mdp\n"
}

#set -x 
# List of required flags
flag=("c_flag" "f_flag" "s_flag" "l_flag" "p_flag" "m_flag")
flagvar=("gro_conf" "listtraj" "tpr_traj" "prot_onelig_tpr" "topol_prot_onelig" "mdp")
extension=("gro" "txt" "tpr" "tpr" "top" "mdp")

# Initialize flags to 0
for k in "${flag[@]}"; do
    eval "$k=0"
done

# Parse command-line options
while getopts "c:s:f:l:p:m:h" option; do
    case "${option}" in
        c) gro_conf=../${OPTARG}; c_flag=1 ;;
        f) listtraj=../${OPTARG}; f_flag=1 ;;
        s) tpr_traj=../${OPTARG}; s_flag=1 ;;
        l) prot_onelig_tpr=../${OPTARG}; l_flag=1 ;;
        p) topol_prot_onelig=../${OPTARG}; p_flag=1 ;;
        m) mdp=../${OPTARG}; m_flag=1 ;;
        h) Help; exit ;;
        *) echo "Flag not recognized, please check with -h"; exit 1 ;;
    esac
done

### INPUT CHECK
for k in "${!flag[@]}"; do
    # Check if the flag is initialized (value is not 0)
    if [ "${!flag[$k]}" -eq "0" ]; then
        echo "${flag[$k]} has not been properly initialized"
        exit 1
    fi
    
    # Extract the last 3 characters of the flag value
    end="${!flagvar[$k]: -3}"
    # Check if the last 3 characters match the expected extension
    if [ "$end" != "${extension[$k]}" ]; then
        echo "Flag ${flag[$k]}: expected extension '${extension[$k]}', but got '$end'"
        exit 1
    fi

done
###

### Variable definition

    ## Workspace definition
        # Name of the output directory
        workdir="test_v5"
        # Limit on parallel tasks
        N=6

    ## Description of the biosystem
        # Ligand name
        ligname=COLC
        #ligand composition (how many beads)
        lig_comp=12

    ## Contacts approach parameters
        # Minimum number of protein-ligand contacts for annotation as "Bound" 
        cont_upcut=1 
        # Maximum number of protein-ligand contacts for annotation as "Unbound" 
        cont_lowcut=0
        #distance for contact calculation
        distcont=0.6
        #long-lived duration threshold (ns)
        ll_thr=100

    ## Events refinement protocol
        #gmx cluster method for initial identification of localized subevents
        method=("gromos")
        #corresponding array for clustering-cutoff (must have same number of elements as method[@])
        thr=("0.4")
        refine_method="gromos"
        refine_thr="0.4"
        #specific duration threshold (ns)
        sd_thr="100"
        #Criteria for assessing if a centroid of a subevent is truly representative of the entire subevent
        #If (avg_rmsd_centroid > $rmsd_thr) -> Warning!, not truly representative
        rmsd_thr="0.4"
        #If (population_centroid > $thr_subclust) -> Warning!, not truly representative
        thr_subclust="0.8"

###

# Create working folder
    if [ -d "$workdir" ];then
        rm -r $workdir
    fi
    echo "Creating work environment in folder $workdir"
    mkdir $workdir
    cd $workdir
    ### Check file existence
    for k in "${!flag[@]}"; do
    # Fetch the actual file path stored in the variable dynamically
    	file_path="${!flagvar[$k]}"
	if [ ! -f "$file_path" ]; then
		echo "Error: File '$file_path' (provided by flag ${flag[$k]}) does not exist."
	        exit 1
	fi
	# -------------------------------------------------------------
    done
    ###

    mkdir log/

# Begin time measurement
    tbegin=$(date +"%T.%3N")

# Print setting for bookkeeoing
    echo -e "gro file=$gro_conf\ntraj file=$listtraj\ntpr file=$tpr_traj\nligname=$ligname\nlig_comp=$lig_comp\ncont_upcut=$cont_upcut\ncont_lowcut=$cont_lowcut\nLL-duration_threshold=$ll_thr">parameters.txt
    echo -e "RefineMethod=$refine_method\nRefineMethodThreshold=$refine_thr\nSP-duration_threshold=$sd_thr\nMaxAvgRMSDdifference=$rmsd_thr\nWarningFractionalPopulationCoverage=$thr_subclust\n">>parameters.txt

#Sanity check on existence of trajectories
    check_traj_existence

#Parse initial structure file conf.gro to determine composition and ligand IDs 
    parse_gro

####################
# Contact analysis #
####################

#HEADER    
    echo "#LignameNR event# begin[ns] end[ns] duration[ns]">report_ll_events.txt
         
#Obtain contact timeseries
    echo "Contact timeseries section..."
    trj=0 #Current trajectory counter
    #Contact analysis
    while read traj; do #Iteration over all replicas	
        ((trj=trj+1))
        if [ "${traj:0:1}" = "/" ];then #already absolute path 
            xtc_traj="$traj"
        else 
            xtc_traj=$(realpath "../$traj") #generate absolute path
        fi         
        (trj_cnt_an "$xtc_traj" "$tpr_traj" "$trj")&
        ((submt=submt+1))
        echo -ne "\rSubmitted [$submt]/[$ntrj]"
        sleep 1
        while [ "$(jobs -r -p | wc -l)" -ge "$N" ]; do #Limit number of parallel tasks
            sleep 1
        done
        echo ""
    done<"${listtraj}" 
    echo ""
    echo -e "\rCheck if there are processes still running before continuing..."
    while [ "$(jobs -r -p | wc -l)" -gt "0" ]; do
            echo -ne "\rWaiting..."
            sleep 1
    done

#Clean-up unnecessary files (might be useful only for debug)
    rm mindist*xvg befehle log_mkndx.log log_mkndx_ev.log

#Generate overall report for long-lived events
    report_events "$ntrj" "$ligname" "lignr_list.dat"

echo -e  "\n#####\n"

####################
# Event refinement #
####################

#Setup for each method
for k in ${!method[@]}; do
    ##
    # Cluster analysis for refinment
    ##
    echo -e "\n#Considering: method: ${method[$k]} ctf: ${thr[$k]}#\n"            
    #Iteration over each replica folder
    for ((j=1;j<=ntrj;j++));do
        echo "Replica [$j]/[$ntrj]"
        #current ligand variable          
        clig=0
        #Iteration over every ligand
        while read lignr; do
            submt="0"
            ((clig=clig+1))
            echo -e "\tLigand [$clig]/[${tot_lig}]"
            lignmID="${ligname}${lignr}"            
            ligdir="rep${j}/${lignmID}"
            if [ -f "${ligdir}/temp_events_ll.txt" ];then 
                llevents=$(< "${ligdir}/temp_events_ll.txt" wc -l)
                #Iterate over each raw long-lived event
                while read line; do
                    ((submt=submt+1))
                    (refine_event "$line" "${ligdir}" "${lignmID}" "${submt}")&
                    sleep 0.25
                    echo -ne "\r\tSubmitted for refinement event [$submt]/[${llevents}]"

                    while [ "$(jobs -r -p | wc -l)" -ge "$N" ]; do
                        sleep 1
                    done

                done<"${ligdir}/temp_events_ll.txt"
                echo ""
            else
                echo "No LL Events for refinement ligand [$clig]/[$tot_lig]"
            fi
            
        done<"lignr_list.dat"
        echo ""
        echo -e "\rCheck if there are processes still running before changing trajectory..."
        while [ "$(jobs -r -p | wc -l)" -gt "0" ]; do
                echo -ne "\rWaiting..."
                sleep 0.5
        done
        echo ""
    done 
    echo ""
    echo -e "\rCheck if there are processes still running before summarising refinement results"
    while [ "$(jobs -r -p | wc -l)" -gt "0" ]; do
            echo -ne "\rWaiting..."
            sleep 0.5
    done

    echo -e  "\n#####\n"

    for ((j=1;j<=ntrj;j++));do
        while read lignr; do
            lignmID="${ligname}${lignr}"
            ligdir="rep${j}/${lignmID}"
            report="report_${method[$k]}.txt"
            if [ ! -f "$ligdir/${report}" ];then
                echo "#LL_Event Cutoff #Clusters" > "$ligdir/${report}"
            fi
            echo "#Event #Method #Cutoff #ClusterID #BeginTime[ns] #EndTime[ns] #Duration[ns]">${ligdir}/sp_event_report.txt
            echo "#Event #Method #Cutoff #Rejected">${ligdir}/sp_event_report_rctd.txt
			if [ -f "${ligdir}/temp_events_ll.txt" ];then
		        while read line;do
		            ev=$(echo "$line" | awk '{print $3}')
		            evdir="ev_${ev}"
		            echo -ne "${ev} ${thr[$k]} ">>"$ligdir/${report}"
		            log="ev_${ev}_log_${method[$k]}_${thr[$k]}.log"
		            grep "Found" $ligdir/${evdir}/$log | awk '{print $2}'>>"$ligdir/${report}"
		            
                    #Determining Refined Long-Lived"
		            new_cluster="clusters_${method[$k]}_${thr[$k]}_rlb.txt"
		            splist="sp_event_${lignmID}_ev_${ev}_${method[$k]}_${thr[$k]}.txt"
		            awk -v thr="$sd_thr" 'FNR==1{print $0} FNR>1{if($4>thr)print$0}' "$ligdir/${evdir}/${new_cluster}" > "${ligdir}/${evdir}/${splist}" 
		            nspev=$(< "${ligdir}/${evdir}/${splist}" wc -l)
		            #If at least one subevent lasted more than $sd_thr, add to list for extraction and further refinement analysis
		            if [ $nspev -gt "1" ];then
		                let nspev-- #Account for header
		                awk -v ev="${ev}" -v m="${method[$k]}" -v ctf="${thr[$k]}" 'FNR>1{ printf("%s %s %s %s\n",ev,m,ctf,$0)}' "${ligdir}/${evdir}/${splist}" >>${ligdir}/sp_event_report.txt 
		            else
		                echo "${ev} ${method[$k]} ${thr[$k]} y">>${ligdir}/sp_event_report_rctd.txt
		            fi
		        done<"${ligdir}/temp_events_ll.txt"
            fi
        done<"lignr_list.dat"
    done
done

echo -e "\n##########\n"

#Extracting subevents
echo -e "##Extracting subevents\n"
for k in ${!method[@]}; do
    (serial_spev["$k"]=0)
    echo "#SP-Serial #Replica #LigID #Event #SubEvent #Begin[ns] #End[ns] #Duration[ns]">"prerefine_sp_events_${method[$k]}_${thr[$k]}_serial.txt"
    for ((j=1;j<=ntrj;j++));do
        while read lignr; do
            lignmID="${ligname}${lignr}"
            ligdir="rep${j}/${lignmID}"
            #ev_to_ext=$(awk '{if ($4)c++}END{print c}' ${ligdir}/sp_event_report.txt )
            if [ -f "${ligdir}/temp_events_ll.txt" ];then
		        while read evline; do
		            ev=$(echo "$evline" | awk '{print $3}')
		            evdir="ev_${ev}"
		            ev_xtc="${lignmID}_ev_${ev}.xtc"
		            splist="sp_event_${lignmID}_ev_${ev}_${method[$k]}_${thr[$k]}.txt"
		            nsub=$(awk 'BEGIN{c=0}{if (substr($0,1,1)=="#")next; c++}END{print c}' "${ligdir}/${evdir}/${splist}")
		            if [ "$nsub" -ge "1" ];then
		                sub=0
		                subevdir="subevents_${method[$k]}_${thr[$k]}"
		                if [ ! -d "${ligdir}/${evdir}/${subevdir}/" ];then
		                    mkdir "${ligdir}/${evdir}/${subevdir}/"
		                fi
		                while read subline;do
		                    if [ "${subline::1}" = "#" ];then
		                        continue
		                    fi
		                    b=$(echo "$subline" | awk '{print $2*1000}')
		                    e=$(echo "$subline" | awk '{print $3*1000}')
		                    ((sub=sub+1))
		                    ((serial_spev["$k"]=serial_spev["$k"]+1))
		                    subdir="sub_${sub}"
		                    if [ ! -d "${ligdir}/${evdir}/${subevdir}/${subdir}/" ];then
		                        mkdir "${ligdir}/${evdir}/${subevdir}/${subdir}/"
		                    fi
		                    #echo "Extracting event [${ev}].[${sub}]"
		                    subevent="${lignmID}_ev_${ev}_sub_${sub}.xtc"
		                    gmx trjconv -f ${ligdir}/${evdir}/${ev_xtc} -b ${b} -e ${e} -o ${ligdir}/${evdir}/${subevdir}/${subdir}/${subevent} > /dev/null 2> log_subevent_extraction.log
		                    rm -rf \#*
		                    b=$(echo "$subline" | awk '{print $2}')
		                    e=$(echo "$subline" | awk '{print $3}')
		                    dur=$(echo "$subline" | awk '{print $4}')                               
		                    echo "${serial_spev[$k]} rep${j} ${lignmID} $ev $sub $b $e $dur ">>"prerefine_sp_events_${method[$k]}_${thr[$k]}_serial.txt" 
		                done<"${ligdir}/${evdir}/${splist}"
		            fi
		        done<"${ligdir}/temp_events_ll.txt"
            fi
        done<"lignr_list.dat"
    done
done

echo -e "\n##########\n"

#Subevent refining
    echo -e "##Refining every subevent\n"
    for k in ${!method[@]}; do
        nsubmt=$(awk '{if (substr($0,1,1)=="#")next; c++}END{print c}' "prerefine_sp_events_${method[$k]}_${thr[$k]}_serial.txt")
        submt=0
        subevdir="subevents_${method[$k]}_${thr[$k]}"
        echo "Method: [${method[$k]}]-[${thr[$k]}]"
        while read line; do
            if [ ${line::1} = "#" ];then
                continue
            fi
            ((submt=submt+1))
            echo -ne "\rSubmitting subevent [$submt]/[$nsubmt]..."
            (subev_clust_analysis "$line" "$subevdir" "$submt") &

            sleep 0.25
            while [ "$(jobs -r -p | wc -l)" -ge "$N" ]; do
                sleep 1
            done

        done<"prerefine_sp_events_${method[$k]}_${thr[$k]}_serial.txt"
        echo -e "\nAll submitted"
        echo -e "\rCheck if there are processes still running before changing method..."
        while [ "$(jobs -r -p | wc -l)" -gt "0" ]; do
                echo -ne "\rWaiting..."
                sleep 1
        done
        echo ""

    done

echo -e "\n##########\n"

#Check if subevent clustering has produced a timeseries for every subevent
    echo -e "###Sanity check on subevent clustering\n"
    for k in ${!method[@]}; do
        nsubmt=$(awk '{if (substr($0,1,1)=="#")next; c++}END{print c}' "prerefine_sp_events_${method[$k]}_${thr[$k]}_serial.txt")
        submt=0
        subevdir="subevents_${method[$k]}_${thr[$k]}"
        echo "Method: [${method[$k]}]-[${thr[$k]}]"
        while read line; do
            if [ ${line::1} = "#" ];then
                continue
            fi
            ((submt=submt+1))
            echo -ne "\rChecking existence of refinement tmsr for centroid [$submt]/[$nsubmt]"
            serial=$(echo "$line" | awk '{print $1}')
            repdir=$(echo "$line" | awk '{print $2}')
            lignmID=$(echo "$line" | awk '{print $3}')
            ev=$(echo "$line" | awk '{print $4}')
            sub=$(echo "$line" | awk '{print $5}')
            evdir="ev_${ev}"
            subdir="sub_${sub}"
            ligdir="${repdir}/${lignmID}"
            tmsr="ev_${ev}_sub_${sub}_tmsr_${refine_method}_${refine_thr}.xvg"
            if [ ! -f "$ligdir/${evdir}/${subevdir}/${subdir}/$tmsr" ] ;then
                echo "Missing $ligdir/${evdir}/${subevdir}/${subdir}/$tmsr"
                exit
            fi
            #clean log of succesfull refinement
            if [ -f "log_cluster_refine_${submt}.log"  ];then
                rm "log_cluster_refine_${submt}.log"
            fi
        done<"prerefine_sp_events_${method[$k]}_${thr[$k]}_serial.txt"
        echo ""
    done

echo -e "\n##########\n"

#Extraction of centroid pdb and tpr
    echo -e "###Extraction of .pdb and .tpr for every subevent centroid\n"
    for k in ${!method[@]}; do
        nsubmt=$(awk '{if (substr($0,1,1)=="#")next; c++}END{print c}' "prerefine_sp_events_${method[$k]}_${thr[$k]}_serial.txt")
        submt=0
        subevdir="subevents_${method[$k]}_${thr[$k]}"
        while read line; do
            if [ ${line::1} = "#" ];then
                continue
            fi
            ((submt=submt+1))
            echo -ne "\rSubmitting [$submt]/[$nsubmt]..."
            (centroid_generation "$line" "$subevdir" "$submt") &

            sleep 0.25
            while [ "$(jobs -r -p | wc -l)" -ge "$N" ]; do
                sleep 1
            done

        done<"prerefine_sp_events_${method[$k]}_${thr[$k]}_serial.txt"
        echo -e "\nAll submitted"
        echo -e "\rCheck if there are processes still running before changing method..."
        while [ "$(jobs -r -p | wc -l)" -gt "0" ]; do
                echo -ne "\rWaiting..."
                sleep 1
        done

    done

echo -e "\n##########\n"

#Sanity check on pdb and tpr and report on fractional coverage
    echo -e "###Sanity check on subevent .pdb .tpr and fractional coverage\n"
    for k in ${!method[@]}; do
        nsubmt=$(awk '{if (substr($0,1,1)=="#")next; c++}END{print c}' "prerefine_sp_events_${method[$k]}_${thr[$k]}_serial.txt")
        submt=0
        subevdir="subevents_${method[$k]}_${thr[$k]}"
        echo "Method: [${method[$k]}]-[${thr[$k]}]"
        echo "#Replica #LigID #ExtMet #ExtMetCtf #RefMet #RefMetCtf #Event #Subevent #FractionalCoverage" >"refine_warning_${method[$k]}_${thr[$k]}.txt"
        while read line; do
            if [ ${line::1} = "#" ];then
                continue
            fi
            ((submt=submt+1))
            echo -ne "\rChecking existence of .pdb and .tpr for centroid [$submt]/[$nsubmt]"
            centroid_existence $line
        done<"prerefine_sp_events_${method[$k]}_${thr[$k]}_serial.txt"
        echo ""
    done

echo -e "\n##########\n"

#Compute centroid consistency
    echo -e "##Assessing centroid consistency\n"
    for k in ${!method[@]}; do
        nsubmt=$(awk '{if (substr($0,1,1)=="#")next; c++}END{print c}' "prerefine_sp_events_${method[$k]}_${thr[$k]}_serial.txt")
        submt=0
        subevdir="subevents_${method[$k]}_${thr[$k]}"
        while read line; do
            if [ ${line::1} = "#" ];then
                continue
            fi
            ((submt=submt+1))
            echo -ne "\rSubmitting [$submt]/[$nsubmt]..."
            (centroid_consistency "$line" "$subevdir" "$submt") &

            sleep 0.25
            while [ "$(jobs -r -p | wc -l)" -ge "$N" ]; do
                sleep 1
            done

        done<"prerefine_sp_events_${method[$k]}_${thr[$k]}_serial.txt"
        echo -e "\nAll submitted"
        echo -e "\rCheck if there are processes still running before changing trajectory..."
        while [ "$(jobs -r -p | wc -l)" -gt "0" ]; do
                echo -ne "\rWaiting..."
                sleep 1
        done
        echo ""
    done

echo -e "\n##########\n"

#Check on centroid consistency
    echo -e "##Sanity check on centroid consistency\n"
    for k in ${!method[@]}; do
        nsubmt=$(awk '{if (substr($0,1,1)=="#")next; c++}END{print c}' "prerefine_sp_events_${method[$k]}_${thr[$k]}_serial.txt")
        submt=0
        subevdir="subevents_${method[$k]}_${thr[$k]}"
        echo "Method: [${method[$k]}]-[${thr[$k]}]"
        all_event_report="general_sp_events_${method[$k]}_${thr[$k]}.txt"
        acpt_event_report="acpt_sp_events_${method[$k]}_${thr[$k]}.txt"
        rfsd_event_report="rfsd_sp_events_${method[$k]}_${thr[$k]}.txt"
        cvg_report="coverage_report_${method[$k]}_${thr[$k]}.txt"
        echo "#SP-serial #Replica #LigID #Ev #Sub #CentroidTime[ps] #Coverage_pc #Avg-RMSD #devst">"$cvg_report"
        echo "#SP-serial #Replica LigID Event SubEvent CentroidTime[ps] Start[ns] End[ns] Duration[ns] #Coverage_pc #Avg-RMSD #devst" >"$acpt_event_report"
        echo "#SP-serial #Replica LigID Event SubEvent CentroidTime[ps] Start[ns] End[ns] Duration[ns] #Coverage_pc #Avg-RMSD #devst" >"$rfsd_event_report"
        echo "#SP-serial #Replica LigID Event SubEvent CentroidTime[ps] Start[ns] End[ns] Duration[ns] #Coverage_pc #Avg-RMSD #devst" >"$all_event_report"
        while read line; do
            if [ ${line::1} = "#" ];then
                continue
            fi
            ((submt=submt+1))
            echo -ne "\rChecking existence of rmsd-consistency for centroid [$submt]/[$nsubmt]"
            centroid_check $line           
        done<"prerefine_sp_events_${method[$k]}_${thr[$k]}_serial.txt"

        echo ""

    done

#Cluster centroids in topographically distinct events-domains
    echo -e "###Clustering centroids to partition in events-domains\n"
    for k in ${!method[@]}; do
        catlist=""
        befehle=""
        i=0
        evlist="acpt_sp_events_${method[$k]}_${thr[$k]}.txt"
        framedir="frames_cntrd_${method[$k]}_${thr[$k]}"
        if [ -d $framedir ]; then
            rm -r "$framedir"
        fi
        mkdir "$framedir"
        if [ ! -f "${framedir}/dict_frames.txt" ];then
            echo "Frame SP-serial">"${framedir}/dict_frames.txt"
        fi
        while read line; do
            if [ "${line::1}" = "#" ];then
                echo "#FrameID #SP-serial">"${framedir}/dict_frames.txt"
                continue
            else
                ((i=i+1000)) #for small i, gmx trjconv -t0 $i gives weird behaviour (e.g. -t0 23 and -t0 24 both produce frame with timestamp 24), expand number for better results
                #SP-serial #Replica LigID Event SubEvent CentroidTime[ps] Start[ns] End[ns] Duration[ns]
                pdblist["$i"]=$( echo "$line" | awk -v mth="${method[$k]}" -v thr="${thr[$k]}" '{printf("%s/%s/ev_%d/subevents_%s_%s/sub_%d/ev_%d_sub_%d_%s_%s_lc.pdb",$2,$3,$4,mth,thr,$5,$4,$5,mth,thr)}')
                if [ ! -f "${pdblist[$i]}" ];then
                    echo "Could not find ${pdblist[$i]}"
                    exit
                else
                    id=$(echo "$line" | awk '{print $1}')
                    out="centroid_${id}.xtc"
                    #Save as time the id of the cluster
                    gmx trjconv -f ${pdblist[$i]} -o "${framedir}/${out}" -t0 "$i" &>> "${framedir}/log_trjconv.log"
                    echo "$i $id">>"${framedir}/dict_frames.txt"
                    catlist+="${framedir}/${out} "
                    befehle+="c\n"
                fi
            fi
        done<"${evlist}"
    
        echo -e "$befehle" | gmx trjcat -f $(echo "$catlist") -o "${framedir}/cat.xtc" -cat  > /dev/null 2> ${framedir}/log_trjcat.log
        gmx check -f "${framedir}/cat.xtc" >&"${framedir}/log_gmxcheck.log"  
        
        if [ ! -d "${framedir}/centroids_xtc" ];then
            mkdir "${framedir}/centroids_xtc"
        fi
        mv ${framedir}/centroid_*.xtc "${framedir}/centroids_xtc"
        echo -e "$ligname" | gmx cluster -method gromos -cutoff 0.8 -f ${framedir}/cat.xtc -clid "${framedir}/tmsr_cntrd_domain.xvg" -s ${prot_onelig_tpr} -nofit -g "${framedir}/gmx_cluster.log" -dist "${framedir}/gmx_cluster_rmsdist.xvg"> /dev/null 2> ${framedir}/log_gmxcluster.log
        
        awk -v outdir="${framedir}" '	
            ARGIND==1 {if (substr($0,1,1)=="#" || substr($0,1,1)=="@") next 
                        l++
                        line[$1]=$0
                        }
            ARGIND==2{ if (substr($0,1,1)=="#" || substr($0,1,1)=="@") next
                            dict[$1]=$2
                    }
            ARGIND==3{ if (substr($0,1,1)=="#" || substr($0,1,1)=="@") next
                            if (gen[$2]=="")
                        {	out=sprintf("%s/events_cluster_domain_%d.txt",outdir,$2)
                                print "#CD #SP-Serial #Rep #LigID #Ev #Sub #CentroidTime[ps] #Start[ns] #End[ns] #Duration[ns]" > out
                            gen[$2]=1
                        }
                    }
            ARGIND==4 { if (substr($0,1,1)=="#" || substr($0,1,1)=="@") next
                            out=sprintf("%s/events_cluster_domain_%d.txt",outdir,$2)
                            printf("%d %s\n",$2,line[dict[int($1)]])>>out
                            k++
                        }
            END     {   outreport=sprintf("%s/awk_error.err",outdir)
                        if (l!=k)
                        {    printf("Mismatch between accepted events read in input [%d] and events printed to cluster-domains [%d]\n",l,k)>outreport
                            exit 1}
                        else
                            exit 0}' ${evlist} "${framedir}/dict_frames.txt" "${framedir}/tmsr_cntrd_domain.xvg" "${framedir}/tmsr_cntrd_domain.xvg"
        if [ "$(echo $?)" -eq "1" ]; then
            echo "AWK Attribution of centroid to centroid domain [${method[$k]}]-[${thr[$k]}] has returned an error"
            exit
        fi
    done

echo -e "\n######\n"

#Analysis of P-bound statistics
    echo "Averaging Pbound"
    binding_statistics $ntrj $tot_lig

#Final clean-up
    mv *.log log/
    echo "END"
    #END time measurement
    tend=$(date +"%T.%3N")
    echo -e "Script began at: $tbegin\nScript ended at: $tend"
    echo -e "Script began at: $tbegin\nScript ended at: $tend">>parameters.txt
