#!/bin/bash
function check_traj_existence {
    ntrj=0   
    while read traj; do  
        if [ ${traj:0:1} = "/" ];then #already absolute path 
            if [ ! -f "$traj" ];then
                echo "$traj could not be found"
                exit 1
            fi
        else 
            if [ ! -f "../$traj" ];then #relative path starting from launch directory, e.g. ../fit.xtc
                echo "$traj could not be found"
                exit 1
            fi
        fi
        let "ntrj+=1"
    done<"${listtraj}"
}

function parse_gro {
#Parse initial structure file conf.gro to determine composition and ligand IDs 
    #Extract all ligands present in simulation box 
        cp ${gro_conf} temp_gro.txt 
        sed -i '1,2d' temp_gro.txt
        grep "$ligname" temp_gro.txt > ligpos.txt
        rm temp_gro.txt			
    #Compute number of ligands
        ((tot_lig=$(< ligpos.txt wc -l)/lig_comp))
        echo "Ligand molecules contained in system: $tot_lig"   
    
    #Read ligpos.txt and obtain list of ligand IDs 
        c=0
        while read -r line; do
            candidate=$(echo "$line" | awk '{print $1}')
            if [ $c -eq "0" ];then
                ((c=lig_comp))
            fi
            if [ $c -eq $lig_comp ]; then
                echo "$candidate">>liglist.dat
                lignr=$(echo "$candidate" | grep -o -E '[0-9]+')
                echo "$lignr">>lignr_list.dat
            fi
            let "c=c-1"
        done<ligpos.txt

	#Sanity check
        ((check=$(wc -l < liglist.dat )-tot_lig))
        if [[ check -ne "0" ]];then
            echo "Expecting $tot_lig ligands, but liglist.dat contains $(wc -l < liglist.dat ) entries: mismatch in ligand identification."
            exit 1
        fi
    
    #generate first index.ndx
        echo -e "q \n" | gmx make_ndx -f ${gro_conf} >&log.log
    
    #Add ligand to index
        while read -r line; do
            candidate=$(echo "$line" | awk '{print $1}') 
            awk -v lig="$ligname" -v resid_target="$candidate" '
            BEGIN {
                print "[ " lig "_r_" resid_target " ]"
            }

            NR > 2 {  # skip header + atom count line

                resid  = substr($0,1,5)
                resname = substr($0,6,5)
                atomnr = substr($0,16,5)

                gsub(/^[ \t]+|[ \t]+$/, "", resid)
                gsub(/^[ \t]+|[ \t]+$/, "", resname)
                gsub(/^[ \t]+|[ \t]+$/, "", atomnr)

                if (resname == lig && resid == resid_target) {
                    printf "%s ", atomnr
                }
            }

            END {
                printf "\n"
            }
            ' "${gro_conf}" >> index.ndx      
        done<lignr_list.dat
        
        # Create all Protein_OneLigand groups for future extraction
        echo "Creating Prot_onelig_selection"
        
        awk 'BEGIN{l=0}{if ($1=="[") {print l " " $2; l++}}' index.ndx > numgroup.txt
        grep "${ligname}_r_" numgroup.txt | awk '{print $1 " " $2}' >numgroup_lig.dat
        
        while read -r line; do
            targetmerge=$(echo "$line" | awk '{print $1}')
            echo -e "1 | $targetmerge \n q \n" | gmx make_ndx -f ${gro_conf} -n index.ndx > /dev/null 2> log_mkndx.log  
            rm \#*
        done<numgroup_lig.dat

        #Needed by gmx mindist contact analysis;  
        awk -v lig="$ligname" 'BEGIN{print "Protein "}{printf("%s_r_%d\n",lig,$1) }' lignr_list.dat >befehle

        rm numgroup.txt numgroup_lig.dat ligpos.txt
}

function trj_cnt_an {
    local xtc_traj=$1
    local tpr_traj=$2
    local trj=$3
    local tgtdir="rep${trj}"
    mkdir ${tgtdir}

    cat befehle |gmx mindist -f "${xtc_traj}" -n index.ndx -d $distcont -dt 1000 -s ${tpr_traj} -on ${tgtdir}/general_cnt_tmsr.xvg -ng $tot_lig -od "mindist_${trj}.xvg" > /dev/null 2> log_mindist.log              
    #Clean output
    #header of timeseries will be: time[ns] #contacts
    local i=0
    local k=0       
    for ((i=1; i <= tot_lig ; i=i+1 )); do
        ((k=i+1))
        #echo "Generating tmsr_${i}.txt"            
        awk -v j=$k '{if (substr($1,1,1) != "#" && substr($1,1,1) != "@"){print $1/1000 " " $j}}' ${tgtdir}/general_cnt_tmsr.xvg >> ${tgtdir}/tmsr_${i}.txt
        
    done
    local line
    i=0
    while read -r line; do #Iteration over all ligands of a replica
        ((i=i+1))
        #echo "$line"
        #echo -e "Replica [${trj}]/[${ntrj}]\tWorking on ligand [$clig]/[$tot_lig]"
        local candidate=$(echo "$line" | awk '{print $1}')
        local lignmID="${ligname}${candidate}"
        local ligdir="${tgtdir}/${lignmID}"
        mkdir $ligdir  
        local targetsel="${ligname}_r_${candidate}"
        #contact timeseries
        local cont_timsr=${lignmID}_cont_tmsr
        mv ${tgtdir}/tmsr_${i}.txt "$ligdir/${cont_timsr}.txt"                
        echo "time[ns] contacts status">"$ligdir/${cont_timsr}_b_u.txt"
        awk -v up=$cont_upcut -v lc=$cont_lowcut '
            BEGIN{status=0}
            {	c=$2;
                if (status==0)
                {	if (c>=up)
                        {	status=1;
                            print $0 " " 1;
                        }
                    else
                        print $0 " " 0
                }
                else 
                {	if (c<=lc) #unbinding
                    {	status=0;
                        print $0 " " 0
                    }
                    else 
                    {	print $0 " " 1;											
                    }
                }                    
            }' $ligdir/${cont_timsr}.txt >>$ligdir/${cont_timsr}_b_u.txt
        #  
        #echo "Will target $ligdir/${cont_timsr}_b_u.txt for individual contact analysis"
        awk 'BEGIN{l=0;c=0}NR>1{l++; if($3)c++;}END{print l " " c " " l-c}' $ligdir/${cont_timsr}_b_u.txt >> $ligdir/report_${lignmID}.txt
        
        #echo "Contact duration"
        echo "Replica LignameNR event# begin[ns] end[ns] duration[ns]">$ligdir/events.txt
        
        #counter begin end duration
        tail -n +2 $ligdir/${cont_timsr}_b_u.txt | awk -v rep="$tgtdir" -v lign="${lignmID}" 'BEGIN{l=0;s=0;c=0;}
        {if($3)
        {   if (s==0)
            {s=1;
            c++;            
            b=$1;}           
        }
        else
            if (s) #unbinding
            {   e=$1;
                printf("%s %s %d %d %d %d\n",rep,lign,c,b,e,e-b)
                s=0;
            }
        }END{if (s)
                {e=$1;
                printf("%s %s %d %d %d %d\n",rep,lign,c,b,e,e-b)}}' >> $ligdir/events.txt
        
        #echo "Extracting longlived binding events"
        echo "Replica LignameNR event# begin[ns] end[ns] duration[ns]">>$ligdir/events_ll.txt
        awk -v thr=$ll_thr 'BEGIN{l=0}{if (l){if ($6>thr) print $0} l++;}' $ligdir/events.txt >> $ligdir/events_ll.txt
        local c=$(awk 'BEGIN{s=0}NR>1{s++}END{print s}' $ligdir/events_ll.txt )
        #Extract LL events
        if [ "$c" -ne "0" ]; then
            tail -n +2  $ligdir/events_ll.txt > $ligdir/temp_events_ll.txt
            local input
            while read input; do
                local ev=$(echo "$input" | awk '{print $3}')
                local b=$(echo "$input" | awk '{print $4}')
                local e=$(echo "$input" | awk '{print $5}')
                local evdir="ev_$ev"
                mkdir "${ligdir}/${evdir}"
                b=$(echo "$b" | awk '{print $1*1000}')
                e=$(echo "$e" | awk '{print $1*1000}')
                output="Protein_${ligname}_r_${candidate}"
                echo "$output" | gmx trjconv -f "${xtc_traj}" -s $tpr_traj -n index.ndx  -b $b -e $e -o ${ligdir}/${evdir}/${lignmID}_ev_${ev}.xtc > /dev/null 2> log_extract.log
                #echo "Extract first frame"
                echo "$output" | gmx trjconv -f "${xtc_traj}" -s $tpr_traj -n index.ndx  -b $b -e $b -o ${ligdir}/${evdir}/${lignmID}_ev_${ev}.pdb > /dev/null 2> log_extract.log
                prot_lig_ndx="prot_${lignmID}.ndx"
                if [ ! -f  $ligdir/${prot_lig_ndx} ]; then
                    echo -e "r ${candidate}\nq" | gmx make_ndx -f ${ligdir}/${evdir}/${lignmID}_ev_${ev}.pdb -o $ligdir/${prot_lig_ndx} > /dev/null 2> log_mkndx_ev.log
                    sed -i "s/r_${candidate}/${lignmID}/g" $ligdir/${prot_lig_ndx}
                fi
            done<"$ligdir/temp_events_ll.txt"

        fi


    done<lignr_list.dat 
}

function report_events {
    local ntrj=$1
    local ligname=$2
    local list=$3
    local i=0
    for ((i=1;i<=ntrj;i++));do
        local tgtdir="rep${i}"
        while read line; do
            local candidate=$(echo "$line" | awk '{print $1}')
            local lignmID="${ligname}${candidate}"
            local ligdir="${tgtdir}/${lignmID}"
            local c=$(awk 'BEGIN{s=0}NR>1{s++}END{print s}' $ligdir/events_ll.txt )
            if [ "$c" -eq "0" ]; then
                echo "$tgtdir ${lignmID}" >>warning_ll.txt
            else
                tail -n +2  $ligdir/events_ll.txt > $ligdir/temp_events_ll.txt
            fi
            tail -n +2  $ligdir/events.txt >>report_events.txt
            tail -n +2  $ligdir/events_ll.txt >>report_ll_events.txt

        done<"$list"
    done
}

function refine_event {
    local line=$1
    local ligdir=$2
    local lignmID=$3
    local submt=$4
    local ev=$(echo "$line" | awk '{print $3}')
    local evdir="ev_$ev"
    local ev_xtc="${lignmID}_ev_${ev}.xtc"
    local log="ev_${ev}_log_${method[$k]}_${thr[$k]}.log"
    local tmsr="ev_${ev}_tmsr_${method[$k]}_${thr[$k]}.xvg"
    local prot_lig_ndx="prot_${lignmID}.ndx"
    echo -e "${lignmID}" | gmx cluster -method ${method[$k]} -f ${ligdir}/${evdir}/${ev_xtc} -s ${prot_onelig_tpr} -n ${ligdir}/${prot_lig_ndx} -g $ligdir/${evdir}/$log -clid $ligdir/${evdir}/$tmsr -nofit -cutoff ${thr[$k]} -o "$ligdir/${evdir}/clust.xpm" -dist "$ligdir/${evdir}/clust.xvg"  > /dev/null 2> log_cluster_ev_${submt}.log
    rm "$ligdir/${evdir}/clust.xpm" "$ligdir/${evdir}/clust.xvg"
    #echo "Reading cluster log to determine centroid positions in time"
    local cluster_dict="cluster_dict_ev_${ev}_${method[$k]}_${thr[$k]}.txt"
    awk 'BEGIN{printf("#ClusterID #TotFrames #Start\n")}{
        #printf("$1 is: %s\t",$1)
        if ( $1 == "cl."){ s=1; next}
        if (s != 1)next
        if (s == 1)
            {	if ($1 == "|" )
                {	#print "Found | , skipping"
                    next
                }
                else
                {	if ($3 >1) # 
                        printf("%d %d %s\n",$1,$3,$6)
                    else
                        printf("%d %d %s\n",$1,$3,$5)
                }
            }} ' $ligdir/${evdir}/$log > $ligdir/${evdir}/${cluster_dict}
    #echo "Cluster duration"
    local cluster_report="ev_${ev}_clusters_${method[$k]}_${thr[$k]}.txt"
    awk 'BEGIN{c=-1}{
        if(substr($0,1,1)=="#" || substr($0,1,1)=="@") next
        t=$1/1000
        if (c==-1)
            {c=$2
            bt=t}
        else
        { if (c!=$2)
            {et=t
            printf("%d %d %d %d\n",c,bt,et,(et-bt))
            bt=t
            c=$2}
        }
        }END{et=t; printf("%d %d %d %d\n",c,bt,et,(et-bt)+1)}' $ligdir/${evdir}/${tmsr} >"$ligdir/${evdir}/${cluster_report}" 

    #echo "Smoothing"
    local temp_cluster="temp_clusters_${method[$k]}_${thr[$k]}_rlb.txt"
    local new_cluster="clusters_${method[$k]}_${thr[$k]}_rlb.txt"
    local new_tmsr="tmsr_${method[$k]}_${thr[$k]}_rlb.txt"
    awk 'BEGIN{c=-1;i=0;printf("ClusterID BeginTime[ns] EndTime[ns] Duration[ns]\n")}
        {	if (c==-1)
            {c=$1
            buffer=1+int(log($4))
            printf("%s\n",$0)}
        else
            if (c==$1)
            {printf("%s\n",$0)}
            else #Challenge section
            { if ($4<=buffer) #New Sections has duration less than buffer
                {	printf("%s %d %d %d\n",c,$2,$3,$4)
                }
            else
                {printf("%s\n",$0)
                c=$1
                buffer=int(log($4))}
            }
            
            }'  "$ligdir/${evdir}/${cluster_report}" > "$ligdir/${evdir}/${temp_cluster}"
        
    #echo "Re-Formatting"
    awk 'BEGIN{printf("#ClusterID #BeginTime[ns] #EndTime[ns] #Duration[ns]\n")}
        NR==2{c=$1;bt=$2;et=$3;}
        NR>2{	if ($1==c)#
                    et=$3
                else
                {	printf("%s %d %d %d\n",c,bt,et,(et-bt))
                    c=$1; 
                    bt=$2;
                    et=$3	
                }
        }END{printf("%s %d %d %d\n",c,bt,et,(et-bt)+1)}' "$ligdir/${evdir}/${temp_cluster}" > "$ligdir/${evdir}/${new_cluster}"
    rm "$ligdir/${evdir}/${temp_cluster}"
    awk 'BEGIN{printf("Time[ns] Cluster\n")} NR>1{c=$1; bt=$2; et=$3; for(t=bt;t<=et;t++) printf("%d %s\n",t,c) } ' "$ligdir/${evdir}/${new_cluster}" > "$ligdir/${evdir}/${new_tmsr}"
    #Sanity check
    local old_dur=$(awk '{sum+=$4}END{print sum}' "$ligdir/${evdir}/${cluster_report}" ) #No header
    local new_dur=$(awk 'NR>1{sum+=$4}END{print sum}' "$ligdir/${evdir}/${new_cluster}" )
    if [ "$old_dur" -ne "$new_dur" ];then
        echo -e  "Mismatch [ev_${ev}] \n Duration "$ligdir/${evdir}/${cluster_report}" : ${old_dur}\nDuration "$ligdir/${evdir}/${new_cluster}" : ${new_dur}\n">>mismatch.txt
        exit 1
    fi
    
}
#line from prerefine
function subev_clust_analysis {
    local line="$1"
    local subevdir="$2"
    local submt="$3"
    local serial=$(echo "$line" | awk '{print $1}')
    local tgtdir=$(echo "$line" | awk '{print $2}')
    local lignmID=$(echo "$line" | awk '{print $3}')
    local ev=$(echo "$line" | awk '{print $4}')
    local sub=$(echo "$line" | awk '{print $5}')
    local ligdir="${tgtdir}/${lignmID}"
    local evdir="ev_${ev}"
    local subdir="sub_${sub}"
    local subevent="${lignmID}_ev_${ev}_sub_${sub}.xtc"
    local prot_lig_ndx="prot_${lignmID}.ndx"
    local log="ev_${ev}_sub_${sub}_log_${refine_method}_${refine_thr}.log"
    local tmsr="ev_${ev}_sub_${sub}_tmsr_${refine_method}_${refine_thr}.xvg"


    #apply cluster analysis to subevent
    echo -e "${lignmID}" | gmx cluster -method ${refine_method} -f ${ligdir}/${evdir}/${subevdir}/${subdir}/${subevent} -s ${prot_onelig_tpr} -n ${ligdir}/${prot_lig_ndx} -g $ligdir/${evdir}/${subevdir}/${subdir}/$log -clid $ligdir/${evdir}/${subevdir}/${subdir}/$tmsr -nofit -cutoff ${refine_thr}  -o "$ligdir/${evdir}/${subevdir}/${subdir}/clust.xpm" -dist "$ligdir/${evdir}/${subevdir}/${subdir}/clust.xvg"> /dev/null 2> log_cluster_refine_${submt}.log
    rm "$ligdir/${evdir}/${subevdir}/${subdir}/clust.xpm" "$ligdir/${evdir}/${subevdir}/${subdir}/clust.xvg"
    local cluster_dict="cluster_dict_ev_${ev}_sub_${sub}_ref_${refine_method}_${refine_thr}.txt"
    awk 'BEGIN{printf("#ClusterID #TotFrames #Start\n")}
        {   if ( $1 == "cl."){ s=1; next}
            if (s != 1)next
            if (s == 1)
                {	if ($1 == "|" )
                    {	next
                    }
                    else
                    {	if ($3 >1) # 
                            printf("%d %d %s\n",$1,$3,$6)
                        else
                            printf("%d %d %s\n",$1,$3,$5)
                    }
                }
        } ' $ligdir/${evdir}/${subevdir}/${subdir}/$log > $ligdir/${evdir}/${subevdir}/${subdir}/temp_dict.txt

    awk '   NR==1{l[FNR]=$0}
            NR>1{l[FNR]=$0; pop[FNR]=$2;sum+=$2}
            END {for(i=1;i<=FNR;i++)
                    if (i==1)
                        printf("%s Pop%s\n",l[i],"%") 
                    else
                        printf("%s %.4f\n",l[i],pop[i]/sum)}' $ligdir/${evdir}/${subevdir}/${subdir}/temp_dict.txt > $ligdir/${evdir}/${subevdir}/${subdir}/${cluster_dict}
}

function centroid_generation {
    local line="$1"
    local subevdir="$2"
    local submt="$3"
    local serial=$(echo "$line" | awk '{print $1}')
    local tgtdir=$(echo "$line" | awk '{print $2}')
    local lignmID=$(echo "$line" | awk '{print $3}')
    local ev=$(echo "$line" | awk '{print $4}')
    local sub=$(echo "$line" | awk '{print $5}')
    local ligdir="${tgtdir}/${lignmID}"
    local evdir="ev_${ev}"
    local subdir="sub_${sub}"
    local subevent="${lignmID}_ev_${ev}_sub_${sub}.xtc"
    local prot_lig_ndx="prot_${lignmID}.ndx"
    #Time reference largest cluster
    local cluster_dict="cluster_dict_ev_${ev}_sub_${sub}_ref_${refine_method}_${refine_thr}.txt"
    local ref_lcl=$(echo "" | awk 'NR==2{print $3}' $ligdir/${evdir}/${subevdir}/${subdir}/${cluster_dict})  
    local subev_pdb="ev_${ev}_sub_${sub}_${method[$k]}_${thr[$k]}_lc.pdb"
    local subev_tpr="ev_${ev}_sub_${sub}_${method[$k]}_${thr[$k]}_lc.tpr"

    echo -e "System\n" | gmx trjconv -f ${ligdir}/${evdir}/${subevdir}/${subdir}/${subevent} -s ${prot_onelig_tpr}  -b "${ref_lcl}" -e "${ref_lcl}" -n ${ligdir}/${prot_lig_ndx} -o ${ligdir}/${evdir}/${subevdir}/${subdir}/${subev_pdb} > /dev/null 2> log_extraction_subcluster_${submt}.log
    #Extraction
    gmx grompp -f ${mdp} -p ${topol_prot_onelig} -c "${ligdir}/${evdir}/${subevdir}/${subdir}/${subev_pdb}" -o "${ligdir}/${evdir}/${subevdir}/${subdir}/${subev_tpr}" > /dev/null 2> log_grompp.log
    if [ ! -f "${ligdir}/${evdir}/${subevdir}/${subdir}/${subev_tpr}" ]; then
        echo "[ERROR]:  ${cntrd_dir}/${subev_tpr} not generated"
        exit
    fi
}

function centroid_consistency {
    local line="$1"
    local subevdir="$2"
    local submt="$3"
    local serial=$(echo "$line" | awk '{print $1}')
    local tgtdir=$(echo "$line" | awk '{print $2}')
    local lignmID=$(echo "$line" | awk '{print $3}')
    local ev=$(echo "$line" | awk '{print $4}')
    local sub=$(echo "$line" | awk '{print $5}')
    local ligdir="${tgtdir}/${lignmID}"
    local evdir="ev_${ev}"
    local subdir="sub_${sub}"
    local subevent="${lignmID}_ev_${ev}_sub_${sub}.xtc"
    local prot_lig_ndx="prot_${lignmID}.ndx"
    
    local cluster_dict="cluster_dict_ev_${ev}_sub_${sub}_ref_${refine_method}_${refine_thr}.txt"
    #population largest cluster
    local pop_lcl=$(echo "" | awk 'NR==2{print $2}' $ligdir/${evdir}/${subevdir}/${subdir}/${cluster_dict})
    #Coverage % largest cluster
    local cov_lcl=$(echo "" | awk 'NR==2{print $4}' $ligdir/${evdir}/${subevdir}/${subdir}/${cluster_dict})
    #Time reference largest cluster 
    local ref_lcl=$(echo "" | awk 'NR==2{print $3}' $ligdir/${evdir}/${subevdir}/${subdir}/${cluster_dict})

    local subev_pdb="ev_${ev}_sub_${sub}_${method[$k]}_${thr[$k]}_lc.pdb"
    local subev_tpr="ev_${ev}_sub_${sub}_${method[$k]}_${thr[$k]}_lc.tpr"
    #Own consistency
    local cntrd_dir="${ligdir}/${evdir}/${subevdir}/${subdir}"
    outxvg="rmsd_cntrd_${serial}_vs_${tgtdir}_${lignmID}_ev${ev}_sub${sub}.xvg"
    echo "${lignmID}" | gmx rms -f "${cntrd_dir}/${subevent}" -s "${cntrd_dir}/${subev_tpr}" -o "${cntrd_dir}/${outxvg}" -n  "${ligdir}/${prot_lig_ndx}" -fit none > /dev/null 2> log_rms_${submt}.log
    avgtxt="rmsd_cntrd_${serial}_vs_${tgtdir}_${lignmID}_ev${ev}_sub${sub}.txt"

    #avg stdev frames
    awk -v rep="$tgtdir" -v lig="$lignmID" -v ev="$ev" -v subev="$sub" -v sr="${serial}" '
    ARGIND==1{if (substr($0,1,1)=="#" || substr($0,1,1)=="@") next; c++;sum+=$2} 
    ARGIND==2 && FNR==1 {avg=sum/c}
    ARGIND==2{if (substr($0,1,1)=="#" || substr($0,1,1)=="@") next; sumq+=((avg-$2)**2)}
    END{printf("%s %s %d %d %f %f %d %d\n",rep,lig,ev,subev,avg,sqrt(sumq/(c-1)),c,sr)}' ${cntrd_dir}/${outxvg} ${cntrd_dir}/${outxvg} > "${cntrd_dir}/${avgtxt}"

}
function centroid_existence {
    serial=$1
    repdir=$2
    lignmID=$3
    ev=$4
    sub=$5
    evdir="ev_${ev}"
    subdir="sub_${sub}"
    ligdir="${repdir}/${lignmID}"
    cntrd_dir="${repdir}/${lignmID}/ev_${ev}/subevents_${method[$k]}_${thr[$k]}/sub_${sub}"
    subev_pdb="ev_${ev}_sub_${sub}_${method[$k]}_${thr[$k]}_lc.pdb"
    subev_tpr="ev_${ev}_sub_${sub}_${method[$k]}_${thr[$k]}_lc.tpr"
    if [ ! -f "${cntrd_dir}/${subev_pdb}" ] ;then
        echo "Missing ${cntrd_dir}/${subev_pdb}"
        exit 1
    fi
    #clean log of succesfull extraction
    if [ -f "log_extraction_subcluster_${submt}.log" ];then
        rm "log_extraction_subcluster_${submt}.log"
    fi
    if [ ! -f "${cntrd_dir}/${subev_tpr}" ]; then
        echo "Missing ${cntrd_dir}/${subev_tpr}"
        exit 1
    fi
            #clean log of succesfull extraction
    if [ -f "log_grompp_${submt}.log" ];then
        rm "log_grompp_${submt}.log" 
    fi        
    #Coverage % largest cluster
    cluster_dict="cluster_dict_ev_${ev}_sub_${sub}_ref_${refine_method}_${refine_thr}.txt"
    cov_lcl=$(echo "" | awk 'NR==2{print $4}' $ligdir/${evdir}/${subevdir}/${subdir}/${cluster_dict}) 
    check=$( awk -v thr="$thr_subclust" 'NR==2{pop=$4;exit}END{if(pop<thr)print "0"; else print "1"}' "$ligdir/${evdir}/${subevdir}/${subdir}/${cluster_dict}")
    if [ $check -eq "0" ]; then
        #echo "Warning: Fractional Coverage [${cov_lcl}] lower than current threshold: [${thr_subclust}]" 
        #CoverageCluster                           
        echo "${repdir} ${lignmID} ${method[$k]} ${thr[$k]} ${refine_method} ${refine_thr} ${ev} ${sub} ${cov_lcl}" >>refine_warning_${method[$k]}_${thr[$k]}.txt
    fi 
}

function centroid_check {
    serial=$1 
    repdir=$2 
    lignmID=$3
    ev=$4 
    sub=$5 
    b=$6 
    e=$7 
    dur=$8 
    evdir="ev_${ev}"
    subdir="sub_${sub}"
    ligdir="${repdir}/${lignmID}"
    cntrd_dir="${repdir}/${lignmID}/ev_${ev}/subevents_${method[$k]}_${thr[$k]}/sub_${sub}"
    outxvg="rmsd_cntrd_${serial}_vs_${repdir}_${lignmID}_ev${ev}_sub${sub}.xvg"
    avgtxt="rmsd_cntrd_${serial}_vs_${repdir}_${lignmID}_ev${ev}_sub${sub}.txt"
    
    #Time reference largest cluster
    cluster_dict="cluster_dict_ev_${ev}_sub_${sub}_ref_${refine_method}_${refine_thr}.txt"
    ref_lcl=$(echo "" | awk 'NR==2{print $3}' $ligdir/${evdir}/${subevdir}/${subdir}/${cluster_dict})
    #Fractional Coverage largest cluster
    cov_lcl=$(echo "" | awk 'NR==2{print $4}' $ligdir/${evdir}/${subevdir}/${subdir}/${cluster_dict})   
    
    if [ ! -f "${cntrd_dir}/${outxvg}" ]; then
        echo -e "\nMissing ${cntrd_dir}/${outxvg}"
        exit
    fi

    if [ ! -f "${cntrd_dir}/${avgtxt}" ] ;then
        echo -e "\nMissing ${cntrd_dir}/${avgtxt}"
        exit
    fi
            #clean log of succesful extraction
    if [ -f "${cntrd_dir}/${avgtxt}" ];then 
        if [ -f "log_rms_${submt}.log" ]; then
            rm "log_rms_${submt}.log"
        fi
    fi
    
    #Read Avg and assess consistency
    avgrmsd=$(awk '{print $5}' "${cntrd_dir}/${avgtxt}")
    devrmsd=$(awk '{print $6}' "${cntrd_dir}/${avgtxt}")
    checkrmsd=$( awk -v rthr="$rmsd_thr" '{if ($5<=rthr) print 1;else print 0}' "${cntrd_dir}/${avgtxt}") 

    if [ "$checkrmsd" -eq "1" ]; then
        #echo "Centroid consistent with itself"
        echo "${serial} ${repdir} ${lignmID} ${ev} ${sub} ${ref_lcl} ${b} ${e} ${dur} ${cov_lcl} ${avgrmsd} ${devrmsd}" >> "$acpt_event_report"
    else
        echo "${serial} ${repdir} ${lignmID} ${ev} ${sub} ${ref_lcl} ${b} ${e} ${dur} ${cov_lcl} ${avgrmsd} ${devrmsd}" >> "$rfsd_event_report"
    fi
    dict_refined="clusters_ev_${ev}_${method[$k]}_${thr[$k]}_ref_${refine_method}_${refine_thr}.txt"
    if [ ! -f "${ligdir}/${evdir}/${dict_refined}" ];then
        echo "#Event #Sub #Start[ns] #End[ns] #Duration #Centroid[ps] #Centroid_Coverage_pc" >${ligdir}/${evdir}/${dict_refined}
    fi
    echo "${ev} ${sub} ${b} ${e} ${dur} ${ref_lcl} ${cov_lcl}" >> ${ligdir}/${evdir}/${dict_refined}

    echo "${serial} ${repdir} ${lignmID} ${ev} ${sub} ${ref_lcl} ${b} ${e} ${dur} ${cov_lcl} ${avgrmsd} ${devrmsd}" >> "$all_event_report"
}

function binding_statistics {
    local ntrj=$1
    local nlig=$2
    #Retrieve p-bound data from every replica/ligand instance
    for ((j=1;j<=ntrj;j++));do
        while read lignr; do
            lignmID="${ligname}${lignr}"
            ligdir="rep${j}/${lignmID}"
            echo "rep${j} ${lignmID} $(tail -1  $ligdir/report_${lignmID}.txt)">>report_onelig.txt
            rm 	"$ligdir/report_${lignmID}.txt"
        done<"lignr_list.dat"
    done
    #Average
    awk -v nlig="$nlig" 'BEGIN{i=0;sumpb=0;j=0}
        { 	i++	
            pb[i]=($4/$3)
            sumpb+=pb[i]
            #printf("Considering i:%d j:%d\n",i,j)
            if(i==nlig)
                { 	avg[j]=sumpb/nlig
                    for(i=1;i<=nlig;i++)
                        sumq[j]+=((avg[j]-pb[i])**2)
                    if (nlig==1)
	        	devst[j]=sqrt(sumq[j]/(nlig))
                    else
		    devst[j]=sqrt(sumq[j]/(nlig-1))
		    #printf("Considering sumpb: %f avg:%f\n",sumpb,avg[j])
                    j++
                    i=0
                    sumpb=0								
                }  
        }
        END{printf("#RepID #Pbound-Avg #Stdev\n")
            for(k=0;k<j;k++)
                sum+=avg[k]
            avgtot=sum/j
            for(k=0;k<j;k++)
                sumqtot+=((avg[k]-avgtot)**2)
            if(j==1)
        	devsttot=sqrt(sumqtot)
        else	
        	devsttot=sqrt(sumqtot/(j-1))
            for(k=0;k<j;k++)
                printf("rep%d %s %s\n",k+1,avg[k],devst[k])
            printf("- - -\nOverall %f %f\n",avgtot,devsttot)}' report_onelig.txt > summary_onelig.txt
}
