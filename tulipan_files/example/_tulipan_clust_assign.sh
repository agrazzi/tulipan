#!/bin/bash
##
# Step 2 of TULIPAN analysis protocol
##
source "/usr/local/gromacs/bin/GMXRC"
source "/path/to/tulipan_files/functions/functions_clustassign.sh"
source "/path/to/tulipan_files/functions/functions_tulipan_histo.sh"


#Files for fitting procedure with python
blankfit="/path/to/tulipan_files/functions/temp_fit_prob.py"
blankfitdens="/path/to/tulipan_files/functions/temp_fit_dens.py"
temphisto="/path/to/tulipan_files/functions/temp_multihisto_v2.py"

#Files for VMD visualization
refbondpdb="/path/to/templates/prot_onelig.pdb"
blank="/path/to/templates/blank_v3.tcl"
blankpdb="/path/to/templates/blankpdb_v3.tcl"

Help(){
	echo -e "\nHints\n- Check path to necessary files\n- Check ligname\n- Check if method[] and thr[] have been set properly"
}

while getopts h option
do
case "${option}"
in
		h) Help;
           exit;;
esac
done

###Variable definition
    ##Workspace definition
        workdir="test_v4"
        N=6
    ##Biosystem
        ligname=COLC
    ##Refinement previously applied
        method=("gromos" )
        thr=("0.4")
    ##Parameters for classification of Binding Modes
        #Maximum average RMSD deviation for assignment to BM
        acc_thr="0.4"
    ##Residence time estimation
        #Minimum number of events for attempting histogram generation and res-time estimation
        min_ev_restime=10
        #Threshold for suspiciously high residence times
        tauthr="10000"
    ## Parameters for visualization of results
        #Minimum number of subevents attributed to the same supracluster for selection
        ext_thr=5
        #Minimum coverage
        mcvg=0.8 # poses covering at least 80% of specifc bound time will be extracted
        #Saving frequency (ps) for final trajectory concatenating all binding events belonging to the same pose
        catdt="5000"
        #Minimum number of poses that will be extracted
        mp=10
        #Maximum number of poses that will be highlighted for visualization with VMD
        extlimit="10"
        #Color ID (VMD-reference) associated to each pose; NOTE: this array must have same number of elements as $extlimit 
        vmdcolors=("0" "1" "2" "3" "4" "7" "9" "10" "11" "12")
###




dt=2





check_files "${refbondpdb}" "${blank}" "${blankpdb}" "${blankfit}" "${blankfitdens}" "${temphisto}"

if [ ! -d "$workdir" ];then
    echo "${workdir}/ not found"
 exit 1
fi

cd $workdir



#Begin time measurement
tbegin=$(date +"%T.%3N")

set -e
#Warning list for python fit
if [ -f "warning_list_fit.err" ];then
    rm "warning_list_fit.err"
fi

if [ -f "warning_list_fit_dns.err" ];then
    rm "warning_list_fit_dns.err"
fi

 

bck=0
bckdns=0
colorlimit="${#vmdcolors[@]}"
defcol="6"

for k in ${!method[@]}; do
    #method dir
    methoddir="attrib_${method[$k]}_${thr[$k]}"
    if [ -d ${methoddir} ];then
        rm -r ${methoddir}
    fi
    mkdir ${methoddir}
    #Top SupraCluster List
    topsc_list="${methoddir}/top_supracluster_list.txt"
    echo "#Domain #SupraclusterID #Events #AvgDuration #StDev #TopEventID #TopEventDuration">"${methoddir}/all_supracl_average_dur.txt"
    echo "#Domain #SupraclusterID #Events #AvgDuration #StDev #TopEventID #TopEventDuration">"${methoddir}/supracl_ll_topevent.txt"
    echo "#Event_domain #SupraclusterID #SP-Events">"${topsc_list}"
    framedir="frames_cntrd_${method[$k]}_${thr[$k]}"
    tot_domains=$(grep "Found" ${framedir}/gmx_cluster.log | awk '{print $2}')
    for ((i=1;i<=tot_domains;i++));do
        tgtdir="${methoddir}/event_domain_${i}"
        mkdir "${tgtdir}"
        echo -ne "${method[$k]}_${thr[$k]} Sanity Check on Event-Domain [$i]/[$tot_domains]...\t"
        ev_dom_list="events_cluster_domain_${i}.txt"
        if [ ! -f "${framedir}/${ev_dom_list}" ];then
            echo "events for domain ${i} not found"
            exit 1
        fi

        #Sanity check on existence of required tpr and trajectorie    
        centroid_exists
        status=$?
        if [[ $status -ne 0 ]]; then
            echo "centroid_exists failed with code $status"
            exit "$status"
        fi
        echo ""
        #Initialize unassigned for a certain domain
        cp "${framedir}/${ev_dom_list}" "${tgtdir}/unassigned.txt"
        ue=$(awk 'BEGIN{l=0}{if(substr($0,1,1)=="#")next;l++}END{print l}' "${tgtdir}/unassigned.txt")
        iter=0
        echo "#CentroidID #filePDB #fileTPR">"${tgtdir}/cntrd_list.txt"
        echo "#CD #Supracluster #Rep #LigID #Time[ps]">"${tgtdir}/retrieval.txt"
        echo "#cntrd_ID #Population" >"${tgtdir}/supracluster_size.log"
        echo "#cntrd_ID #Accepted #Refused">"${tgtdir}/tmsr_attribution.log"
        while [ "$ue" -gt "0" ]; do
            ((iter=iter+1))
            echo "Event attribution, iteration [$iter]; Unassigned event remaining: $ue"
            #Centroid for iter
            cntrd_tpr=$( awk -v mth="${method[$k]}" -v thr="${thr[$k]}" 'NR==2{printf("%s/%s/ev_%d/subevents_%s_%s/sub_%d/ev_%d_sub_%d_%s_%s_lc.tpr",$3,$4,$5,mth,thr,$6,$5,$6,mth,thr)}' ${tgtdir}/unassigned.txt )
            if [ ! -f "${cntrd_tpr}" ];then
                echo "Missing tpr: ${cntrd_tpr}"
                exit
            fi
            cntrd_pdb=$( awk -v mth="${method[$k]}" -v thr="${thr[$k]}" 'NR==2{printf("%s/%s/ev_%d/subevents_%s_%s/sub_%d/ev_%d_sub_%d_%s_%s_lc.pdb",$3,$4,$5,mth,thr,$6,$5,$6,mth,thr)}' ${tgtdir}/unassigned.txt )
            echo "$iter ${cntrd_pdb} ${cntrd_tpr}" >>"${tgtdir}/cntrd_list.txt"
            #Drop Info for retrieval
            awk -v iter="$iter" 'NR==2{printf("%d %d %s %s %s\n",$1,iter,$3,$4,$7)}' ${tgtdir}/unassigned.txt >>${tgtdir}/retrieval.txt
            
            submt=0
            while read input; do
                if [ "${input::1}" = "#" ];then
                    continue
                fi
                (event_rms "$input" "${cntrd_tpr}" "$dt")  &
                ((submt=submt+1))
                echo -ne "\rSubmitted [$submt]/[$ue]"
                sleep 0.1
                while [ "$(jobs -r -p | wc -l)" -ge "$N" ]; do
                    sleep 0.25
                done
            done<"${tgtdir}/unassigned.txt"

            echo ""
            echo -e "\rCheck if there are process still running before continuing..."

            while [ "$(jobs -r -p | wc -l)" -gt "0" ]; do
                    echo -ne "\rWaiting parallel tasks to complete..."
                    sleep 0.25
            done
            echo ""
            #Reorganize distances to table
            tabularize

            #Sanity check on table
            tbln=$(awk 'BEGIN{l=0}{if(substr($0,1,1)=="#")next;l++}END{print l}' "${tgtdir}/${table}")
            if [ "$ue" -ne "$tbln" ];then
                echo "Mismatch between unassinged event [$ue] and results stored in table:[$tbln]"
                exit
            fi

            #Parse the table and determine which events can be attributed to the candidate supracluster and which should be rejected
            read_table_distances
            
            #Regenerate unassigned and recompute ue
            ue=$(awk 'BEGIN{l=0}{if(substr($0,1,1)=="#")next;l++}END{print l}' "${tgtdir}/unassigned.txt" )
            cntrd_dir="cntrd_${iter}"
            mkdir ${tgtdir}/${cntrd_dir}
            mv ${tgtdir}/rmsd_*xvg "${tgtdir}/${cntrd_dir}/"
            if [ "$rfsdln" -eq "0" ]; then
                echo "No events remaining"
                break
            fi
        done
        #All sub-events within the domain have been attributed to a supracluster
        #Resorting
        rawfile="${tgtdir}/supracluster_size.log"
        resortfile="${tgtdir}/resort_supracluster_size.log"
        { head -n1 "$rawfile"; tail -n +2 "$rawfile" | sort -k2,2nr; } > ${resortfile}
        #
        #Extracting most important clusters from each domain
        awk -v thr="$ext_thr" -v domain="$i" '{if (substr($1,1,1)=="#") next; if ($2>=thr) printf("%d %d %d \n",domain,$1,$2)}' ${resortfile} >> ${topsc_list}
        
        #Computing time statistcs for supraclusters of a given domain
        avg_dur_supracluster

        #Clean-up         
        for ((j=1;j<=iter;j++)); do
            cntrdir="cntrd_${j}"
            if [ ! -f "${tgtdir}/supracluster_${j}.txt" ];then
                echo "Missing ${tgtdir}/supracluster_${j}.txt"
                exit 1
            else
                mv "${tgtdir}/supracluster_${j}.txt" "${tgtdir}/cntrd_${j}/supracluster_${j}.txt"
            fi
            
            if [ ! -f "${tgtdir}/rejected_${j}.txt" ];then
                                echo "Missing ${tgtdir}/rejected_${j}.txt"
                exit 1
            else
                mv "${tgtdir}/rejected_${j}.txt" "${tgtdir}/cntrd_${j}/rejected_${j}.txt"
            fi
            
            if [ ! -f "${tgtdir}/table_rmsd_dist_centroid_${j}.txt" ];then
                echo "Missing ${tgtdir}/table_rmsd_dist_centroid_${j}.txt"
                exit 1
            else
                mv "${tgtdir}/table_rmsd_dist_centroid_${j}.txt" "${tgtdir}/cntrd_${j}/table_rmsd_dist_centroid_${j}.txt"
            fi
        done
        if [ "$(< "${tgtdir}/unassigned.txt" wc -l)" -eq "1" ];then
            rm "${tgtdir}/unassigned.txt"
        else
            echo "Expecting ${tgtdir}/unassigned.txt to have no events, but still contains $(tail -n +2 "${tgtdir}/unassigned.txt" ) "
            exit 1
        fi
        awk 'NR>1{print $0}' "${tgtdir}/${dom_dur_tb}" >>"${methoddir}/all_supracl_average_dur.txt"
        topdur_ext_thr="400"
        awk -v thr="${topdur_ext_thr}" 'NR>1{if ($7>=thr) print $0}' "${tgtdir}/${dom_dur_tb}" >> "${methoddir}/supracl_ll_topevent.txt"
        #Extract Supraclusters with interesting duration
        #iteration over domains 
    done
    #iteration over methods
done

    echo "Pbound statistics"
    supracluster_coverage


echo "Residence Time Estimation"

if [ -d "error_logs/" ];then
    rm -r error_logs
    mkdir error_logs/
else
    mkdir error_logs/
fi


for k in ${!method[@]}; do
	methoddir="attrib_${method[$k]}_${thr[$k]}"
	if [ ! -f "${methoddir}/table_domains_cvg_resorted.txt" ];then
		echo "Missing table resorted"
		exit 1
	fi
    echo "#Domain #Supracluster #Events #AvgDur #StdDev #EventsSK #AvgDur #StdDev">"${methoddir}/avgreport.txt"
    echo "#Domain #Supracluster #Events #AvgDur #StdDev #EventsSK #AvgDur #StdDev #SturgesBins #SturgesBinW #DoaneBins #DoaneBinW #FreedmanBins #FreedmanBinW"> "${methoddir}/historeportdns.txt"
	echo "#Domain #Supracluster #Events #AvgDur #StdDev #EventsSK #AvgDur #StdDev #SturgesBins #SturgesBinW #DoaneBins #DoaneBinW #FreedmanBins #FreedmanBinW"> "${methoddir}/historeport.txt"
    echo "$(head -n 1 "${methoddir}/historeport.txt" ) #A_ST #E_A_ST #Tau_ST #E_Tau_ST #RMSE_ST #Estimator_ST #A_D #E_A_D #Tau_D #E_Tau_D #RMSE_D #Estimator_D #A_FD #E_A_FD #Tau_FD #E_Tau_FD #RMSE_FD #Estimator_FD"> "${methoddir}/fitreport.txt"
    echo "$(head -n 1 "${methoddir}/historeport.txt" ) #A_ST #E_A_ST #Tau_ST #E_Tau_ST #RMSE_ST #Estimator_ST #A_D #E_A_D #Tau_D #E_Tau_D #RMSE_D #Estimator_D #A_FD #E_A_FD #Tau_FD #E_Tau_FD #RMSE_FD #Estimator_FD"> "${methoddir}/fitreport_dns.txt"
    echo "#Domain #Supracluster #Events #SturgesNZBins #DoaneNZBins #FreedmanNZBins" >"${methoddir}/nzreport.txt"
    echo "#Domain #Supracluster #Events #SturgesNZBins #DoaneNZBins #FreedmanNZBins" >"${methoddir}/nzreport_dns.txt"
    ncluster=$(awk 'BEGIN{c=0}{if (substr($0,1,1)=="#")next; c++}END{print c}' "${methoddir}/table_domains_cvg_resorted.txt" )
	submt=0
	echo "Working on method [${method[$k]}] - ctf [${thr[$k]}] "
	while read -r evdmnID scid nev oth; do
		if [ "${evdmnID::1}" = "#" ];then
			continue
		else
			((submt=submt+1))
			echo -ne "\rFitting Residence Time for Binding Mode [$submt]/[$ncluster]"
			
            if [ "$nev" -lt "$min_ev_restime" ];then
                #Too few binding events for histogram generation
                continue
			else
                			
                supraclfile="${methoddir}/event_domain_${evdmnID}/cntrd_${scid}/supracluster_${scid}.txt"
				outpath="${methoddir}/event_domain_${evdmnID}/cntrd_${scid}"

				if [ ! -f "$supraclfile" ];then
					echo "Missing $supraclfile"
					exit 1
				fi
                
                ( histo_fit $submt $evdmnID $scid $supraclfile $outpath "$temphisto" "$tauthr" )&

                sleep 0.5
                while [ "$(jobs -r -p | wc -l)" -ge "$N" ]; do
                    sleep 0.5
                done       
                              
			fi			
		fi
	done<"${methoddir}/table_domains_cvg_resorted.txt"
    echo -e "\rCheck if there are processes still running before changing method..."
    while [ "$(jobs -r -p | wc -l)" -gt "0" ]; do
            echo -ne "\rWaiting..."
            sleep 1
    done
    submt=0
    echo -e "\nSummarising Results"
    while read evdmnID scid nev oth;do
        if [ "${evdmnID::1}" = "#" ];then
                continue
        else
            ((submt=submt+1))
            echo -ne "\r\tBinding Mode [$submt]/[$ncluster]"
            if [ "$nev" -lt "10" ];then
                continue
            else                              
                outpath="${methoddir}/event_domain_${evdmnID}/cntrd_${scid}"
                if [ ! -f "${outpath}/avg.dat" ];then
                    echo "Missing ${outpath}/avg.dat "
                    exit 1
                fi
                if [ ! -f "${outpath}/report_histo.txt" ];then
                    echo "Missing ${outpath}/report_histo.txt "
                    exit 1
                fi
                if [ ! -f "${outpath}/fitresults.txt" ];then
                    echo "Missing ${outpath}/fitresults.txt"
                    exit 1
                fi
                if [ ! -f "${outpath}/report_dns_histo.txt" ];then
                    echo "Missing ${outpath}/report_dns_histo.txt "
                    exit 1
                fi
                if [ ! -f "${outpath}/fitresults_dns.txt" ];then
                    echo "Missing ${outpath}/fitresults_dns.txt"
                    exit 1
                fi
                tail -n 1 "${outpath}/avg.dat" >> "${methoddir}/avgreport.txt"
                echo "$(tail -n 1 "${outpath}/avg.dat") $(tail -n 1 "${outpath}/report_histo.txt")" >> "${methoddir}/historeport.txt"
                echo "$(tail -n 1 "${outpath}/avg.dat") $(tail -n 1 "${outpath}/report_dns_histo.txt")" >> "${methoddir}/historeportdns.txt"
                echo "$(tail -n 1 "${methoddir}/historeport.txt") $(tail -n 1 "${outpath}/fitresults.txt")" >> "${methoddir}/fitreport.txt"
                echo "$(tail -n 1 "${methoddir}/historeportdns.txt") $(tail -n 1 "${outpath}/fitresults_dns.txt")" >> "${methoddir}/fitreport_dns.txt"
                echo "$evdmnID $scid $nev $(tail -n 1 "${outpath}/nzbins.txt")" >> "${methoddir}/nzreport.txt"
                echo "$evdmnID $scid $nev $(tail -n 1 "${outpath}/nzbins_dns.txt")" >> "${methoddir}/nzreport_dns.txt"
                if [ -f "error_${submt}.log" ];then #if error_$.log is not present, the python fit has not been performed, no need to look for Traceback
                    if grep -q "Traceback" "error_${submt}.log"; then
                        ((bck=bck+1))
                        mv  error_${submt}.log error_logs/error_${submt}_bck_${bck}.log
                        mv fit_${submt}.py     error_logs/fit_${submt}_${bck}.py
                        mv multihisto_${submt}.py error_logs/multihisto_${submt}_${bck}.py
                        echo "Warning for error_${submt}.log ;BackedUp as error_${submt}_bck_${bck}.log" >>"warning_list_fit.err"
                    else
                        rm "error_${submt}.log" fit_${submt}.py multihisto_${submt}.py
                    fi
                else
                    if [ -f "multihisto_${submt}.py" ];then
                        rm multihisto_${submt}.py
                    fi
                fi
                if [ -f "error_${submt}_dns.log" ];then #if error_$.log is not present, the python fit has not been performed, no need to look for Traceback
                    if grep -q "Traceback" "error_${submt}_dns.log"; then
                        ((bckdns=bckdns+1))
                        mv  error_${submt}_dns.log error_logs/error_${submt}_dns_bck_${bckdns}.log
                        mv fit_${submt}_dens.py     error_logs/fit_${submt}_dens_${bckdns}.py
                        echo "Warning for error_${submt}.log ;BackedUp as error_${submt}_bck_${bck}.log" >>"warning_list_fit_dns.err"
                    else
                        rm "error_${submt}_dns.log" fit_${submt}_dens.py 
                    fi
                fi
            fi
        fi
    done<"${methoddir}/table_domains_cvg_resorted.txt"
    echo ""
    #Sanity check on size of estimator error

    echo "#Domain #Supracluster #Events #TotalDuration #Coverage #CumulativeCoverage">"${methoddir}/table_domains_cvg_rstd_cum.txt"
    awk 'BEGIN{sum=0}NR>1{ sum+=$7; print $1 OFS $2 OFS $3 OFS $4 OFS $7 OFS sum }' "${methoddir}/table_domains_cvg_resorted.txt">>"${methoddir}/table_domains_cvg_rstd_cum.txt"
    #Choice of top rule for histogram
    top_rule_histo
    
    echo "#Domain #Supracluster #Events #Duration #Cvg #Cvg-Cum #AvgSk #DevSk #Rule #Tau #ErrTau #ExtractionID"> "${methoddir}/extraction_selection.txt"

    awk -v mp="$mp" -v mcvg="$mcvg" 'BEGIN{c=0}{if(substr($0,1,1)=="#")next;c++; if (c<=mp || $7 <= mcvg) printf("%s%s%d\n",$0,OFS,c)}' "${methoddir}/table_domains_cvg_rstd_cum_RT.txt" >> "${methoddir}/extraction_selection.txt"
	awk -v el="$extlimit" 'BEGIN{el++; } NR==1{print $0 " #Label"}NR>1 && NR<=el {printf("%s D%dS%d\n",$0,$1,$2)}' "${methoddir}/extraction_selection.txt" > "${methoddir}/extraction_selection_short.txt"
    
done

echo "Extracting poses"

for k in ${!method[@]}; do
    echo -e "\nWorking on method [${method[$k]}] - ctf [${thr[$k]}]"
    methoddir="attrib_${method[$k]}_${thr[$k]}"
    if ls "${methoddir}"/ext_bm_* &>/dev/null;then
        rm -r ${methoddir}/ext_bm_*
    fi
    count="0"
    nextractions=$(awk 'BEGIN{c=0}{if(substr($0,1,1)=="#")next;c++;}END{print c}' "${methoddir}/extraction_selection.txt" )
    while read line;do
        if [ "${line::1}" = "#" ];then
            continue
        else
            evdmnID=$(echo "$line" | awk '{print $1}')
            scid=$(echo "$line" | awk '{print $2}')
            outid=$(echo "$line" | awk '{print $12}')
            outdir="${methoddir}/ext_bm_${outid}"
            echo -ne "\rWorking on BM [$outid]/[$nextractions]"
            if [ ! -d "$outdir" ]; then
                mkdir $outdir
            else
                echo "Trying to reuse same outid; check input requests"
                exit
            fi
            
            evdomain="event_domain_${evdmnID}"
            sclist="supracluster_${scid}.txt"
            cntrd_dir="cntrd_${scid}"
            path="${methoddir}/${evdomain}/${cntrd_dir}"
            
            if [ ! -f "${path}/${sclist}" ];then
                echo "Missing ${path}/${sclist}"
                exit
            fi
            while read -r rep lignmID ev sub oth; do
                if [ "${rep::1}" = "#" ];then
                    continue
                else
                    xtc_dir="${rep}/${lignmID}/ev_${ev}/subevents_${method[$k]}_${thr[$k]}/sub_${sub}"
                    xtc_subev="${lignmID}_ev_${ev}_sub_${sub}.xtc"
                    if [ -f "${xtc_dir}/${xtc_subev}" ];then
                        filepath=$(realpath "${xtc_dir}/${xtc_subev}")
                        echo -n  "$filepath " >>listcat.txt
                        echo "c" >>befehle.txt
                    else
                        echo "${xtc_dir}/${xtc_subev} not found"
                        exit
                    fi
                fi
            done<"${path}/${sclist}"
        fi
        
        cat "befehle.txt" | gmx trjcat -cat -settime -f $(cat "listcat.txt") -dt "$catdt" -o ${outdir}/cat.xtc > /dev/null 2> log_cat_extraction.log
        rm listcat.txt
        rm befehle.txt
        cntrdlist="${methoddir}/${evdomain}/cntrd_list.txt"
        pdbfile=$( awk -v scid="$scid" '{if (substr($0,1,1)=="#")next; if ($1==scid)print $2}' ${cntrdlist})
        if [ "$pdbfile" = "" ];then
            echo "PDB not read correctly"
            exit
        fi
        if [ ! -f "${pdbfile}" ]; then
            echo "PDB not found"
            exit
        fi
        
        cp "$pdbfile" ${outdir}/
        pdbfile=$(realpath ${outdir}/*pdb)
        xtcfile=$(realpath ${outdir}/*xtc)      
        sed -i "s/ENDMDL//g" "${pdbfile}"
	    #Retrieve CONECT record from reference pdb
        grep "CONECT" ${refbondpdb} >> "${pdbfile}"

        if [ "$count" -lt "$colorlimit" ];then
        	pair[$count]="{ \"${pdbfile}\" \"${xtcfile}\" \"${vmdcolors[$count]}\" }"
        else
        	pair[$count]="{ \"${pdbfile}\" \"${xtcfile}\" \"${defcol}\" }"
        fi
        ((count=count+1))
    done<"${methoddir}/extraction_selection.txt"

    echo -e  "\nPreparing VMD .tcl file"
 
    cp ${blank} ${methoddir}/load_extraction.tcl
    sed -i "s/LIGNAME/$ligname/g" ${methoddir}/load_extraction.tcl
    for ((i=0;i<count;i++)); do
        if [ "$i" -lt $((count-1)) ];then
        
            sed -i "s|TEMP|${pair[$i]}\nTEMP|g" ${methoddir}/load_extraction.tcl 
        else
            sed -i "s|TEMP|${pair[$i]}|g" ${methoddir}/load_extraction.tcl
        fi
    done
    unset pair[@]

done

################
echo "Preparing selection of top-${extlimit} poses"

for k in ${!method[@]}; do
    echo -e "\nWorking on method [${method[$k]}] - ctf [${thr[$k]}]"
    methoddir="attrib_${method[$k]}_${thr[$k]}"
    nextractions=$(awk 'BEGIN{c=0}{if(substr($0,1,1)=="#")next;c++;}END{print c}' "${methoddir}/extraction_selection_short.txt" )
	count=0
    while read line;do
        if [ "${line::1}" = "#" ];then
            continue
        else
            outid=$(echo "$line" | awk '{print $12}')
            echo -ne "\rWorking on BM [$outid]/[$nextractions]"     
            outdir="${methoddir}/ext_bm_${outid}"
        fi
        pdbfile=$(realpath ${outdir}/*pdb)
        xtcfile=$(realpath ${outdir}/*xtc)
        echo -ne ""
        if [ "$count" -lt "$colorlimit" ];then
        	pair[$count]="{ \"${pdbfile}\" \"${xtcfile}\" \"${vmdcolors[$count]}\" }"
            pairpdb[$count]="{ \"${pdbfile}\" \"${vmdcolors[$count]}\" }"
        else
        	pair[$count]="{ \"${pdbfile}\" \"${xtcfile}\" \"${defcol}\" }"
            pairpdb[$count]="{ \"${pdbfile}\" \"${defcol}\" }"
        fi
        ((count=count+1))
	done<"${methoddir}/extraction_selection_short.txt"

    echo -e  "\nPreparing VMD .tcl file" 
    cp ${blank} ${methoddir}/load_extraction_selection.tcl
    cp ${blankpdb} ${methoddir}/load_extraction_selection_pdb.tcl
    sed -i "s/LIGNAME/$ligname/g"  ${methoddir}/load_extraction_selection.tcl
    sed -i "s/LIGNAME/$ligname/g" ${methoddir}/load_extraction_selection_pdb.tcl
    for ((i=0;i<count;i++)); do
        if [ "$i" -lt $((count-1)) ];then
            sed -i "s|TEMP|${pair[$i]}\nTEMP|g" "${methoddir}/load_extraction_selection.tcl" 
            sed -i "s|TEMP|${pairpdb[$i]}\nTEMP|g" "${methoddir}/load_extraction_selection_pdb.tcl"
        else
            sed -i "s|TEMP|${pair[$i]}|g" "${methoddir}/load_extraction_selection.tcl"
            sed -i "s|TEMP|${pairpdb[$i]}|g" "${methoddir}/load_extraction_selection_pdb.tcl"
        fi
    done
    unset pair[@]
    unset pairpdb[@]

done

echo -e "\nNORMAL TERMINATION"

cd .. 

#END time measurement
tend=$(date +"%T.%3N")
echo -e "\nScript began at: $tbegin\nScript ended at: $tend"
