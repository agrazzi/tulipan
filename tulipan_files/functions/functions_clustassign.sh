#!/bin/bash

function check_files {
if [ "$#" -eq 0 ];then
    echo "No arguments passed to check function"
    exit 1
else
    for file in "$@"; do
        if [ ! -f "$file" ];then
            echo -e "\nMissing file: $file"
            exit 1
        fi
        
    done
fi
}

function centroid_exists {
    while read -r cd serial rep lignmID ev sub clref; do
        prot_lig_ndx="${rep}/${lignmID}/prot_${lignmID}.ndx"
        if [ "${cd::1}" = "#" ];then
            continue
        fi
        cntrd_tpr="${rep}/${lignmID}/ev_${ev}/subevents_${method[$k]}_${thr[$k]}/sub_${sub}/ev_${ev}_sub_${sub}_${method[$k]}_${thr[$k]}_lc.tpr"
        
        cntrd_pdb="${rep}/${lignmID}/ev_${ev}/subevents_${method[$k]}_${thr[$k]}/sub_${sub}/ev_${ev}_sub_${sub}_${method[$k]}_${thr[$k]}_lc.pdb"
        subevent_xtc="${rep}/${lignmID}/ev_${ev}/subevents_${method[$k]}_${thr[$k]}/sub_${sub}/${lignmID}_ev_${ev}_sub_${sub}.xtc"
        if [ ! -f "${cntrd_tpr}" ];then
            echo "Missing tpr: ${cntrd_tpr}"
            return 1
        fi
        if [ ! -f "${cntrd_pdb}" ];then
            echo "Missing pdb: ${cntrd_pdb}"
            return 2
        fi
        if [ ! -f "${subevent_xtc}" ];then
            echo "Missing xtc: ${subevent_xtc}"
            return 3
        fi
    done<"${framedir}/${ev_dom_list}"
    
}

#event_analysis $line ${cntrd_tpr}
function event_rms {     
                local serial=$(echo "$1" | awk '{print $2}')
                local rep=$(echo "$1" | awk '{print $3}')
                local lignmID=$(echo "$1" | awk '{print $4}')
                local ev=$(echo "$1" | awk '{print $5}')
                local sub=$(echo "$1" | awk '{print $6}')
                local prot_lig_ndx="${rep}/${lignmID}/prot_${lignmID}.ndx"
                local outxvg="rmsd_cntrd_${iter}_vs_${rep}_${lignmID}_ev${ev}_sub${sub}.xvg"
                local subevent_xtc="${rep}/${lignmID}/ev_${ev}/subevents_${method[$k]}_${thr[$k]}/sub_${sub}/${lignmID}_ev_${ev}_sub_${sub}.xtc"
                echo "${lignmID}" | gmx rms -f ${subevent_xtc} -s $2 -o ${tgtdir}/${outxvg} -n  ${prot_lig_ndx} -fit none > /dev/null 2> log_rms.log
                local avgtxt="rmsd_cntrd_${iter}_vs_${rep}_${lignmID}_ev${ev}_sub${sub}.txt"
                awk -v rep="$rep" -v lig="$lignmID" -v ev="$ev" -v subev="$sub" -v sr="$serial" '
                ARGIND==1{if (substr($0,1,1)=="#" || substr($0,1,1)=="@") next; c++;sum+=$2} 
                ARGIND==2 && FNR==1 {avg=sum/c}
                ARGIND==2{if (substr($0,1,1)=="#" || substr($0,1,1)=="@") next; sumq+=((avg-$2)^2)}
                END{printf("%s %s %d %d %f %f %d %d\n",rep,lig,ev,subev,avg,sqrt(sumq/(c-1)),c,sr)}' ${tgtdir}/${outxvg} ${tgtdir}/${outxvg} > "${tgtdir}/${avgtxt}"
                if [ ! -f "${tgtdir}/${avgtxt}" ]; then
                    echo "inside Missing ${tgtdir}/${avgtxt}"
                    exit
                fi
}

function tabularize {
    table="table_rmsd_dist_centroid_${iter}.txt"
    echo "#Replica #Ligand #Event #Subevent #Avg #Stdev #Frames #Serial">${tgtdir}/${table}
    while read -r cd sps rep lignmID ev sub oth; do
        if [ "${cd::1}" = "#" ];then
            continue
        fi
        avgtxt="rmsd_cntrd_${iter}_vs_${rep}_${lignmID}_ev${ev}_sub${sub}.txt"
        if [ -f "${tgtdir}/${avgtxt}" ];then
            cat "${tgtdir}/${avgtxt}" >> ${tgtdir}/${table}
            rm  "${tgtdir}/${avgtxt}"
        else
            echo "Missing ${tgtdir}/${avgtxt}"
            exit
            if [ ! -f "${tgtdir}/lost.txt" ];then
                echo "Centroid#ID Requested_unsuccesfully">${tgtdir}/lost.txt
            fi
            echo "${iter} ${tgtdir}/${avgtxt}" >> ${tgtdir}/lost.txt
            exit
        fi 
    done<"${tgtdir}/unassigned.txt"
    }
    
function read_table_distances {
#Parse the table and determine which events can be attributed to the candidate supracluster and which should be rejected
awk -v thr="$acc_thr" -v iter="$iter" -v dir="${tgtdir}" '
        BEGIN{  out_acpt=sprintf("%s/supracluster_%d.txt",dir,iter);
            
            out_rfsd=sprintf("%s/rejected_%d.txt",dir,iter)
        }
        NR==1 {print $0 >>out_acpt;print $0 >>out_rfsd }
        {   if(substr($0,1,1)=="#")next;
            if($5<=thr)
                print $0 >>out_acpt;
            else
                print $0 >>out_rfsd;}' ${tgtdir}/${table}
    acpt_list="supracluster_${iter}.txt"
    rfsd_list="rejected_${iter}.txt"
    acptln=$(awk 'BEGIN{l=0}{if(substr($0,1,1)=="#")next;l++}END{print l}' "${tgtdir}/${acpt_list}")
    rfsdln=$(awk 'BEGIN{l=0}{if(substr($0,1,1)=="#")next;l++}END{print l}' "${tgtdir}/${rfsd_list}")
    checkawk=$(echo "" | awk -v acpt="$acptln" -v rfsd="$rfsdln" -v ue="$ue" '{sum=acpt+rfsd; if (sum!=ue)print 1;else print 0}')
    
    if [ "$checkawk" -eq "1" ];then
        echo "Mismatch between considered UnassignedEvents [$ue] and classified (accepted[$acptln]; refused[$rfsdln])"
        exit
    fi

    echo "${iter} ${acptln}">>"${tgtdir}/supracluster_size.log"
    echo "${iter} ${acptln} ${rfsdln}">>"${tgtdir}/tmsr_attribution.log"


    awk '   BEGIN{dump=0;c=0;}
            ARGIND==2 && FNR==1 {print $0}
            ARGIND==1{if (substr($1,1,1)=="#")next; l++; sr[l]=$8;}
            ARGIND==2 && FNR>1{      if (substr($1,1,1)=="#")next; 
                            if (c==l) #already recovered all events attributed, now dump everything
                                            dump=1
                            if (dump == 1)
                                    print $0;
                            else
                            {       for(j=1;j<=l;j++)
                                            {   if (sr[j]==$2)#found attributed
                                                    {c++;next}}
                                    #not among attributed
                                    print $0                                        
                                }
        }' "${tgtdir}/${acpt_list}" ${tgtdir}/unassigned.txt > ${tgtdir}/temp.txt
    mv ${tgtdir}/temp.txt ${tgtdir}/unassigned.txt
}

function avg_dur_supracluster {
    dom_dur_tb="table_supracl_dur_domain_$i.txt"
    echo "#Domain #SupraclusterID #Events #AvgDuration #StDev #TopEventID #TopEventDuration">"${tgtdir}/${dom_dur_tb}"
    for ((j=1;j<=iter;j++)); do
        acpt_list="supracluster_${j}.txt"
        gawk -v dmn="$i" -v scid="$j" -v dir="$tgtdir" 'BEGIN{maxdur=0;topid="";sum=0;sumq=0;outfile=sprintf("%s/table_supracl_dur_domain_%d.txt",dir,dmn)}
            {if(substr($0,1,1)=="#") next
            l++
            dur[l]=($7-1) #-1 to account for difference between frames and ns
            sum+=($7-1)
            if (dur[l]>maxdur)
                    {maxdur=dur[l]
                    topid=$8}
            }
        END{    avg=sum/l; 
            if (l>1)
            {       for(j=1;j<=l;j++)
                            sumq+=((avg-dur[j])^2)
                    dev=sqrt(sumq/(l-1))
                    printf("%s %s %d %.2f %.2f %d %d\n",dmn,scid,l,avg,dev,topid,maxdur)>>outfile
            }
            else
            printf("%s %s %d %f 0 %d %d\n",dmn,scid,l,avg,topid,maxdur)>>outfile
            }' "${tgtdir}/${acpt_list}"
        
    done
            }
            
function supracluster_coverage {
    #Determine n-trajectory from number of rep/ folders
    ntrj=$( find . -maxdepth 1 -type d -name "rep*" | wc -l)
    echo "Assuming $ntrj rep/ folders"
    if [ -f "listdir.txt" ];then
        rm listdir.txt
    fi
    for ((i=1;i<=ntrj;i++));do
        if [ ! -d "rep$i" ];then
            echo "Unable to locate rep$i"
            exit
        fi
        echo "rep$i">>listdir.txt
    done

    for k in ${!method[@]}; do
        evlist="acpt_sp_events_${method[$k]}_${thr[$k]}.txt"    
        if [ ! -f "$evlist" ];then
            echo "Missing $evlist"
            exit 1
        fi
        
        gawk '
            ARGIND==1 {rep[NR]=$1;totrep++}
            ARGIND==2 {    if (substr($0,1,1)=="#")next;
                                for (i=1;i<=totrep;i++)
                                    {   if ($1==rep[i])
                                            {   repsumbound[i]+=$4
                                                repsumparsed[i]+=$3
                                                break
                                            }   
                                    }
                        }
            ARGIND==3{if (substr($0,1,1)=="#")next;
                    for (i=1;i<=totrep;i++)
                        {   if ($2==rep[i])
                                {   repsumsp[i]+=$9
                                    break
                                }   
                        }
                        }
            END{
            printf("#Replica #TotFramesParsed #TotFramesAspBound #TotFramesSPBound #SPBoundTotCoverage #SPBoundBoundCoverage\n")
            for (i=1;i<=totrep;i++)
                {   
                    totcvg[i]=(repsumsp[i]/repsumparsed[i])
                    boundcvg[i]=(repsumsp[i]/repsumbound[i])
                    printf("%s %d %d %d %.4f %.4f\n",rep[i],repsumparsed[i],repsumbound[i],repsumsp[i],totcvg[i],boundcvg[i])          
                }
            for (i=1;i<=totrep;i++)
                {   sumtotcvg+=totcvg[i]
                    sumboundcvg+=boundcvg[i]

                }
            avgtotcvg=sumtotcvg/totrep
            avgboundcvg=sumboundcvg/totrep
            if (totrep>1)
                {   for (i=1;i<=totrep;i++)
                        {   sumqtotcvg+=((totcvg[i]-avgtotcvg)^2)
                            sumqboundcvg+=((boundcvg[i]-avgboundcvg)^2)
                        }
                    devtotcvg=sqrt(sumqtotcvg/(totrep-1))
                    devboundcvg=sqrt(sumqboundcvg/(totrep-1))
                    printf("#Overall #Avg-SPBTCVg #DevSt  #Avg-SPBBCvg #DevSt\n--> %.4f %.4f %.4f %.4f\n",avgtotcvg,devtotcvg,avgboundcvg,devboundcvg)
                }
            else
                printf("#Overall #Avg-SPBTCVg #DevSt  #Avg-SPBBCvg #DevSt\n--> %.4f - %.4f -\n",avgtotcvg,avgboundcvg) } ' listdir.txt  report_onelig.txt "$evlist" >"spbound_coverage_${method[$k]}_${thr[$k]}.txt"


        methoddir="attrib_${method[$k]}_${thr[$k]}"

        table_dmn_cvg="table_domains_cvg.txt"
        table_dmn_cvg_srtd="table_domains_cvg_resorted.txt"
        echo "#Domain #Supracluster #Events #TotalDuration #Avg #DevSt #Coverage">"${methoddir}/${table_dmn_cvg}"
        echo "#Domain #Supracluster #Events #TotalDuration #Avg #DevSt #Coverage">"${methoddir}/${table_dmn_cvg_srtd}"
        set -e
        framedir="frames_cntrd_${method[$k]}_${thr[$k]}"
        tot_domains=$(grep "Found" ${framedir}/gmx_cluster.log | awk '{print $2}')

        if [ -f "${methoddir}/temp_dmn_cvg.txt" ];then
            rm "${methoddir}/temp_dmn_cvg.txt"
        fi
        for ((i=1;i<=tot_domains;i++));do
            tgtdir="${methoddir}/event_domain_${i}"
            if [ ! -d "${tgtdir}" ];then
                    exit
            fi
            echo -ne "\r${method[$k]}_${thr[$k]} Working on domain [$i]/[$tot_domains]"

            nsupra=$(awk 'BEGIN{l=0}{if(substr($0,1,1)=="#")next;l++}END{print l}' "${tgtdir}/cntrd_list.txt")
            if [ "$nsupra" -eq "0" ];then
                exit
            fi
            if [ -f "${tgtdir}/temp_supracluster_duration.txt" ];then
                rm "${tgtdir}/temp_supracluster_duration.txt"
            fi 

            for ((j=1;j<=nsupra;j++));do
                cntrd_dir="cntrd_${j}"
                supralist="supracluster_${j}.txt"
                if [ ! -f "${tgtdir}/${cntrd_dir}/$supralist" ];then
                    exit
                fi
                awk -v i="$j" 'BEGIN{sum=0;ev=0}{if(substr($0,1,1)=="#")next;sum+=($7-1);ev++}END{print i OFS ev OFS sum}' "${tgtdir}/${cntrd_dir}/$supralist" >>"${tgtdir}/temp_supracluster_duration.txt"
            done
                awk   -v dmn=$i 'BEGIN{l=0}
                    ARGIND==1{if(substr($0,1,1)=="#")next;ev[$1]=$2;sum[$1]=$3;l++}
                    ARGIND==2{if(substr($0,1,1)=="#")next;avg[$2]=$4;dev[$2]=$5}
                    END {   printf("#Domain #Supracluster #Events #TotalDuration #Avg #DevSt\n")
                            for (i=1;i<=l;i++)
                            {   printf("%d %d %d %d %s %s\n",dmn,i,ev[i],sum[i],avg[i],dev[i])}
                        }' "${tgtdir}/temp_supracluster_duration.txt" "${tgtdir}/table_supracl_dur_domain_${i}.txt" >"${tgtdir}/table_domain_duration.txt"
                        tail -n +2 "${tgtdir}/table_domain_duration.txt">>"${methoddir}/temp_dmn_cvg.txt"
                #compare how the supracluster cover the specific events refined of each method
                
        done
        awk 'BEGIN{sum=0}
            ARGIND==1   {if(substr($0,1,1)=="#")next; sum+=$9}
            ARGIND==2   {if(substr($0,1,1)=="#")next; printf("%s %.4f\n",$0,$4/sum)}' ${evlist}  "${methoddir}/temp_dmn_cvg.txt" >> "${methoddir}/${table_dmn_cvg}"
        #Resorting
            rawfile="${methoddir}/${table_dmn_cvg}"
            resortfile="${methoddir}/${table_dmn_cvg_srtd}"
            { head -n1 "$rawfile"; tail -n +2 "$rawfile" | sort -k7,7nr; } > ${resortfile}
    
        echo ""
    done
}

function top_rule_histo {
awk 'BEGIN{c=0;dict[1]=17; dict[2]=23; dict[3]=29; rule[1]="Sturges" ; rule[2]="Doane"; rule[3]="FD";print "#Domain #Supracluster #Events #Duration #Cvg #Cvg-Cum #AvgSk #DevSk #Rule #Tau #ErrTau"}
        ARGIND==1 && FNR >1 {  
            top[FNR]=0; max=0
            for(i=4;i<=NF;i++)
                {if($i>=max)
                    {   max=$i
                        if (top[FNR]==0) #initialize backup
                            {   bck[FNR]=(i-3)
                                bckrule[FNR]=rule[(i-3)] }
                        else
                            {   bck[FNR]=top[FNR]
                                bckrule[FNR]=toprule[FNR]   }
                        top[FNR]=(i-3)
                        toprule[FNR]=rule[(i-3)]
                    } 
                }
        }
        ARGIND==2 && FNR>1{
            j=dict[top[FNR]]; 
            k=j+1;
            p=j+3;
            if ($j==-1)
                {tau[FNR]="-"
                ertau[FNR]="-";
                }
            else 
                {   if ($p<1 && $p>-1) #Error on estimator acceptable
                        {tau[FNR]=$(j);
                        ertau[FNR]=$k;}
                    else    #Check if backup rule yields better results
                        {   #Redefine target columns
                            j=dict[bck[FNR]]; 
                            k=j+1;
                            p=j+3;
                            if ($j==-1)
                                {   #Fit has not been performed with the backup rule
                                    tau[FNR]="-"
                                    ertau[FNR]="-";                            
                                }
                            else #fit available for backup, check magnitude on error of estimator
                            {   if ($p<1 && $p>-1)
                                    {   #backup is better; update info
                                        tau[FNR]=$(j);
                                        ertau[FNR]=$k;
                                        toprule[FNR]=bckrule[FNR]
                                    }
                                else
                                    {   #Even backup is wrong
                                        tau[FNR]="-"
                                        ertau[FNR]="-";
                                    }
                            }
                        }
                }
            avgsk[FNR]=$7
            devsk[FNR]=$8
            }        
        ARGIND==3 && FNR>1 {
            if ($3>=10)
                {l++;tmp=l+1
                printf("%s %.1f %.1f %s %s %s\n",$0,avgsk[tmp],devsk[tmp],toprule[tmp],tau[tmp],ertau[tmp]) }
            else
                printf("%s - - - - -\n",$0)
        }' "${methoddir}/nzreport.txt" "${methoddir}/fitreport.txt" "${methoddir}/table_domains_cvg_rstd_cum.txt" > "${methoddir}/table_domains_cvg_rstd_cum_RT.txt"
}

function histo_fit {
    local submt=$1
    local evdmnID=$2
    local scid=$3
    local supraclfile=$4
    local outpath=$5
    local blank_mlthst=$6

    awk -v evdmn="$evdmnID" -v scid="$scid" 'BEGIN{sum=0;sk_sum=0;c=0}{if (substr($1,1,1)=="#" || substr($1,1,1)=="@")next;
                data[NR]=$7
                sum+=$7
                c++}
    END{smpavg=sum/c; #simple avg
        asort(data)
        if (c>20)
            {   k=c*0.05
                if (k-int(k)>0.5)
                    skip=(int(k)+1)
                else
                    skip=int(k)
            }
        else
            { skip=1
            }
        ncsk=c-2*skip
        for (i=1+skip;i<=(c-skip);i++)
            { sk_sum+=data[i]
            }
        skavg=sk_sum/ncsk #skip avg
        if (c>1)
        {for (i=1;i<=c;i++)
            sumq+=((smpavg-data[i])^2)
        smpdev=sqrt(sumq/(c-1))}
        else
            smpdev=0
        sumq=0
        if (ncsk>1)
        {for (i=1+skip;i<=(c-skip);i++)
            sumq+=((skavg-data[i])^2)
        skdev=sqrt(sumq/(ncsk-1))}
        else
            skdev=0
        
        printf("%s %s %d %.1f %.1f %d %.1f %.1f\n",evdmn,scid,c,smpavg,smpdev,ncsk,skavg,skdev)
        }' ${supraclfile} >"${outpath}/avg.dat"
    local smpAvg=$(cat "${outpath}/avg.dat" | awk '{print $4}')
                    
    cp  ${blank_mlthst} multihisto_${submt}.py
    
    local sedoutpath="${methoddir}/event_domain_${evdmnID}/cntrd_${scid}"
    
    sed -i "s|OUTPATH|${sedoutpath}|g" multihisto_${submt}.py
    python3 multihisto_${submt}.py ${supraclfile} 7

    #Fit different histograms 
    local histlist=("histo_stu.txt" "histo_doane.txt" "histo_fd.txt")
    local histlistdns=("histo_dns_stu.txt" "histo_dns_doane.txt" "histo_dns_fd.txt")
    local fitplt=("fit_stu.png" "fit_doane.png" "fit_fd.png")
    local fitpltdns=("fit_dns_stu.png" "fit_dns_doane.png" "fit_dns_fd.png")
    local fitprm=("fit_stu.txt" "fit_doane.txt" "fit_fd.txt")
    local fitprmdns=("fit_dns_stu.txt" "fit_dns_doane.txt" "fit_dns_fd.txt")
    local fitdata=""
    local fitdatadns=""
    local k
    echo "#A_ST #E_A_ST #Tau_ST #E_Tau_ST #RMSE_ST #Estimator_ST #A_D #E_A_D #Tau_D #E_Tau_D #RMSE_D #Estimator_D #A_FD #E_A_FD #Tau_FD #E_Tau_FD #RMSE_FD #Estimator_FD" >"${outpath}/fitresults.txt"
    echo "#A_ST #E_A_ST #Tau_ST #E_Tau_ST #RMSE_ST #Estimator_ST #A_D #E_A_D #Tau_D #E_Tau_D #RMSE_D #Estimator_D #A_FD #E_A_FD #Tau_FD #E_Tau_FD #RMSE_FD #Estimator_FD" >"${outpath}/fitresults_dns.txt"
    for k in ${!histlist[@]};do
        if [ ! -f ${outpath}/${histlist[$k]} ];then
            echo "Missing ${outpath}/${histlist[$k]}"
            exit 1
        fi
        local nzhist=$(echo "${histlist[$k]}" | awk '{printf("%s_nz.txt",substr($1,1,length($1)-4))}')
        local nzhistdns=$(echo "${histlistdns[$k]}" | awk '{printf("%s_nz.txt",substr($1,1,length($1)-4))}')
        #Generate Non-Zero histogram
        awk '{if ($2) print $0}' ${outpath}/${histlist[$k]} > ${outpath}/${nzhist}
        awk '{if ($2) print $0}' ${outpath}/${histlistdns[$k]} > ${outpath}/${nzhistdns}

        #Check size of Non-Zero histogram
        local nl=$(< "${outpath}/${nzhist}" wc -l)
        local nldns=$(< "${outpath}/${nzhistdns}" wc -l)
        local nznl[$k]="$nl"
        if [ "$nl" -ge "4" ];then #Enough non-zero bins to attempt fit
            cp ${blankfit} fit_${submt}.py
            local sedoutpath="${methoddir}/event_domain_${evdmnID}/cntrd_${scid}"
            sed -i "s|OUTPATH|${sedoutpath}|g" fit_${submt}.py
            sed -i "s|HISTO|${nzhist}|g" fit_${submt}.py
            sed -i "s|PLOT|${fitplt[$k]}|g" fit_${submt}.py
            sed -i "s|FITPRM|${fitprm[$k]}|g" fit_${submt}.py

            python3 fit_${submt}.py > /dev/null 2> error_${submt}.log    
        elif [ "$nl" -ge "1" ]; then
            echo -e "Only $nl non-zero bins detected: fit not performed\nPrinting filler value for the parameters\nA = -1 - -1\nTau = -1 - -1" >  ${outpath}/${fitprm[$k]}   
        else
            echo -e "WARNING $nl non-zero bins detected: fit not performed\nPrinting error value for the parameters\nA = -2 - -2\nTau = -2 - -2\nSomething has gone wrong" >  ${outpath}/${fitprm[$k]} 
        fi
        if [ "$k" -ne "0" ];then
            fitdata+=" "
        fi  
        fitdata+=$(grep "A " ${outpath}/${fitprm[$k]} | awk '{printf("%.4f %.4f",$3,$5)}')
        fitdata+=$(grep "Tau " ${outpath}/${fitprm[$k]} | awk '{printf(" %.1f %.1f",$3,$5)}')
        
        local AEst=$(grep "A " ${outpath}/${fitprm[$k]} | awk '{printf("%.4f",$3)}')
        local TauEst=$(grep "Tau " ${outpath}/${fitprm[$k]} | awk '{printf("%.1f",$3)}')
        local tauthr=$7
        local taucheck=$(echo "$TauEst" | awk -v tauthr="$tauthr" '{ if ($1>tauthr) print 1; else print 0}') 
        if [ "$taucheck" -eq "1" ]; then #Suspiciously high restime, assume something has gone wrong with the fit
            fitdata+=" -1 20"
        else
        
        
            local test=$(echo "$AEst" | awk '{if ($1 != -1) print 1; else print 0;}')
            if [ "$test" -eq "1" ];then
                fitdata+=$(grep "RMSE " ${outpath}/${fitprm[$k]} | awk '{printf(" %.4f",$3)}')
                #Estimator for discrete distribution
                #no header 
                local estimator=$( echo "" | awk -v ae="$AEst" -v te="$TauEst" -v smpavg="$smpAvg" 'BEGIN{sum=0}{sum+=(ae*$1*exp(-$1/te))}END{print log(sum/smpavg)/log(10)}' ${outpath}/${nzhist} )
                fitdata+=" $estimator"
            else
                #No fit, filler value for RMSE and MagnitudeEstimator 
                fitdata+=" -1 10"
            fi
        fi
        
        local nznldns[$k]="$nldns"
        if [ "$nldns" -ge "4" ];then #Enough non-zero bins to attempt fit
            sedoutpath="${methoddir}/event_domain_${evdmnID}/cntrd_${scid}"
            cp ${blankfitdens} fit_${submt}_dens.py
            
            sed -i "s|OUTPATH|${sedoutpath}|g" fit_${submt}_dens.py
            sed -i "s|HISTO|${nzhistdns}|g" fit_${submt}_dens.py
            sed -i "s|PLOT|${fitpltdns[$k]}|g" fit_${submt}_dens.py
            sed -i "s|FITPRM|${fitprmdns[$k]}|g" fit_${submt}_dens.py
            
            
            python3 fit_${submt}_dens.py > /dev/null 2> error_${submt}_dns.log
        elif [ "$nldns" -ge "1" ]; then
            echo -e "Only $nl non-zero bins detected: fit not performed\nPrinting filler value for the parameters\nA = -1 - -1\nTau = -1 - -1" >  ${outpath}/${fitprmdns[$k]}   

        else
            echo -e "WARNING $nl non-zero bins detected: fit not performed\nPrinting error value for the parameters\nA = -2 - -2\nTau = -2 - -2\nSomething has gone wrong" >  ${outpath}/${fitprmdns[$k]} 
        fi
        if [ "$k" -ne "0" ];then
            fitdatadns+=" "
        fi
        fitdatadns+=$(grep "A " ${outpath}/${fitprmdns[$k]} | awk '{printf("%.4f %.4f",$3,$5)}')
        fitdatadns+=$(grep "Tau " ${outpath}/${fitprmdns[$k]} | awk '{printf(" %.1f %.1f",$3,$5)}')
        

        local AEst=$(grep "A " ${outpath}/${fitprmdns[$k]} | awk '{printf("%.4f",$3)}')
        local TauEst=$(grep "Tau " ${outpath}/${fitprmdns[$k]} | awk '{printf("%.1f",$3)}')
        
        local taucheck=$(echo "$TauEst" | awk -v tauthr="$tauthr" '{ if ($1>tauthr) print 1; else print 0}' )
        if [ "$taucheck" -eq "1" ]; then #Suspiciously high restime, assume something has gone wrong with the fit
            fitdatadns+=" -1 20"
        else
            #Sanity check on fit
            local test=$(echo "$AEst" | awk '{if ($1 != -1) print 1; else print 0;}')
            if [ "$test" -eq "1" ];then
                #Estimator for continuos distribution
                fitdatadns+=$(grep "RMSE " ${outpath}/${fitprmdns[$k]} | awk '{printf(" %.4f" ,$3)}')
                local estimator=$(echo "$AEst $TauEst $smpAvg" | awk '{print log($1*($2^2)/($3))/log(10)}')
                fitdatadns+=" $estimator"
            else
                #No fit, filler value for RMSE and MagnitudeEstimator 
                fitdatadns+=" -1 10"
            fi
        fi

    done
    echo "$fitdata" >> "${outpath}/fitresults.txt"
    echo "$fitdatadns" >> "${outpath}/fitresults_dns.txt"
    echo -e "#SturgesNZBins #DoaneNZBins #FreedmanNZBins\n${nznl[@]}">"${outpath}/nzbins.txt"
    echo -e "#SturgesNZBins #DoaneNZBins #FreedmanNZBins\n${nznldns[@]}">"${outpath}/nzbins_dns.txt"
                
}


