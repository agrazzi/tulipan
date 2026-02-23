#!/bin/bash
function histo_fit {
    local submt=$1
    local evdmnID=$2
    local scid=$3
    local supraclfile=$4
    local outpath=$5
    local blank_mlthst=$6
    #Compute simple average and trimmed average (skip top and bottom 5% of entries)
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
                local estimator=$(echo "$AEst $TauEst $smpAvg" | awk '{print log($1*($2**2)/($3))/log(10)}')
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



