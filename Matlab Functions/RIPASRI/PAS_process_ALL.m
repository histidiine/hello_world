function PAS_routine

gp_ID=0;ini=0;
while strcmp(gp_ID,'HV')==0 && strcmp(gp_ID,'CRPS')==0 && strcmp(gp_ID,'CRPSD')==0 && strcmp(gp_ID,'FHD')==0 && strcmp(gp_ID,'FUNCT')==0
    if ini == 1
        disp('SORRY ! That is not a valid group name, check out what I suggested in the parenthesis down here v')
    end
    gp_ID = input('Please enter the code for the group you want to process into single quotes (HV,CRPS,CRPSD,FHD,FUNCT):   ');
    ini=1;
end


[data_goodrecs, resp] = PAS_process_MAR2017(gp_ID,sub_ID);


end