function data_cond=order_by_cond_JUL2017(data_matrix,conditions,n_trial)

count1=1;count2=1;count3=1;count4=1;count5=1;count6=1;count7=1;count8=1;count9=1;count10=1;count11=1;count12=1;count13=1;count14=1;count15=1;
for j = 1:2 %step (baseline, post)
    count1=1;count2=1;count3=1;count4=1;count5=1;count6=1;count7=1;count8=1;...
        count9=1;count10=1;count11=1;count12=1;count13=1;count14=1;count15=1;
    temp_data = squeeze(data_matrix(:,:,:,j)); % timepoints, muscle, trials [5000, 2, 10]
    for i = 1:n_trial %trials
        if conditions(i,1,j)==10
            if conditions(i,2,j)==0
                data_cond(:,:,count4,4,j) =  temp_data(:,:,i); % timepoints, muscle, trials (per cond), cond, step
                count4=count4+1;
            else if conditions(i,2,j)==-0.5
                    data_cond(:,:,count3,3,j) =  temp_data(:,:,i);
                    count3=count3+1;
                else if conditions(i,2,j)==-1
                        data_cond(:,:,count2,2,j) =  temp_data(:,:,i);
                        count2=count2+1;
                    end
                end
            end
        else if conditions(i,1,j)==200 % SINGLE PULSE
                data_cond(:,:,count1,1,j) =  temp_data(:,:,i);
                count1=count1+1;
            else if conditions(i,2,j)==-0.5
                    data_cond(:,:,count5,5,j) =  temp_data(:,:,i);
                    count5=count5+1;
                else if conditions(i,2,j)==-1
                        data_cond(:,:,count6,6,j) =  temp_data(:,:,i);
                        count6=count6+1;
                    else if conditions(i,2,j)==-2
                            data_cond(:,:,count7,7,j) =  temp_data(:,:,i);
                            count7=count7+1;
                        else if conditions(i,2,j)==-5
                                data_cond(:,:,count8,8,j) =  temp_data(:,:,i);
                                count8=count8+1;
                            else if conditions(i,2,j)==-7.5
                                    data_cond(:,:,count9,9,j) =  temp_data(:,:,i);
                                    count9=count9+1;
                                else if conditions(i,2,j)==-10
                                        data_cond(:,:,count10,10,j) =  temp_data(:,:,i);
                                        count10=count10+1;
                                    else if conditions(i,2,j)==-25
                                            data_cond(:,:,count11,11,j) =  temp_data(:,:,i);
                                            count11=count11+1;
                                        else if conditions(i,2,j)==-50
                                                data_cond(:,:,count12,12,j) =  temp_data(:,:,i);
                                                count12=count12+1;
                                            else if conditions(i,2,j)==-75
                                                    data_cond(:,:,count13,13,j) =  temp_data(:,:,i);
                                                    count13=count13+1;
                                                else if conditions(i,2,j)==-100
                                                        data_cond(:,:,count14,14,j) =  temp_data(:,:,i);
                                                        count14=count14+1;
                                                    else if conditions(i,2,j)==-200
                                                            data_cond(:,:,count15,15,j) =  temp_data(:,:,i);
                                                            count15=count15+1;                                                            
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                
            end
        end
    end
    clear count1 count2 count3 count4 count5 count6 count7 count8 count9 count10 count11 count12 count13 count14 count15
end

end