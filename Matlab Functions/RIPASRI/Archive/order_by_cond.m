function output=order_by_cond(data_matrix,log,num)
% log = conditions, num = number of trials

count1=1;count2=1;count3=1;count4=1;count5=1;count6=1;count7=1;count8=1;count9=1;count10=1;count11=1;count12=1;count13=1;count14=1;count15=1;
for i = 1:num
    if log(i,3)==10
        if log(i,5)==0
            FCR_cond_ms(:,count3,3) =  data_matrix(:,i);
            count3=count3+1;
        else if log(i,5)==-0.5
                FCR_cond_ms(:,count2,2) =  data_matrix(:,i);
                count2=count2+1;
            else if log(i,5)==-1
                    FCR_cond_ms(:,count1,1) =  data_matrix(:,i);
                    count1=count1+1;
                end
            end
        end
    else if log(i,3)==200
            FCR_cond_ms(:,count15,15) =  data_matrix(:,i);
            count15=count15+1;
        else if log(i,5)==-0.5
                FCR_cond_ms(:,count4,4) =  data_matrix(:,i);
                count4=count4+1;
            else if log(i,5)==-1
                    FCR_cond_ms(:,count5,5) =  data_matrix(:,i);
                    count5=count5+1;
                else if log(i,5)==-2
                        FCR_cond_ms(:,count6,6) =  data_matrix(:,i);
                        count6=count6+1;
                    else if log(i,5)==-5
                            FCR_cond_ms(:,count7,7) =  data_matrix(:,i);
                            count7=count7+1;
                        else if log(i,5)==-7.5
                                FCR_cond_ms(:,count8,8) =  data_matrix(:,i);
                                count8=count8+1;
                            else if log(i,5)==-10
                                    FCR_cond_ms(:,count9,9) =  data_matrix(:,i);
                                    count9=count9+1;
                                else if log(i,5)==-25
                                        FCR_cond_ms(:,count10,10) =  data_matrix(:,i);
                                        count10=count10+1;
                                    else if log(i,5)==-50
                                            FCR_cond_ms(:,count11,11) =  data_matrix(:,i);
                                            count11=count11+1;
                                        else if log(i,5)==-75
                                                FCR_cond_ms(:,count12,12) =  data_matrix(:,i);
                                                count12=count12+1;
                                            else if log(i,5)==-100
                                                    FCR_cond_ms(:,count13,13) =  data_matrix(:,i);
                                                    count13=count13+1;
                                                else if log(i,5)==-200
                                                        FCR_cond_ms(:,count14,14) =  data_matrix(:,i);
                                                        count14=count14+1;
                                                        
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
output=FCR_cond_ms;
end