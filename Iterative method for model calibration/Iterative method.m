clc
clear all
close all
S_1=[9.69474,7.92501,7.331,6.91571,6.69071];
S_2=[8.41,8.14,7.4,6.24,5.02];
S_3=[0.81,1.18,2.14,3.45,4.86];
S_4=[0.9193,1.00432,1.23303,1.32664,1.36613];
[S1,S2,S3,S4] = f(0.2,2,100,0.5,5,0.6,0.2,0.85,0.03,0.09);
error1 = sum(sqrt(S1-S_1))+sum(sqrt(S2-S_2))+sum(sqrt(S3-S_3))+sum(sqrt(S4-S_4));
New1 = S1;
New2 = S2;
New3 = S3;
New4 = S4;
for K_O = 0.05:0.01:5
    for K_S = 0.5:0.1:5
        for K_UAP = 50:10:300
            for K_NO = 0.1:0.1:1
                for k_STO = 0.5:0.5:10
                    for n_NO = 0.2:0.1:5
                        for b_HO = 0.1:0.1:1
                            for Y_STOO = 0.5:0.5:5
                                for i_NSS = 0.01:0.01:1
                                    for k_BAPO = 0.01:0.01:2
                                        [New_S_1,New_S_2,New_S_3,New_S_4]=f(K_O,K_S,K_UAP,K_NO,k_STO,n_NO,b_HO,Y_STOO,i_NSS,k_BAPO);
                                        error2 = sum(sqrt(New_S_1-S_1))+sum(sqrt(New_S_2-S_2))+sum(sqrt(New_S_3-S_3))+sum(sqrt(New_S_4-S_4));
                                        if abs(error2)<abs(error1)
                                            New1 = New_S_1;
                                            New2 = New_S_2;
                                            New3 = New_S_3;
                                            New4 = New_S_4;
                                            error1 = error2;
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




