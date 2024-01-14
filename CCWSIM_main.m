function [Grid_Sim] = CCWSIM_main(TI, hd, T , OL, CT, fc, prop, cand, wavelet_level, SG_size)

Grid_Sim = NaN(SG_size);
sizeout = [SG_size SG_size];

%----------- Calculating the Discrete Wavelet Transform (DWT) for TI ------
% The Haar mother wavelet is employed for decomposition

[A1_TI,H1_TI,V1_TI,D1_TI] = dwt2(TI,'haar');
[c,s]= wavedec2(TI,wavelet_level,'haar');
A3_TI = appcoef2(c,s,'haar',wavelet_level);

temp = ones([OL/(2^wavelet_level) T/(2^wavelet_level)]);
CCtop = xcorr2(A3_TI.^2, temp);

temp = ones([T/(2^wavelet_level) OL/(2^wavelet_level)]);
CCside = xcorr2(A3_TI.^2,temp);

temp = ones([(T-OL)/(2^wavelet_level) OL/(2^wavelet_level)]);
CCsidesmall = xcorr2(A3_TI.^2,temp);

cntr = 0;
a = numel(1:T-OL:Grid_Sim(1)-T+1);
b = numel(1:T-OL:Grid_Sim(2)-T+1);
LOC = NaN(a*b,2);

for i=[1:T-OL:sizeout(1)-T, sizeout(1)-T+1]
    for j=[1:T-OL:sizeout(2)-T, sizeout(2)-T+1]
        cntr = cntr+1;
        selectrows = ceil(i/(2^wavelet_level)):min(size(hd,1),ceil(i/(2^wavelet_level))+CT(1)*T/(2^wavelet_level)-1);
        selectcols = ceil(j/(2^wavelet_level)):min(size(hd,2),ceil(j/(2^wavelet_level))+CT(2)*T/(2^wavelet_level)-1);
        
        hd_dev = hd(ceil(i/(2^wavelet_level)):ceil(i/(2^wavelet_level))+T/(2^wavelet_level)-1,ceil(j/(2^wavelet_level)):ceil(j/(2^wavelet_level))+T/(2^wavelet_level)-1);
        hd_dev = reshape(hd_dev,1,size(hd_dev,1)*size(hd_dev,2));
        hd_indicator = sum(isfinite(hd_dev));
        hd_index = find(~isnan(hd_dev));
        
        HD_dev = hd(selectrows, selectcols);
        HD_dev(1:T/(2^wavelet_level),1:T/(2^wavelet_level)) = NaN;
        HD_indicator = sum(isfinite(HD_dev(:)));
        
        dev_init = Grid_Sim(i:i+T-1, j:j+T-1);
        
        if (i > 1) && (j > 1)
            shared = Grid_Sim(i:i+OL-1,j:j+T-1);
            [c,s] = wavedec2(shared,wavelet_level,'haar');
            A3_shared = appcoef2(c,s,'haar',wavelet_level);
            CC = CCtop - 2 * xcorr2(A3_TI, A3_shared) + sum(A3_shared(:).^2);
            CC = CC(OL/2^wavelet_level:end-(T/2^wavelet_level)+1,T/2^wavelet_level:end-T/2^wavelet_level+1);
            
            shared = Grid_Sim(i+OL:i+T-1,j:j+OL-1);
            [c,s] = wavedec2(shared,wavelet_level,'haar');
            A3_shared = appcoef2(c,s,'haar',wavelet_level);
            CC2 = CCsidesmall - 2 * xcorr2(A3_TI, A3_shared) + sum(A3_shared(:).^2);
            CC = CC + CC2(T/(2^wavelet_level):end-(T/2^wavelet_level)+(OL/2^wavelet_level)+1, OL/2^wavelet_level:end-(T/2^wavelet_level)+1);
            
        elseif i > 1
            shared = Grid_Sim(i:i+OL-1,j:j+T-1);
            [c,s] = wavedec2(shared,wavelet_level,'haar');
            A3_shared = appcoef2(c,s,'haar',wavelet_level);
            CC = CCtop - 2 * xcorr2(A3_TI, A3_shared) + sum(A3_shared(:).^2);
            CC = CC(OL/(2^wavelet_level):end-T/(2^wavelet_level)+1,T/(2^wavelet_level):end-T/(2^wavelet_level)+1);
            
        elseif j > 1
            shared = Grid_Sim(i:i+T-1,j:j+OL-1);
            [c,s] = wavedec2(shared,wavelet_level,'haar');
            A3_shared = appcoef2(c,s,'haar',wavelet_level);
            CC = CCside - 2 * xcorr2(A3_TI, A3_shared) + sum(A3_shared(:).^2);
            CC = CC(T/(2^wavelet_level):end-T/(2^wavelet_level)+1,OL/(2^wavelet_level):end-T/(2^wavelet_level)+1);
        else
            CC = rand(size(A3_TI,1)-T/(2^wavelet_level),size(A3_TI,2)-T/(2^wavelet_level));
        end
        
        %----------  UNCONDITIONAL SIMULATION  ----------------------------
        if hd_indicator==0 && HD_indicator==0
            [~,loc] = sort(CC(:));
            [ibest, jbest] = ind2sub(size(CC),loc(1:cand,1));
            if fc~=0
                [c, ~] = hist_cat(TI, Grid_Sim, T, OL, fc, ibest, jbest, i, j);
            else
                c = ceil(rand * length(ibest));
            end
            pos = [ibest(c) jbest(c)];
            X_final = pos(1); Y_final = pos(2);
        else
            %------------- CONDITIONAL SIMULATION  ------------------------
            [~, loc] = sort(CC(:));
            [loc(:,2),loc(:,3)] = ind2sub(size(CC),loc(:));
            now=0; difh=1E+5; difH=1E+5; DifT=1E+5;
            
            while (DifT~=0) && (now+1<=ceil(prop*size(loc,1)))
                now = now+1;
                x = loc(now,2);
                y = loc(now,3);
                 
                target = A3_TI(x:x+T/(2^wavelet_level)-1,y:y+T/(2^wavelet_level)-1);
                target = round(target);
                
                target(target < 2^(wavelet_level)/2) = 0;
                target(target > 2^(wavelet_level)/2 & target <= 2^(wavelet_level)) = 1;
                target(target > 2^(wavelet_level) & target <= 2*2^(wavelet_level)) = 2;

                X = x:min(size(TI,1)/(2^wavelet_level),x+CT(1)*T/(2^wavelet_level)-1);
                Y = y:min(size(TI,2)/(2^wavelet_level),y+CT(2)*T/(2^wavelet_level)-1);
                
                TARGET = A3_TI(X, Y);
                TARGET = round(TARGET);
                
                TARGET(TARGET <(2^(wavelet_level)/2)) = 0;
                TARGET(TARGET > 2^(wavelet_level)/2 & TARGET <= 2^(wavelet_level)) = 1;
                TARGET(TARGET > 2^(wavelet_level) & TARGET <= 2*2^(wavelet_level)) = 2;
            
               
                HD_dev2 = HD_dev(1:min(size(TARGET,1),size(HD_dev,1)),1:min(size(TARGET,2),size(HD_dev,2)));
                HD_dev2(1:T/(2^wavelet_level),1:T/(2^wavelet_level)) = NaN;
                HD_dev2 = reshape(HD_dev2,1,size(HD_dev2,1)*size(HD_dev2,2));
                HD_index2 = find(~isnan(HD_dev2));
                
                Difh = sum(abs(hd_dev(hd_index)-target(hd_index)));
                DifH = sum(abs(HD_dev2(HD_index2)-TARGET(HD_index2)));
                DifT = sum(Difh + DifH);
                
                if Difh <= difh
                    difh = Difh;
                    X_final = x; Y_final = y;
                    if Difh < difh
                        difH = 1E+5;
                    end
                    if DifH < difH
                        difH = DifH;
                        X_final = x; Y_final = y;
                    end
                end
            end
        end
        
        LOC (cntr,1)= X_final; LOC (cntr,2)= Y_final;
        
        Target_A = A1_TI((2^wavelet_level)/2*X_final-(2^wavelet_level/2-1):...
            (2^wavelet_level)/2*X_final+T/2-(2^wavelet_level)/2,(2^wavelet_level)/2*Y_final...
            -(2^wavelet_level/2-1):(2^wavelet_level)/2*Y_final+T/2-(2^wavelet_level)/2);
        Target_H = H1_TI((2^wavelet_level)/2*X_final-(2^wavelet_level/2-1):...
            (2^wavelet_level)/2*X_final+T/2-(2^wavelet_level)/2,(2^wavelet_level)/2*Y_final...
            -(2^wavelet_level/2-1):(2^wavelet_level)/2*Y_final+T/2-(2^wavelet_level)/2);
        Target_V = V1_TI((2^wavelet_level)/2*X_final-(2^wavelet_level/2-1):...
            (2^wavelet_level)/2*X_final+T/2-(2^wavelet_level)/2,(2^wavelet_level)/2*Y_final...
            -(2^wavelet_level/2-1):(2^wavelet_level)/2*Y_final+T/2-(2^wavelet_level)/2);
        Target_D = D1_TI((2^wavelet_level)/2*X_final-(2^wavelet_level/2-1):...
            (2^wavelet_level)/2*X_final+T/2-(2^wavelet_level)/2,(2^wavelet_level)/2*Y_final...
            -(2^wavelet_level/2-1):(2^wavelet_level)/2*Y_final+T/2-(2^wavelet_level)/2);
        
        %---------------- Calculation of Inverse-DWT-----------------------
        Target = idwt2(Target_A,Target_H,Target_V,Target_D,'haar');
        
        M = mincut_func(Target, Grid_Sim(i:i+T-1, j:j+T-1), T, OL, i, j);
        
        dev_index = isnan(dev_init);
        dev_init (dev_index) = 0   ;
        
        Grid_Sim(i:i+T-1,j:j+T-1) = combine_2D(dev_init,Target, M);
    end
end