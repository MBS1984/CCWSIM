
function [C] = CCWSIM_2D(ti, hd, T, OL, CT, prop, cand, fc, mrp, T_vibration, wavelet_level, real_numb)
SG_size = size(hd);
%------------- GENERATE MULTIPLE RANDOM PATH FLAGS ------------------------
HD1_1 = hd; ti1_1 = ti; lag = 1;
if mrp == 1
    lag = 8;
    HD2_1 = flip(HD1_1,1); ti2_1 = flip(ti1_1,1);
    HD3_1 = flip(HD1_1,2); ti3_1 = flip(ti1_1,2);
    HD4_1 = flip(flip(HD1_1,1),2); ti4_1 = flip(flip(ti1_1,1),2);
    HD5_1 = HD1_1'; ti5_1 = ti1_1';
    HD6_1 = HD3_1'; ti6_1 = ti3_1';
    HD7_1 = HD2_1'; ti7_1 = ti2_1';
    HD8_1 = HD4_1'; ti8_1 = ti4_1';
end

C = zeros(numel(hd),1);

h = waitbar(0, sprintf('CCWSIM is running ... %i realization(s)', real_numb), ...
    'CreateCancelBtn', 'setappdata(gcbf,''canceling'',1)');
setappdata(h,'canceling',0)

R1 = 1:lag:real_numb;
if mrp == 1
    R2 = 2:lag:real_numb;
    R3 = 3:lag:real_numb;
    R4 = 4:lag:real_numb;
    R5 = 5:lag:real_numb;
    R6 = 6:lag:real_numb;
    R7 = 7:lag:real_numb;
    R8 = 8:lag:real_numb;
end
if (real_numb > 1) && (T_vibration == 1)
    T_new = ones(real_numb,1);
    T_new(1:ceil(real_numb/3),1) = T;
    T_new(1+ceil(real_numb/3):2*ceil(real_numb/3),1) = T + 4;
    T_new(1+2*ceil(real_numb/3):end,1) = T - 4;
    T_new(:,2) = rand(size(T_new,1),1);
    T_new = sortrows(T_new,2);
    T_new = T_new(:,1);
else
    T_new = repmat(T,[real_numb 1]);
end
error_location = zeros(numel(find(~isnan(hd))),real_numb);

% Initialize an array to store the tEnd values
tEndArray = zeros(real_numb, 1);
% Initialize an array to store the percentage of mismatched HD

misMatched = zeros(real_numb, 1);


%------------INITIATE THE SIMULATION USING THE CCWSIM ALGORITHM -----------
for i = 1:real_numb
    if any(R1==i)
        hd0 = HD1_1; ti0 = ti1_1;
        
    elseif any(R2==i)
        hd0 = HD2_1; ti0 = ti2_1;
        
    elseif any(R3==i)
        hd0 = HD3_1; ti0 = ti3_1;
        
    elseif any(R4==i)
        hd0 = HD4_1; ti0 = ti4_1;
        
    elseif any(R5==i)
        hd0 = HD5_1; ti0 = ti5_1;
        
    elseif any(R6==i)
        hd0 = HD6_1; ti0 = ti6_1;
        
    elseif any(R7==i)
        hd0 = HD7_1; ti0 = ti7_1;
        
    elseif any(R8==i)
        hd0 = HD8_1; ti0 = ti8_1;
        
    end
    
    tStart = tic;
    
    % Resize conditional data based on DWT level
    hd0 = hd_resize_2D(hd0,[size(hd,1)/2.^wavelet_level  size(hd,2)/2.^wavelet_level]);
%     TI=ti0;
    
    % simulate based on DWT level
    [Grid_Sim] = CCWSIM_main(ti0, hd0, T, OL, CT, fc, prop, cand, wavelet_level, SG_size) ;
    
    tEnd = toc(tStart);
    % ---------------  END OF THE SIMULATION  -----------------------------
    
    if any(R2==i)
        Grid_Sim = flip(Grid_Sim,1);
        
    elseif any(R3==i)
        Grid_Sim = flip(Grid_Sim,2);
        
    elseif any(R4==i)
        Grid_Sim = flip(flip(Grid_Sim,1),2);
        
    elseif any(R5==i)
        Grid_Sim = Grid_Sim';
        
    elseif any(R6==i)
        Grid_Sim = Grid_Sim';
        Grid_Sim = flip(Grid_Sim,2);
        
    elseif any(R7==i)
        Grid_Sim = Grid_Sim';
        Grid_Sim = flip(Grid_Sim,1);
        
    elseif any(R8==i)
        Grid_Sim = Grid_Sim';
        Grid_Sim = flip(flip(Grid_Sim,1),2);
    else
        Grid_Sim = Grid_Sim(1:size(hd,1),1:size(hd,2));
    end
    
    C(:,i) = Grid_Sim(:);
   
    waitbar(i / real_numb, h, sprintf('CCWSIM is running...Please wait...%ith realization completed', i))
    
    disp(['********  CPU time for the grid size of ',num2str(size(hd,1)),'x',...
        num2str(size(hd,2)), ' is ', num2str(tEnd),...
        '  (s) ********'])
    tEndArray(i) = tEnd;  % Store tEnd for every realization
    % Save tEndArray to the workspace
    assignin('base', 'tEndArray', tEndArray);
    meanTime = mean(tEndArray,1);
    assignin('base', 'meanTime', meanTime);
    
    % Show the simulation result
    subplot(1,2,1); imagesc(ti);title TI; axis equal tight xy;
    colormap copper
    subplot(1,2,2); imagesc(Grid_Sim); title CCWSIM; axis equal tight xy;
    colormap copper
    hold on
    [zero_element, non_zero_element] = Display_cnd_data(hd);
end
delete(h)
end