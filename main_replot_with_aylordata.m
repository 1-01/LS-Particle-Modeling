
%Flux calculation, and particles NOT conserved, taken out when they leave
%the model
%Measuring normalized flux of particles at x-gates downwind
%Flux starts ~ 1, and falls to 0 as we go downwind
%There is a weird peak above 1 at first xgate, not sure why

%edit log:


function ls_1_26_21_replot_with_aylordata()

close all
clear all
beep on
tic

%% GET the experi
exp = input('Enter Experiment Number: ');

switch exp
    case 1
        [xmodel_upper ymodel_upper xexp_upper yexp_upper xmodel_lower ymodel_lower xexp_lower yexp_lower] = aylorexp1;
        xlimits = [0 0.6];
    case 2
        [xmodel_upper ymodel_upper xexp_upper yexp_upper xmodel_lower ymodel_lower xexp_lower yexp_lower] = aylorexp2;
        xlimits = [0 0.3];
    case 3
        [xmodel_upper ymodel_upper xexp_upper yexp_upper xmodel_lower ymodel_lower xexp_lower yexp_lower] = aylorexp3;
        xlimits = [0 0.3];
    case 4
        [xmodel_upper ymodel_upper xexp_upper yexp_upper xmodel_lower ymodel_lower xexp_lower yexp_lower] = aylorexp4;
        xlimits = [0 0.3];
    case 5
        [xmodel_upper ymodel_upper xexp_upper yexp_upper xmodel_lower ymodel_lower xexp_lower yexp_lower] = aylorexp5;
        xlimits = [0 0.3];
    case 6
        [xmodel_upper,ymodel_upper,xexp_upper,yexp_upper,xmodel_lower,ymodel_lower,xexp_lower,yexp_lower] = aylorexp6;
        xlimits = [0 0.3];
    otherwise
        error('not an existing exp')
end


[h, d, h0, z0, zw, L, ustar, v_s, k, np, zmin, zmax, tmax, Q_LS, t0, deltstart, sigu_zw, sigw_zw, sigu_h, sigw_h, ex_h, psi_h, Ubar_h, T_Lh, T_L0, beta,x_sense] = param(exp);


% CREATE SPATIAL ARRAYS FOR PARTICLE COUNTS
% xgates and xgate_edges are defined differently than zgates and
% zgate_midpoints because all that matters for xgates is "did a particle
% pass through that particular vertical line, whereas for zgates, the
% vertical space matters.
xgatestart = x_sense; %first xgate
xgatedelta = x_sense/1; %space btw each xgate
xgates = [xgatestart:xgatedelta:x_sense]; %holds the xgates, each index is for the space between the gates including the right edge
xgatenum = length(xgates);
xgate_edges = [xgatestart - xgatedelta/2:xgatedelta:x_sense+xgatedelta/2]; %for pdf histogram functions, so that "xgates" is at the midpoint btw elements in this array
zgates = [zmin:.001:zmax] ; %holds the zgate edges, each index defines the space between the gates including the top edge
zgatenum = length(zgates);
zmidpoints = zgates(1:(zgatenum-1)) + diff(zgates)./2;
zgatedelta = diff(zgates);


runs = 5;      %run the experiment # times to get averages


xexp = xexp_lower;
yexp = yexp_lower;
xmodel = xmodel_lower;
ymodel = ymodel_lower;

plotname = {'lower' 'upper'};



for a = 1:length(h0)   %do lower first, then upper
    C_ls = zeros(runs,zgatenum-1);
    for b = 1:runs
        
        %INITIALIZE CONCENTRATION COUNTS
        count_sum = zeros(xgatenum,zgatenum-1);
        overu_sum = zeros(xgatenum,zgatenum-1);
        u_sum = zeros(xgatenum,zgatenum-1);
        up_sum = zeros(xgatenum,zgatenum-1);
        avgu = zeros(xgatenum,zgatenum-1);
        avgup = zeros(xgatenum,zgatenum-1);
        avgw = zeros(xgatenum,zgatenum-1);
        avgwp = zeros(xgatenum,zgatenum-1);
        adv = 1;
        
        c_raw = zeros(xgatenum,zgatenum-1);
        c_ls = zeros(xgatenum,zgatenum-1);
        in_air = zeros(xgatenum,1);
        pdfx = NaN(xgatenum,1);
        pdfe = NaN(xgatenum,1);
        pdft = NaN(xgatenum,1);
        numpart = 0;
        tops = 0;
        deps = 0;
        ends = 0;
        
        %CREATE ARRAYS FOR PARTICLE TRAJECTORY DATA
        numberofrands = 10^6; 
        M = randn(numberofrands,1); %create large array of random numbers, faster this way
        rnd = 1; %keeps track of index of the random number (in random array M) that we're on. We should increment by 1 every time we use a number from M
        X = zeros(10^4,np);
        Z = zeros(10^4,np);
        T = zeros(10^4,np);
        DT = zeros(10^4,np);
        U = zeros(10^4,np);
        
        %RELEASING ONE PARTICLE AT A TIME LOOP
        for j = 1:np                            %for each particle released
            inair_numpart = zeros(xgatenum,1);  %count of num particles still in air at each xgate, regardless of time
            numpart = numpart+1;                %keeping track of # particles released
            m = 1;                              %first step forward of this particle, increments inside while loop
            z = h0(a);                          %initialize z to release height for this experiment
            x = 0;                              %initialize x to zero
            t0 = t0 + deltstart;                %staggering release time of each particle by deltstart (in 'param.m')
            t = t0;                             %initialize t to t0
            
            
            %GET FIRST WIND STATS
            [sigu sigw dsigu2 dsigw2 uw duw tau Ubar] = windstats(z,L,ustar,z0,k,v_s,h,d,zw, Ubar_h, sigu_h, sigw_h, T_Lh, T_L0, beta, sigu_zw, sigw_zw);
            
            % INITIALIZE VELOCITY FLUCTATIONS (from his code)
            Correl_uw = uw/(sigu*sigw);
            TermCorr_uw = Correl_uw*M(rnd) + (1-Correl_uw^2)^.5*M(rnd+1);
            up = sigu*TermCorr_uw;
            wp = sigw*M(rnd+2);
            rnd = rnd + 3; %increment the random number count because we used 3 of them
            
            % GET FIRST TIME INCREMENT
            dt = calc_dt(Ubar,sigu,tau);
            
            
            % INITIALIZE CONCENTRATION CALCULATIONS
            xprev = x;
            zprev = z;
            uprev = up + Ubar;
            xprev_space = int16(find((xgates-x)>= 0,1)-1);  %find the index of the closest gate to the left of the x location
            if isempty(xprev_space)                             %if its after the last gate (empty), assign it to final gate num
                xprev_space = int16(xgatenum);
            end
            
            % INITIALIZE PARTICLE TRAJECTORY DATA
            X(1,j) = xprev;
            Z(1,j) = zprev;
            U(1,j) = uprev;
            T(1,j) = t;
            DT(1,j) = dt;
            
            % INCREMENTING PARTICLE POSITION LOOP
            while z > zmin && z < zmax && x < (x_sense + 0.5)  %increment particles within these bounds
                m = m+1;
                [sigu sigw dsigu2 dsigw2 uw duw tau Ubar] = windstats(z,L,ustar,z0,k,v_s,h,d,zw, Ubar_h, sigu_h, sigw_h, T_Lh, T_L0, beta, sigu_zw, sigw_zw);
                
                [a_u, a_w, b_u, b_w] = coeffs(sigu, sigw, dsigu2, dsigw2, uw, duw, tau, up, wp);
                
                
                dup = a_u*dt + b_u*M(rnd)*sqrt(dt); %a_u nonturb, b_u turb
                rnd = rnd+1;
                dwp = a_w*dt + b_w*M(rnd)*sqrt(dt); %a_w nonturb, b_w turb
                rnd = rnd+1;
                
                up = up + dup;          %increment horizontal wind vel
                u = up + Ubar;
                wp = wp + dwp;
                
                dx = u*dt;              %delta x
                dz = (wp - v_s)*dt;     %delta z (and he assumes Vavg = 0)
                t = t + dt;             %increment time by dt
                x = x + dx;
                z = z + dz;             %update z position of particl
                
                
                % SAVING particle info in each timestep to plot paths directly
                % at each incremented step "m" and each particle "numpart"
                U(m,numpart) = u;
                X(m,numpart) = x;
                Z(m,numpart) = z;
                T(m,numpart) = t;
                DT(m,numpart) = dt;
                
                
                
                
                
                %% CONCENTRATION CALCULATIONS
                x_space = int16(find((xgates-x)>= 0,1)-1); %find the index of the closest gate to the left of the x location
                if isempty(x_space) %if its after the last gate, assign it to the final gate num
                    x_space = int16(xgatenum);
                end
                
                %COUNT NUMBER OF PARTICLES THAT HIT THE GROUND AND LEAVE
                if z <= zmin                %count number of particles that have been deposited
                    deps = deps + 1;
                    pdfx(deps) = x;
                    if x_space > 0
                        if x - xgates(x_space) < xgatedelta/2  %ie if it's closer to the gate on the left of where it deposits, don't count it as "inair" at that gate. If it's closer to the gate on the right, "array" calc will count it properly.
                            inair_numpart(x_space) = 0;
                        end
                    end
                elseif z >= zmax                %count number of particles that have left top of model
                    tops = tops + 1;
                    pdft(tops) = x;
                elseif x >= (x_sense+.5)        %count number of particles that exit right edge of the model
                    ends = ends + 1;
                    pdfe(ends) = z;
                else                            %else figure out which gates the particle crossed while in air
                    
                    
                    dir = x_space-xprev_space;
                    if dir ~= 0 %if it passed an xgate
                        if dir > 0  %if it moved forward
                            array = xprev_space+int16(1):int16(1):x_space; %account for jumping multiple gates, moving forward, start at gate# to the left of xprev_space
                        else
                            array = xprev_space:int16(-1):x_space+int16(1); %if particle crosses one or more xgates and is moving BACKWARD, end at gate# to the right of x_space
                        end
                        jump = x-xprev; %x distance particle jumped
                        for l = array %cycle through indices of each gate particle passed through (l is the index of xgates)
                            travel = xgates(l) - xprev; %distance between gate that was crossed and previous x value
                            z_cross = (z-zprev)/jump*travel + zprev; %linear interpolation, find z value at which particle crossed x-gate
                            if z_cross >= zgates(1) && z_cross < zgates(zgatenum)   %if it crossed at a height greater than the first zgate and less than the last zgate
                                u_cross = (u-uprev)/jump*travel + uprev; %linear interpolation, find horizontal vel at which particle crossed x-gate
                                z_index = find((zgates-z_cross)>= 0, 1)-1; %index of z midpoint this particle resides in, includes lower gate
                                count_sum(l,z_index) = count_sum(l,z_index) + 1; %sum particles that have crossed this xgate (l) at this particular height (z_index)
                                overu_sum(l,z_index) = overu_sum(l,z_index) + 1/u_cross; %sum inverse of velocities as particle crosses this gate
                            end
                            inair_numpart(l) = 1;   %denoting whether this particle is still in the air at each xgate
                        end
                    end
                    
                    
                    
                    
                    xprev = x;
                    xprev_space = x_space;
                    zprev = z;
                    uprev = u;
                    
                    
                    %get more random numbers if out of them
                    if rnd >= numberofrands-5
                        rnd = 1;
                        M = randn(10^6,1);
                    end
                    
                    %recalculate dt and increment
                    dt = calc_dt(Ubar,sigu,tau);
                    
                end
                
            end
            in_air = in_air + inair_numpart; %summing how many particles are still in the air at each xgate
            
            
        end
        
        c_raw = Q_LS./(np.*zgatedelta).*overu_sum; 
        c_ls = smoothdata(c_raw, 2, 'gaussian',500)*ustar/Q_LS;
        C_ls(b,:) = c_ls;
        
        
    end
    stdev = std(C_ls);
    cvar = var(C_ls);
    avg = mean(C_ls);
    subplot(1,length(h0),a)
    plot(avg, zmidpoints/h,'b','LineWidth',2)
    axis square
    hold on
    plot(avg-stdev, zmidpoints/h,'r','HandleVisibility','off')
    plot(avg+stdev, zmidpoints/h,'r')
    ylim([0 3])
    plot(xmodel,ymodel/h,'-m','LineWidth',2)
    plot(xexp  ,yexp/h, '*-', 'LineWidth',2)
    yline(zw/h, '-', 'roughness sublayer');
    yline(1,'-','canopy top', 'LabelVerticalAlignment','bottom');
    legend('mean concentration','standard dev','aylor model','aylor experim')
    ylabel('Height/canopyheight')
    xlabel('C*ustar/Q_{LS}')
    xlim(xlimits);
    ylim([0 3])
    title(plotname{a})
    %set(gca, 'FontSize',20)
    
    xexp = xexp_upper;
    yexp = yexp_upper;
    xmodel = xmodel_upper;
    ymodel = ymodel_upper;
end


%sgtitle(['Normalized concentrations - Aylor 2001 Lycopodium Experiment ', num2str(exp),'    (' ,num2str(runs), ' runs, ', num2str(np), ' particles each,  x-sense = ',num2str(x_sense), ' m)'])
set(gcf, 'Position',  [100, 100, 1000, 600])
set(findall(gcf,'-property','FontSize'),'FontSize',18)
beep
toc

end

