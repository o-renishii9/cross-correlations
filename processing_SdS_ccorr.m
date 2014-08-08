clear all;
close all;
clc;

%% NOTE TO SELF
% the combinations of models are:
% 1 --> prem vs prem+topo
% 2 --> s20 vs s20+topo
% 3 --> prem+topo vs s20+topo
% 4 --> prem vs s20+topo

%% read all transverse component seismograms
tic
h = waitbar(0,'Please wait . . .');

addpath('/net/home/koroni/tz/processing/seis/prem')
addpath('/net/home/koroni/tz/processing/seis/prem_topo')
addpath('/net/home/koroni/tz/processing/seis/s20')
addpath('/net/home/koroni/tz/processing/seis/s20_topo')

D = dir(['/net/home/koroni/tz/processing/seis/prem/','.*sac']);
D2 = dir(['/net/home/koroni/tz/processing/seis/prem_topo/','*.sac']);
D3 = dir(['/net/home/koroni/tz/processing/seis/s20/','*.sac']);
D4 = dir(['/net/home/koroni/tz/processing/seis/s20_topo/','*.sac']);

v = 6211;
v1 = length(D(not([D.isdir])));
v2 = length(D2(not([D2.isdir])));
v3 = length(D3(not([D3.isdir])));
v4 = length(D4(not([D4.isdir])));

%preallocation for speed

mid_arr_times = zeros(v);

times_ss = zeros(v);
times_s400s = zeros(v);
times_s670s = zeros(v);
dh_pred_410 = rand(v);
dh_pred_670 = rand(v);

for j=1:v;
    if j<=v;
        %% read the data for each model
        
        [t,amp_prem] = readsac(['seis/prem/prem' num2str(j) '.sac']);
        [t,amp_prem_topo] = readsac(['seis/prem_topo/prem_topo' num2str(j) '.sac']);
        [t1,amp_s20] = readsac(['seis/s20/s20' num2str(j) '.sac']);
        [t1,amp_s20_topo] = readsac(['seis/s20_topo/s20_topo' num2str(j) '.sac']);
        
        [prem,npts1,delta1,b1,dist1,az1,baz1,gcarc1]=loadsac(['seis/prem/prem' num2str(j) '.sac']);
        [prem_topo,npts2,delta2,b2,dist2,az2,baz2,gcarc2]=loadsac(['seis/prem_topo/prem_topo' num2str(j) '.sac']);
        [s20,npts,delta3,b3,dist3,az3,baz3,gcarc3]=loadsac(['seis/s20/s20' num2str(j) '.sac']);
        [s20_topo,npts4,delta4,b4,dist4,az4,baz4,gcarc4]=loadsac(['seis/s20_topo/s20_topo' num2str(j) '.sac']);
        
        %% test results with filtered data, if you use that change the name of the input data in windowing section
        %
        %         wn2 = 0.05; % low end of bandpass filter (Hz)/100s (Gu&Dziewonski2002-they use half-welch though)
        %         wn1 = 0.005; % high end of bandpass filter (Hz)/  40s
        %         wn = [wn1,wn2];
        %         wnlp = 0.025;% if you use the lowpass filter
        %
        %         [z,p] = butter(2,wn,'bandpass');
        %         filtp = filter(z,p,prem);
        %         premF = transpose(filtp);
        %
        %         filtpt = filter(z,p,prem_topo);
        %         prem_topoF = transpose(filtpt);
        %
        %         filts20 = filter(z,p,s20);
        %         s20F = transpose(filts20);
        %
        %         filts20t = filter(z,p,s20_topo);
        %         s20_topoF = transpose(filts20t);
        
        %% interpolate  these two data vectors using cubic splines and then
        % extrapolate the deltas for your data points to get the corresponding
        % arrival time
        % This is done for all three phases SS, S400S, S660S
        
        X2= load('deltas_s400s.txt');
        Y2= load('tt_s400s.txt');
        
        X1= load('deltas_ss.txt');
        Y1= load('tt_ss.txt');
        
        X3= load('deltas_s660s.txt');
        Y3= load('tt_s660s.txt');
        
        d_all = load('deltas.dat');
        
        times_ss = interp1(X1,Y1,d_all,'spline');
        times_s400s= interp1(X2,Y2,d_all,'spline');
        times_s660s = interp1(X3,Y3,d_all,'spline');
        
        %% Quality check 1: if the arrival times of s400s and s670s differ
        % by more than 120s then do not use these signals
        
        if times_s400s(j)-times_s660s(j) <= 120;
            
            % lat,long of midpoints
            
            lt1 = load('midlats.dat');
            lg1 = load('midlons.dat');
            len_lt = length(lt1);
            len_lg = length(lg1);
            
            %% window around the precursors, length of window is 100s +-50s from the predicted arrival
            
            dt= abs(t(2)-t(1));    % sampling rate
            dtt1 = abs(t1(2)-t1(1));
            t0 = t(1);              % beginning time
            len = length(t);        % original time array
            len1 = length(t1);
            Fn = 1/(2*dt);          % Nyquist frequency
            
            ts = times_s400s(j)-50; % starting time of window - earliest S670S (i.e for 110deg)
            te = times_s400s(j)+50; % ending time of window   - latest SS (i.e for 170deg)
            
            % find indices for windows
            
            is = (ts-t0)/dt+1;
            ie = (te-t0)/dt+1;
            i1 = floor(is);
            i2 = floor(ie);
            
            wlen1 = i2-i1;	% length of the window
            
            [seis_win_p,win_seis_p] = windoweb(transpose(prem),len,i1,wlen1,10);
            [seis_win_pt,win_seis_pt] = windoweb(transpose(prem_topo),len,i1,wlen1,10);
            [seis_win_s20,win_seis_s20] = windoweb(transpose(s20),len1,i1,wlen1,10);
            [seis_win_s20t,win_seis_s20t] = windoweb(transpose(s20_topo),len1,i1,wlen1,10);
            
            %% compute the distance at each midpoint
            % use then these deltas to interpolate to the taup curves and thus to find
            % the arrival time at that point
            
            % mid_arr_times = interp1(X1,Y1,mid_dist,'spline');
            
            x1 = seis_win_p;      % windowed around S400S signal 1
            x2 = seis_win_pt;     % windowed arounnd S400S signal 2
            x3 = seis_win_s20;
            x4 = seis_win_s20t;
            
            xl1 = length(x1);
            xl2 = length(x2);
            xl3 = length(x3);
            xl4 = length(x4);
         
            % correlation coefficients for all 4 combinations
            
            r_ppt = corrcoef(x1,x2);
            R_ppt(j) = r_ppt(1,2);
            R_ppt = transpose(R_ppt);
            
            r_s20s20t = corrcoef(x3,x4);
            R_s20s20t(j) = r_s20s20t(1,2);
            R_s20s20t = transpose(R_s20s20t);
            
            r_pt_s20t = corrcoef(x2,x4(1:xl2));
            R_pt_s20t(j) = r_pt_s20t(1,2);
            R_pt_s20t=transpose(R_pt_s20t);
            
            r_p_s20t = corrcoef(x1,x4(1:xl1));
            R_p_s20t(j) = r_p_s20t(1,2);
            R_p_s20t=transpose(R_p_s20t);
            
            % gives positions of signals which correlate well (>=0.5)
            suff1 = find(abs(R_ppt)>=0.5);
            suff2 = find(abs(R_s20s20t)>=0.5);
            suff3 = find(abs(R_pt_s20t)>=0.5);
            suff4 = find(abs(R_p_s20t)>=0.5);
            
%% Perform cross correlation under some particular requirements
% NOTE1: WHen the signals are anti-correlated, i.e r<0 then if r>-0.5 take
% estimate the time delay at the lag which corresponds to the minimum correlation value
maxLag400 = round(times_s400s(j)./dt);
maxLag660 = round(times_s660s(j)./dt);

            if abs(R_ppt(j))>= 0.5 ;
                [c_ppt,lags1] = xcorr(x1,x2,maxLag400,'coeff');
                [cmax1,icmax1] = max(abs(c_ppt));
                [cmin1,icmin1] = min(abs(c_ppt));
%                 if R_ppt(j) > 0;
                    icmx1(j) = (icmax1).*dt;
                    delay_eva1= (-1).*lags1(icmax1);
                    actual_del1 = delay_eva1.*dt;
                    dt1(j) = actual_del1;
%                 elseif R_ppt(j) < 0;
%                     icmn1(j) = (icmin1).*dt;
%                     delay_eva1= (1).*lags1(icmin1);
%                     actual_del1 = delay_eva1.*dt;
%                     dt1(j) = actual_del1;
%                 end
            end
            
            if abs(R_s20s20t(j)) >= 0.5;
                [c_s20s20t,lags2] = xcorr(x3,x4,maxLag400,'coeff');
                [cmax2,icmax2] = max(abs(c_s20s20t));
                [cmin2,icmin2] = min(abs(c_s20s20t));
%                 if R_s20s20t(j) > 0 ;
                    icmx2(j) = (icmax2).*dt;
                    delay_eva2= (-1).*lags2(icmax2);
                    actual_del2 = delay_eva2.*dt;
                    dt2(j) = actual_del2;
%                 elseif R_s20s20t(j) < 0 ;
%                     icmn2(j) = (icmin2).*dt;
%                     delay_eva2= (1).*lags2(icmn2);
%                     actual_del2 = delay_eva2.*dt;
%                     dt2(j) = actual_del2;
%                 end
                
            end
            
            if abs(R_pt_s20t(j)) >= 0.5;
                [c_pt_s20t,lags3] = xcorr(x2,x4(1:xl2),maxLag400,'coeff');
                [cmax3,icmax3] = max(abs(c_pt_s20t));
                [cmin3,icmin3] = min(abs(c_pt_s20t));
%                 if R_pt_s20t(j) < 0 ;
%                     icmn3(j) = (icmin3).*dt;
%                     delay_eva3= (1).*lags3(icmin3);
%                     actual_del3 = delay_eva3.*dt;
%                     dt3(j) = actual_del3;
%                 elseif R_pt_s20t(j) > 0;
                    icmx3(j) = (icmax3).*dt;
                    delay_eva3= (-1).*lags3(icmax3);
                    actual_del3 = delay_eva3.*dt;
                    dt3(j) = actual_del3;
%                 end
            end
            
            if abs(R_p_s20t(j)) >= 0.5;
                [c_p_s20t,lags4] = xcorr(x1,x4(1:xl1),maxLag400,'coeff');
                [cmax4,icmax4] = max(abs(c_p_s20t));
                [cmin4,icmin4] = min(abs(c_p_s20t));
%                 if R_p_s20t(j) < 0;
%                     icmn4(j) = (icmin4).*dt;
%                     delay_eva4= (1).*lags4(icmin4);
%                     actual_del4 = delay_eva4.*dt;
%                     dt4(j) = actual_del4;
%                 elseif R_p_s20t(j) > 0;
                    icmx4(j) = (icmax4).*dt;
                    delay_eva4= (-1).*lags4(icmax4);
                    actual_del4 = delay_eva4.*dt;
                    dt4(j) = actual_del4;
%                 end
            end
            
            
            %% same for S660S
            
            ts660 = times_s660s(j)-50; % starting time of window
            te660 = times_s660s(j)+50; % ending time of window
            
            % find indices for windows
            
            is660 = (ts660-t0)/dt+1;
            ie660 = (te660-t0)/dt+1;
            i1_660 = floor(is660);
            i2_660 = floor(ie660);
            
            wlen1_660 = i2_660-i1_660;	% length of the window
            
            [seis_win_p660,win_seis_p] = windoweb(transpose(prem),len,i1_660,wlen1_660,10);
            [seis_win_pt660,win_seis_pt] = windoweb(transpose(prem_topo),len,i1_660,wlen1_660,10);
            [seis_win_s20660,win_seis_s20] = windoweb(transpose(s20),len1,i1_660,wlen1_660,10);
            [seis_win_s20t660,win_seis_s20t] = windoweb(transpose(s20_topo),len1,i1_660,wlen1_660,10);
            
            %% compute the distance at each midpoint
            % use then these deltas to interpolate to the taup curves and thus to find
            % the arrival time at that point
            
            % mid_arr_times = interp1(X1,Y1,mid_dist,'spline');
            
            x11 = seis_win_p660;      % windowed around S400S signal 1
            x22 = seis_win_pt660; % windowed arounnd S400S signal 2
            x33 = seis_win_s20660;
            x44 = seis_win_s20t660;
            
            xl11 = length(x11);
            xl22 = length(x22);
            xl33 = length(x33);
            xl44 = length(x44);
            
            %         fx1 = fft(x11);
            %         fx2 = fft(x22);
            %         fx3 = fft(x33);
            %         fx4 = fft(x44);
            
            % correlation coefficients for all 4 combinations
            
            r_ppt660 = corrcoef(x11,x22);
            R_ppt660(j) = r_ppt660(1,2);
            
            r_s20s20t660 = corrcoef(x33,x44);
            R_s20s20t660(j) = r_s20s20t660(1,2);
            
            r_pt_s20t660 = corrcoef(x22,x44(1:xl22));
            R_pt_s20t660(j) = r_pt_s20t660(1,2);
            
            r_p_s20t660 = corrcoef(x11,x44(1:xl11));
            R_p_s20t660(j) = r_p_s20t660(1,2);
            
            % Quality check 1: keep only correlations higher or equal to 0.5
            
            suff11 = find(abs(R_ppt660)>=0.5);
            suff22 = find(abs(R_s20s20t660)>=0.5);
            suff33 = find(abs(R_pt_s20t660)>=0.5);
            suff44 = find(abs(R_p_s20t660)>=0.5);
            
            if abs(R_ppt660(j)) >= 0.5 ;
                
                [c_ppt660,lags11] = xcorr(x11,x22,maxLag660,'coeff');
                [cmax11,icmax11]= max(abs(c_ppt660));
                [cmin11,icmin11] = min(abs(c_ppt660));
%                 if R_ppt660(j) > 0;
                    icmx11(j) =  (icmax11).*dt;
                    delay_eva11= lags11(icmax11);
                    actual_del11 = (-1).*delay_eva11.*dt;
                    dt11(j) = actual_del11;
%                 elseif  R_ppt660(j) < 0;
%                     icmn11(j) =  (icmin11).*dt;
%                     delay_eva11= (1).*lags11(icmin11);
%                     actual_del11 = delay_eva11.*dt;
%                     dt11(j) = actual_del11;
%                 end
            end
            
            if abs(R_s20s20t660(j)) >= 0.5;
                
                [c_s20s20t660,lags22] = xcorr(x33,x44,maxLag660,'coeff');
                
                [cmax22,icmax22] = max(abs(c_s20s20t660));
                [cmin22,icmin22] = min(abs(c_s20s20t660));
%                 if R_s20s20t660(j) > 0;
                    icmx22(j) =  (icmax22).*dt;
                    delay_eva22= (-1).*lags22(icmax22);
                    actual_del22 = delay_eva22.*dt;
                    dt22(j) = actual_del22;
%                 elseif R_s20s20t660(j) < 0;
%                     icmn22(j) =  (icmin22).*dt;
%                     delay_eva22= (1).*lags22(icmin22);
%                     actual_del22 = delay_eva22.*dt;
%                     dt22(j) = actual_del22;
%                 end
            end
            
            if abs(R_pt_s20t660(j)) >= 0.5;
                [c_pt_s20t660,lags33] = xcorr(x22,x44(1:xl22),maxLag660,'coeff');
                [cmax33,icmax33] = max(abs(c_pt_s20t660));
                [cmin33,icmin33] = min(abs(c_pt_s20t660));
%                 if R_pt_s20t660(j) > 0;
                    icmx33(j) =  (icmax33).*dt;
                    delay_eva33= (-1).*lags33(icmax33);
                    actual_del33 = delay_eva33.*dt;
                    dt33(j) = actual_del33;
%                 elseif R_pt_s20t660(j) < 0;
%                     icmn33(j) =  (icmin33).*dt;
%                     delay_eva33= (1).*lags33(icmin33);
%                     actual_del33 = delay_eva33.*dt;
%                     dt33(j) = actual_del33;
%                 end
            end
            
            
            if abs(R_p_s20t660(j)) >= 0.5;
                [c_p_s20t660,lags44] = xcorr(x11,x44(1:xl11),maxLag660,'coeff');
                [cmax44,icmax44] = max(abs(c_p_s20t660));
                [cmin44,icmin44] = min(abs(c_p_s20t660));
%                 if R_p_s20t660(j) > 0;
                    icmx44(j) = (icmax44).*dt;
                    delay_eva44= (-1).*lags44(icmax44);
                    actual_del44 = delay_eva44.*dt;
                    dt44(j) = actual_del44;
%                 elseif R_p_s20t660(j) < 0;
%                     icmn44(j) = (icmin44).*dt;
%                     delay_eva44= (1).*lags44(icmin44);
%                     actual_del44 = delay_eva44.*dt;
%                     dt44(j) = actual_del44;
%                 end
            end
            
            
            %%
            % compute the predicted (i.e according to ray theory) travel time differences
            % this won't change for all regressions
            Vs_400 = 4.76; % km/s
            Vs_660 = 5.94; % km/s
            
            dh_400 = load('dh400.dat');
            dh_660 = load('dh660.dat');
            
            dt_pred_400(j) = ((-2)*dh_400(j))./Vs_400; % s
            dt_pred_660(j) = ((-2)*dh_660(j))./Vs_660; % s
            
        end
    end
        waitbar(j/10000,h)
end

delays{:,1:8} = {dt1./1000,dt2,dt3./1000,dt4./1000,dt11,dt22,dt33,dt44};
% 400
dt1_indx = find(dt1);
dt2_indx = find(dt2);
dt3_indx = find(dt4);
dt4_indx = find(dt4);
% 660
dt11_indx = find(dt11);
dt22_indx = find(dt22);
dt33_indx = find(dt33);
dt44_indx = find(dt44);
delays_indx{:,1:8} = {dt1_indx,dt2_indx,dt3_indx,dt4_indx,dt11_indx,dt22_indx,dt33_indx,dt44_indx};

%% perform linear regression to compare predicted and measured by cross correlation travel times
% note that this has to be adopted for all 4 combinations of models
% this is just a preface

% preallocate x,y
% x400 = zeros(6211);
% y400 = zeros(6211);
% x660 = zeros(6211);
% y660 = zeros(v);
for k=1:4; % 4 combinations of models used now
    if k<=4;
        %         y400 = zeros(length(delays_indx{:,1}{:,k}));
        
             
        % 400 km discontinuity
         for x=1:length(delays_indx{:,1}{:,k});

           l = delays_indx{:,1}{:,k};
           x400 = (delays{:,1}{:,k}(l));
           y400 = dt_pred_400(l);
        end
        
        coeff400 = polyfit(x400,y400,1);
        P400 = polyval(coeff400,x400);
        % plot data
        figure()
        plot(x400,y400,'*',x400,P400,':') ;
        xlabel('measured travel time differences');
        ylabel('predicted travelt ime differences');
        title('Linear regreesion of predicted and measured dt (400 km)');
        
    end
    
end
for kk= 5:8;
    if kk<=8;
        % 660 km discontinuity
        for xx=1:length(delays_indx{:,1}{:,kk});
            
            ll = delays_indx{:,1}{:,kk};
            x660 = (delays{:,1}{:,kk}(ll));
            y660 = dt_pred_660(ll);
        end
        coeff660 = polyfit(x660,y660,1);
        P660 = polyval(coeff660,x660)
        % plot data
        figure()
        plot(x660,y660,'*',x660,P660,':') 
        xlabel('measured travel time differences');
        ylabel('predicted travelt ime differences');
        title('Linear regreesion of predicted and measured dt (660 km)');
        
        
    end
    
end





close(h)
toc
