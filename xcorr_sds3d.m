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

addpath('/net/home/koroni/tz/processing/seis/s20')
addpath('/net/home/koroni/tz/processing/seis/s20_topo')

D3 = dir(['/net/home/koroni/tz/processing/seis/s20/','*.sac']);
D4 = dir(['/net/home/koroni/tz/processing/seis/s20_topo/','*.sac']);

v = 6211;

v3 = length(D3(not([D3.isdir])));
v4 = length(D4(not([D4.isdir])));

%preallocation for speed

% mid_arr_times = zeros(v);
%
% times_ss = zeros(v);
% times_s400s = zeros(v);
% times_s670s = zeros(v);
% dh_pred_410 = rand(v);
% dh_pred_670 = rand(v);
%
for j=1:v;
    if j<=v;
        %% read the data for each model
        
        [t1] = readsac(['seis/s20/s20' num2str(j) '.sac']);
        
        %% compute the epicentral distance
        
        [s20,npts,delta3,b3,dist3,az3,baz3,gcarc3]=loadsac(['seis/s20/s20' num2str(j) '.sac']);
        [s20_topo,npts4,delta4,b4,dist4,az4,baz4,gcarc4]=loadsac(['seis/s20_topo/s20_topo' num2str(j) '.sac']);
        
        [delta,mlat,mlon] = get_delt_mid(['seis/s20_topo/s20_topo' num2str(j) '.sac']);
        
        mlt(j) = mlat;
        mln(j) = mlon;
        deltas(j) = delta;
        if mln(j) > 360 ;
            mln(j) = mln(j) -360;
        end
        
        mln = transpose(mln);
        mlt = transpose(mlt);
        %% test results with filtered data, if you use that change the name of the input data in windowing section
        %         %
        wn2 = 0.05; % low end of bandpass filter (Hz)/100s (Gu&Dziewonski2002-they use half-welch though)
        wn1 = 0.005; % high end of bandpass filter (Hz)/  40s
        wn = [wn1,wn2];
        wnlp = 0.025;% if you use the lowpass filter
        
        [z,p] = butter(2,wn,'bandpass');
        %                         filtp = filter(z,p,prem);
        %                         premF = transpose(filtp);
        %
        %                         filtpt = filter(z,p,prem_topo);
        %                         prem_topoF = transpose(filtpt);
        
        filts20 = filter(z,p,s20);
        s20F = transpose(filts20);
        
        filts20t = filter(z,p,s20_topo);
        s20_topoF = transpose(filts20t);
        
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
            
            dtt1 = abs(t1(2)-t1(1));
            t0 = t1(1);              % beginning time
            %             len = length(t);        % original time array
            len1 = length(t1);
            Fn = 1/(2*dtt1);          % Nyquist frequency
            
            ts = times_s400s(j)-50; % starting time of window - earliest S670S (i.e for 110deg)
            te = times_s400s(j)+50; % ending time of window   - latest SS (i.e for 170deg)
            
            % find indices for windows
            
            is = (ts-t0)/dtt1+1;
            ie = (te-t0)/dtt1+1;
            i1 = floor(is);
            i2 = floor(ie);
            
            wlen1 = i2-i1;	% length of the window
            
            [seis_win_s20,win_seis_s20] = windoweb(s20F,len1,i1,wlen1,10);
            [seis_win_s20t,win_seis_s20t] = windoweb(s20_topoF,len1,i1,wlen1,10);
            
            %% compute the distance at each midpoint
            % use then these deltas to interpolate to the taup curves and thus to find
            % the arrival time at that point
            
            % mid_arr_times = interp1(X1,Y1,mid_dist,'spline');
            x3 = seis_win_s20;
            x4 = seis_win_s20t;
            
            xl3 = length(x3);
            xl4 = length(x4);
            
            % correlation coefficients for all 4 combinations
            
            
            r_s20s20t = corrcoef(x4,x3);
            R_s20s20t(j) = r_s20s20t(1,2);
            R_s20s20t = transpose(R_s20s20t);
            
            % gives positions of signals which correlate well (>=0.5)
            suff = find(abs(R_s20s20t)>=0.5);
            
            %% Perform cross correlation under some particular requirements
            % NOTE1: WHen the signals are anti-correlated, i.e r<0 then if R>=0.5 take
            % estimate the time delay at the lag which corresponds to the
            % max correlation
            maxlg = 150;  % because max time delay cannot be larger than 42s
            
            if (R_s20s20t(j)) >= 0.5;
                xx3 = fft(x3);
                xx4 = fft(x4);
                [ccs20s20t,lags] = xcorr(x3,x4,maxlg);
                [ccmax,maxlag] = max(abs(ccs20s20t));
                
                delay2 = lags(maxlag);
                dt2(j) = delay2*dtt1;
                
            end
            
            
            %% same for S660S
            
            ts660 = times_s660s(j)-50; % starting time of window
            te660 = times_s660s(j)+50; % ending time of window
            
            % find indices for windows
            
            is660 = (ts660-t0)/dtt1+1;
            ie660 = (te660-t0)/dtt1+1;
            i1_660 = floor(is660);
            i2_660 = floor(ie660);
            
            wlen1_660 = i2_660-i1_660;	% length of the window
            
            [seis_win_s20660,win_seis_s20] = windoweb(s20F,len1,i1_660,wlen1_660,10);
            [seis_win_s20t660,win_seis_s20t] = windoweb(s20_topof,len1,i1_660,wlen1_660,10);
            
            
            x33 = transpose(seis_win_s20660);
            x44 = transpose(seis_win_s20t660);
            
            xl33 = length(x33);
            xl44 = length(x44);
            %
            %             fx1 = fft(x11);
            %             fx2 = fft(x22);
            %             fx3 = fft(x33);
            %             fx4 = fft(x44);
            %
            % correlation coefficients for all 4 combinations
            
            r_s20s20t660 = corrcoef(x44,x33);
            R_s20s20t660(j) = r_s20s20t660(1,2);
            
            
            % Quality check 1: keep only correlations higher or equal to 0.5
            
            suff22 = find(abs(R_s20s20t660)>=0.5);
            
            if (R_s20s20t660(j)) >= 0.5;
                xx33 = fft(x33);
                xx44 = fft(x44);
                [ccs20s20t660,lags660] = xcorr(x33,x44,maxlg);
                [ccmax2,maxlag2] = max(abs(ccs20s20t660));
                
                delay22 = lags660(maxlag2);
                dt22(j) = delay22*dtt1;
                
            end
            
        end
        %%
        % compute the predicted (i.e according to ray theory) travel time differences
        % this won't change for all regressions
        Vs_400 = 4.76; % km/s
        Vs_660 = 5.94; % km/s
        
        dh_400 = load('dh400');
        dh_660 = load('dh660');
        r4 = 5961; % radius to 400
        r6 = 5701; % radius to 670
        pe = 0.107 ;% p is the ray parameter equal to 11.8 s/deg;
        eta4 = r4./Vs_400;
        eta6 = r6./Vs_660;
        
        dt_pred_400(j) = (((-2)*dh_400(j))*sqrt(eta4.^2-pe.^2))./r4; % s
        dt_pred_660(j) = (((-2)*dh_660(j))*sqrt(eta6.^2-pe.^2))./r6; % s
        
        dt_pred_400t = transpose(dt_pred_400); % faster to work in columns
        dt_pred_660t = transpose(dt_pred_660);
    end
    waitbar(j/10000,h)
    
end


% 400

dt2 = transpose(dt2);
dt22 = transpose(dt22);

delays{:,:} = {dt2,dt22}

dt2_indx = transpose(find(dt2));
dt22_indx = transpose(find(dt22));

delays_indx{:,:} = {dt2_indx,dt22_indx};

%% perform linear regression to compare predicted and measured by cross correlation travel times
% note that this has to be adopted for all 4 combinations of models
% this is just a preface


cc2 = xcorr(dt2,dt_pred_400t(1:length(dt2)),'coeff');
ccf2 = corrcoef(dt2,dt_pred_400t(1:length(dt2)));

cc22 = xcorr(dt22,dt_pred_660t(1:length(dt22)),'COEFF');
ccf22 = corrcoef(dt22,dt_pred_660t(1:length(dt22)));

%
for k=1; % 4 combinations of models used now // need this to run through all dt vectors
    
    % 400 km discontinuity
    for x= 1:length(delays_indx{:,:}{:,k});
        
        l = delays_indx{:,:}{:,k}(x);
        
        x400(x) = (delays{:,:}{:,k}(l));
        
        y400(x) = dt_pred_400(l);
    end
    
    A=dt2;
    
    B=dt_pred_400t;
    
    figure
    scatter(B,A);
    hold on
    coef=polyfit(B,A,1);  % linear regression model
    m1=coef(1);
    n1=coef(2);
    Abest=m1*B+n1;
    mse=sum((A-Abest).^2)/length(A);    % mean square error
    
    plot(B,Abest,B,A,'o')
    title ('Linear regression at 400','fontsize',14)
    ylabel('Measured dt (s)','fontsize',12)
    xlabel('Calculated dt (s)','fontsize',12)
    xlim([-30 30])
    ylim([-30 30])
    grid on;
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1, ' ' ,'HorizontalAlignment','center','VerticalAlignment', 'top','fontsize',18);
    %legend(strcat('dt_pred=m2*dt_meas+n2', num2str(someValue)))
    
    orient landscape;
    print('-dpdf','LR400.pdf')
    
    
end

for kk= 1;
    
    % 660 km discontinuity
    for xx=1:length(delays_indx{:,1}{:,kk});
        
        ll = delays_indx{:,:}{:,kk}(xx);
        x660(xx) = (delays{:,:}{:,kk}(ll));
        y660(xx) = dt_pred_660(ll);
    end
    
    C=dt22;
    D=dt_pred_660t;
    
    figure
    scatter(D,C)
    hold on
    coef=polyfit(D,C,1);  % linear regression model
    m11=coef(1);
    n11=coef(2);
    Cbest=m11*D+n11;
    mse1=sum((C-Cbest).^2)/length(C);    % mean square error
    
    plot(D,Cbest,D,C,'o')
    title ('Linear regression at 660','fontsize',14)
    ylabel('Measured dt (s)','fontsize',12)
    xlabel('Calculated dt (s)','fontsize',12)
    xlim([-30 30])
    ylim([-30 30])
    grid on;
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1, ' ' ,'HorizontalAlignment','center','VerticalAlignment', 'top','fontsize',18);
    %legend(strcat('dt_pred=m2*dt_meas+n2', num2str(someValue)))
    
    orient landscape;
    print('-dpdf','LR660.pdf')
    
    
end

%%  Obtain the topography from the measured by cross correlation times

% 400 km discontinuity

% dh_meas400 = -2*dt2





















close(h)
toc
