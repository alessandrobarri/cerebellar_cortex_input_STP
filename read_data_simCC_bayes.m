function out=read_MULTI_bayes_code_release(name)

Iini=[0.025, 0.05, 0.1, 0.2, 0.3];
Iend=[0.15, 0.2, 0.3, 0.4, 0.5];
Rini=[0.015, 0.025, 0.05, 0.1, 0.2];
Rend=[0.2, 0.3, 0.4, 0.5, 0.6];

%% read files
% specify the folder of the simulation results
directory='code_release/output/';
MFrates=read_check_file(directory,name,'/err_final.dat');

path_PC_025150=[directory name '/PC_bayes_025150.dat'];
path_PC_050200=[directory name '/PC_bayes_050200.dat'];
path_PC_100300=[directory name '/PC_bayes_100300.dat'];
path_PC_200400=[directory name '/PC_bayes_200400.dat'];
path_PC_300500=[directory name '/PC_bayes_300500.dat'];
PC(1).delay=dlmread(path_PC_025150);
PC(2).delay=dlmread(path_PC_050200);
PC(3).delay=dlmread(path_PC_100300);
PC(4).delay=dlmread(path_PC_200400);
PC(5).delay=dlmread(path_PC_300500);

PCtarget=read_check_file(directory,name,'/htarget_BINS.dat');
J=read_check_file(directory,name,'/J.dat');

time=read_check_file(directory,name,'/time.dat');
timePC=read_check_file(directory,name,'/time_learn.dat');
timeGCtrans=read_check_file(directory,name,'/timeGCtrans.dat');

MF=read_check_file(directory,name,'/MF.dat');
MLI=read_check_file(directory,name,'/MLI.dat');
GoC=read_check_file(directory,name,'/GoC.dat');

GC=read_check_file(directory,name,'/GC.dat');
GCtrans=read_check_file(directory,name,'/GCtrans.dat');

CF=read_check_file(directory,name,'/CF.dat');

STP=read_check_file(directory,name,'/MFGC_STP.dat');
timeSTP=STP(:,1);
STP=STP(:,2:end);
STPpara=read_check_file(directory,name,'/STPpara.dat');

parameters=read_check_file(directory,name,'/learnparam.dat');
flags=read_check_file(directory,name,'/flags.dat');
stats=read_check_file(directory,name,'/stats.dat');

%% sizes
nt=size(PC(1).delay,2);     % number of time steps
ndelays=length(PC);

n=unique(parameters(:,8),'stable');
ll = size(MFrates, 1);         % number of trials
llc  = ll - mod(ll, n);

%% compute DN responses
dt=timePC(2)-timePC(1);
pidx=timePC>=0;
for j=1:ndelays
    for k=1:ll
        PCmean=mean(PC(j).delay(k,:));
        DN(j).delay(k,pidx)=dt*cumsum(PCmean-PC(j).delay(k,pidx));
    end
end

%% compute TRACE estimates and bayesian theory for naive w=0.1
for j=1:ndelays
    idxt=find(timePC>Rini(j) & timePC<Rend(j));
    [te,ts,tm,te_var]=BLS_timing(Iini(j),Iend(j),Rini(j),Rend(j),'dt',dt);
    TRACE(j).te=te;
    TRACE(j).ts=ts;
    TRACE(j).tm=tm;
    TRACE(j).te_var=te_var;
    TRACE(j).time=timePC(idxt);
    for k=1:ll
        DNtemp=DN(j).delay(k,:);
        trace=(Iend(j)-Iini(j))*(DNtemp-min(DNtemp))/(max(DNtemp)-min(DNtemp)) +Iini(j);
        TRACE(j).TR(k,:)=trace(idxt);
    end
end

%% average over different realisations
for k=1:ndelays
    PC_mean(k).delay=[];
    for i=1:nt
        y  = reshape(PC(k).delay(1:llc,i), n, []);
        Result = transpose(sum(y, 1) / n);
        PC_mean(k).delay=[PC_mean(k).delay Result]; %#ok<*AGROW>
    end
end

for k=1:ndelays
    DN_mean(k).delay=[];
    for i=1:nt
        y  = reshape(DN(k).delay(1:llc,i), n, []);
        Result = transpose(sum(y, 1) / n);
        DN_mean(k).delay=[DN_mean(k).delay Result]; %#ok<*AGROW>
    end
end

for k=1:ndelays
    TRACE_mean(k).TR=[];
    for i=1:length(TRACE(k).TR)
        y  = reshape(TRACE(k).TR(1:llc,i), n, []);
        Result = transpose(sum(y, 1) / n);
        TRACE_mean(k).TR=[TRACE_mean(k).TR Result]; %#ok<*AGROW>
    end
end

%% output
out.MFrates=MFrates;
if(~isempty(GC))
    out.basis.time=time;
    out.basis.timeSTP=timeSTP;
    out.basis.MF=MF;
    out.basis.GC=GC;
    out.basis.GCtrans=GCtrans';
    out.basis.timeGCtrans=timeGCtrans;
    out.basis.MLI=MLI;
    out.basis.GoC=GoC;
    out.basis.CF=CF;
    out.basis.STP=STP;
    out.basis.STPpara=STPpara;
end
out.time=timePC;
out.PCtarget=PCtarget;
out.PC=PC;
out.PC_mean=PC_mean;
out.DN=DN;
out.DN_mean=DN_mean;
out.TRACE=TRACE;
out.TRACE_mean=TRACE_mean;
out.J=J';

%% sort through simulation parameters
out.parameters.Ntrials=unique(parameters(:,1),'stable');
out.parameters.Nreal=unique(parameters(:,8),'stable');
out.parameters.Npara=unique(parameters(:,9),'stable');
out.parameters.Npara2=unique(parameters(:,10),'stable');

out.parameters.CFsp=unique(parameters(:,2),'stable');
out.parameters.CFscale=unique(parameters(:,3),'stable');
out.parameters.errWeight=unique(parameters(:,5),'stable');
out.parameters.alpha_learn=unique(parameters(:,6),'stable');
out.parameters.L2weight=unique(parameters(:,32),'stable');

out.parameters.nfrac_GC=unique(parameters(:,15),'stable');
out.parameters.GC_avrg_target=unique(parameters(:,16),'stable');
out.parameters.GC_normfac=unique(parameters(:,17),'stable');
out.parameters.GC_thr_reduc_fac=unique(parameters(:,18),'stable');

out.parameters.JI=unique(parameters(:,7),'stable');
out.parameters.avrgJIE=unique(parameters(:,12),'stable');
out.parameters.avrgJEI=unique(parameters(:,13),'stable');
out.parameters.gain_global=unique(parameters(:,14),'stable');

out.parameters.MF_U_corr=unique(parameters(:,19),'stable');
out.parameters.Ncl=unique(parameters(:,20),'stable');
out.parameters.pact=unique(parameters(1,21:25),'stable','rows');
out.parameters.grouplims=unique(parameters(:,26:31),'stable','rows');

out.parameters.RNseed=unique(parameters(:,11),'stable');

%% sort through stats
out.stats.hGC_avrg=unique(stats(:,1),'stable');
out.stats.hGC_var=unique(stats(:,2),'stable');
out.stats.Thresh_avrg=unique(stats(:,3),'stable');
out.stats.gain_avrg=unique(stats(:,4),'stable');
out.stats.CL_avrg=unique(stats(:,5),'stable');

%% sort through flags
out.flags.MF_STP=unique(flags(:,1),'stable');
out.flags.MF_patterns=unique(flags(:,3),'stable');
out.flags.MF_groups=unique(flags(:,6),'stable');
out.flags.GC_threshold=unique(flags(:,2),'stable');
out.flags.GC_gain=unique(flags(:,4),'stable');
out.flags.drv_supp_struct=unique(flags(:,7),'stable');
out.flags.STP_transients=unique(flags(:,8),'stable');
out.flags.U_distrib=unique(flags(:,9),'stable');
out.flags.mode=unique(flags(:,10),'stable');
out.flags.GoC=unique(flags(:,11),'stable');
out.flags.Tsignal=unique(flags(:,12),'stable');
out.flags.syn_cut=unique(flags(:,13),'stable');
out.flags.GC_cut=unique(flags(:,14),'stable');

%% identify flags
if(out.flags.MF_STP==0)
    out.flags.MF_STP='off';
elseif(out.flags.MF_STP==1)
    out.flags.MF_STP='on';
elseif(out.flags.MF_STP==2)
    out.flags.MF_STP='all STP identical';
end

if(out.flags.GC_threshold==0)
    out.flags.GC_threshold='fixed';
elseif(out.flags.GC_threshold==1)
    out.flags.GC_threshold='adjusted individually';
elseif(out.flags.GC_threshold==2)
    out.flags.GC_threshold='weights renormalised';
elseif(out.flags.GC_threshold==4)
    out.flags.GC_threshold='adjusted globaly';
end

if(out.flags.GC_gain==0)
    out.flags.GC_gain='adjusted globaly';
elseif(out.flags.GC_gain==1)
    out.flags.GC_gain='adjusted individually (based on GC average)';
elseif(out.flags.GC_gain==2)
    out.flags.GC_gain='adjusted individually (based on GC average when active)';
elseif(out.flags.GC_gain==3)
    out.flags.GC_gain='fixed';
end

if(out.flags.MF_patterns==0)
    out.flags.MF_patterns='thresholded Gaussian';
elseif(out.flags.MF_patterns==1)
    out.flags.MF_patterns='lognormal';
elseif(out.flags.MF_patterns==2)
    out.flags.MF_patterns='exponential';
elseif(out.flags.MF_patterns==3)
    out.flags.MF_patterns='gamma';
elseif(out.flags.MF_patterns==4)
    out.flags.MF_patterns='multiple uniform distributions';
elseif(out.flags.MF_patterns==5)
    out.flags.MF_patterns='single, seamsless uniform distribution';
elseif(out.flags.MF_patterns==6)
    out.flags.MF_patterns='truncated Gaussian';
end

if(out.flags.MF_groups==2)
    out.flags.MF_groups='5 MF types and 5 synapse types with corr=1 (native)';
elseif(out.flags.MF_groups==5)
    out.flags.MF_groups='2 MF types and 2 synapse types with corr=MFUcorr (simplified SD model)';
elseif(out.flags.MF_groups==55)
    out.flags.MF_groups='MF distrb. with 5 peaks and Ncl MF/SYN types with corr=MFUcorr (fig5)';
end

if(out.flags.drv_supp_struct==0)
    out.flags.drv_supp_struct='random MF assignment to GC';
elseif(out.flags.drv_supp_struct==1)
    out.flags.drv_supp_struct='at least one driver MF per GC';
elseif(out.flags.drv_supp_struct==2)
    out.flags.drv_supp_struct='for SD model: 2 drivers + 2 supporters';
end

if(out.flags.STP_transients==0)
    out.flags.STP_transients='off';
elseif(out.flags.STP_transients==1)
    out.flags.STP_transients='on';
end

if(out.flags.U_distrib==0)
    out.flags.U_distrib='all U within group identical';
elseif(out.flags.U_distrib==2)
    out.flags.U_distrib='uniform distribution';
end

if(out.flags.mode==6)
    out.flags.mode='scan over correlation between MF rate & rel. prob.';
elseif(out.flags.mode==7)
    out.flags.mode='scan over cI';
elseif(out.flags.mode==77)
    out.flags.mode='bayesian learning';
elseif(out.flags.mode==88)
    out.flags.mode='scan over MF rates';
elseif(out.flags.mode==99)
    out.flags.mode='sample eye-lid simulations';
end

if(out.flags.GoC==0)
    out.flags.GoC='off';
elseif(out.flags.GoC==1)
    out.flags.GoC='on';
end

if(out.flags.Tsignal==0)
    out.flags.Teaching_signal='Gauss';
elseif(out.flags.Tsignal==1)
    out.flags.Teaching_signal='delta';
end

if(out.flags.syn_cut==0)
    out.flags.syn_cut='no cut';
elseif(out.flags.syn_cut==2)
    out.flags.syn_cut='cut supporters';
elseif(out.flags.syn_cut==3)
    out.flags.syn_cut='cut drivers';
elseif(out.flags.syn_cut==4)
    out.flags.syn_cut='cut slow pools';
elseif(out.flags.syn_cut==5)
    out.flags.syn_cut='cut fast pools';
end

if(out.flags.GC_cut==0)
    out.flags.GC_cut='no cut';
elseif(out.flags.GC_cut==1)
    out.flags.GC_cut='cut above';
elseif(out.flags.GC_cut==2)
    out.flags.GC_cut='cut below';
end

out.Teaching_signal_type=out.flags.Teaching_signal;
out.mode=out.flags.mode;

end

%% auxiliary functions
function [out,stemp]=read_check_file(direc,name,fname)
out=[];
path=[direc name fname];
if isfile(path)
    stemp=dir(path);
    if(stemp.bytes~=0)
        out=dlmread(path);
    end
end
end

function [te,ts,tm,te_var]=BLS_timing(tmin,tmax,tpmin,tpmax,varargin)
% tmin, tmax are the bounds of the uniform interval of delays in RSG task
% tpmin, tpmax are the bounds of interval for which the BLS estimator is to
% be calculated;

w=0.1;

dt=0.002;
timeflag=0;
if(nargin>4)
    v1 = 1;
    while v1<=length(varargin)
        if strcmp(varargin{v1}, 'dt')
            dt=varargin{v1+1};
            v1=v1+2;
        elseif strcmp( varargin{v1}, 'time' )
            time = varargin{v1+1};
            timeflag=1;
            v1=v1+2;
        elseif strcmp( varargin{v1}, 'w' )
            w = varargin{v1+1};
            v1=v1+2;
        else
            v1=v1+1;
        end
    end
end

if(timeflag)
    tm=time;
    ts=tm;
else
    tm=linspace(tpmin*0.9,tpmax*1.1,(tpmax*1.1-tpmin*0.9)/dt);
    ts=linspace(tpmin,tpmax,(tpmax-tpmin)/dt);
end

prior=@(ts,tmin,tmax)unifpdf(ts,tmin,tmax);
likelihood=@(tm,ts,w) exp(-(ts-tm).^2 ./(2*(w.*ts).^2) )./(sqrt(2*pi)*w*ts);

%% te vs tm
integrand_ts=@(ts,tm,w,tmin,tmax) ts.*prior(ts,tmin,tmax).*likelihood(tm,ts,w);
integrand_ts2=@(ts,tm,w,tmin,tmax) ts.*ts.*prior(ts,tmin,tmax).*likelihood(tm,ts,w);
integrand_norm=@(ts,tm,w,tmin,tmax) prior(ts,tmin,tmax).*likelihood(tm,ts,w);

te=zeros(1,length(tm));
te_var=zeros(1,length(tm));
for i=1:length(tm)   
    normfac=integral(@(x)integrand_norm(x,tm(i),w,tmin,tmax),tmin,tmax);
    te(i)=integral(@(x)integrand_ts(x,tm(i),w,tmin,tmax),tmin,tmax)/normfac;
    te_var(i)=integral(@(x)integrand_ts2(x,tm(i),w,tmin,tmax),tmin,tmax)/normfac - te(i)*te(i);
end

end