function out=read_MULTI_code_release(name)

shift=[0.025 0.05 0.1 0.2 0.3 0.5 0.7];

% dt=0.005;
% Tpre=0.1;

%% read files

% specify the folder of the simulation results
direc='code_release/output/';

MFrates=read_check_file(direc,name,'/err_final.dat');

PC(1).delay=read_check_file(direc,name,'/PC_final_025.dat');
PC(2).delay=read_check_file(direc,name,'/PC_final_050.dat');
PC(3).delay=read_check_file(direc,name,'/PC_final_100.dat');
PC(4).delay=read_check_file(direc,name,'/PC_final_200.dat');
PC(5).delay=read_check_file(direc,name,'/PC_final_300.dat');
PC(6).delay=read_check_file(direc,name,'/PC_final_500.dat');
PC(7).delay=read_check_file(direc,name,'/PC_final_700.dat');

PCtarget=read_check_file(direc,name,'/htarget_BINS.dat');
J=read_check_file(direc,name,'/J.dat');

time=read_check_file(direc,name,'/time.dat');
timePC=read_check_file(direc,name,'/time_learn.dat');
timeGCtrans=read_check_file(direc,name,'/timeGCtrans.dat');

MF=read_check_file(direc,name,'/MF.dat');
MLI=read_check_file(direc,name,'/MLI.dat');
GoC=read_check_file(direc,name,'/GoC.dat');

GC=read_check_file(direc,name,'/GC.dat');
GCtrans=read_check_file(direc,name,'/GCtrans.dat');

CF=read_check_file(direc,name,'/CF.dat');

STP=read_check_file(direc,name,'/MFGC_STP.dat');
if(isempty(STP))
    timeSTP=[];
else
    timeSTP=STP(:,1);
    STP=STP(:,2:end);
end
STPpara=read_check_file(direc,name,'/STPpara.dat');

parameters=read_check_file(direc,name,'/learnparam.dat');
flags=read_check_file(direc,name,'/flags.dat');
stats=read_check_file(direc,name,'/stats.dat');

%% sizes
nt=size(PC(1).delay,2);
ndelays=length(PC);

n=unique(parameters(:,8),'stable');
ntrials = size(MFrates, 1);
llc  = ntrials - mod(ntrials, n);

%% average over different realisations
for k=1:ndelays
    PC_mean(k).delay=[];
    for i=1:nt
        y  = reshape(PC(k).delay(1:llc,i), n, []);
        %     Result = transpose(nansum(y, 1) / n);
        Result = transpose(sum(y, 1) / n);
        PC_mean(k).delay=[PC_mean(k).delay Result]; %#ok<*AGROW>
    end
end

%% get statistics from final PC responses
for j=1:ndelays
    temp=eye_MIN_TE_FWHM(timePC,PC(j).delay,shift(j));
    Merr(j).mean=[];
    Merr(j).std=[];
    for i=1:4
        y  = reshape(temp(1:llc,i), n, []);
        Result = transpose(nansum(y, 1) / n);
        Merr(j).mean=[Merr(j).mean Result];
        
        Result= transpose(nanstd(y, 0, 1));
        Merr(j).std=[Merr(j).std Result];
    end
end

Nd=length(Merr);
temp=[];
Merr_single=[];
for j=1:Nd
    temp=[temp sum(abs(Merr(j).mean(:,1:3)),2)];
end
for j=1:Nd
    Merr_single(:,j)=sum(abs(Merr(j).mean(:,1:3)),2);
end
Merr_tot=sum(temp(:,1:7),2);
Merr_tot_late=sum(temp(:,5:7),2);
Merr_tot_early=sum(temp(:,1:4),2);

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
out.J=J';
out.Merr=Merr;
out.Merr_single=Merr_single;
out.Merr_tot=Merr_tot;
out.Merr_tot_early=Merr_tot_early;
out.Merr_tot_late=Merr_tot_late;

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

end
