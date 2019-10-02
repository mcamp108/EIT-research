function out= testing(seq, imdl)

% Description:
%   
%
% Parameters:
%
%
% Returns:
%
%
% Author:
%   Mark Campbell
%   Carleton University
%   markacampbell@cmail.carleton.ca

% 23.Sep.2019

% once bad elecs are identified
badelecs = [1];
kk=meas_icov_rm_elecs(imdl,badelecs);
ee = find(diag(kk)~=1);
mm = find(msel);
ee = mm(ee);


% fmdl= ng_mk_cyl_models(2,[32,1],[0.05]);
[stim,msel] = mk_stim_patterns(32,1,[0,5],[0,5],{'no_meas_current_next2','no_rotate_meas'},1);
% fmdl.meas_select = msel;
% fmdl = mdl_normalize(fmdl,0);
% fmdl.stimulation = stim;
% img = mk_image(fmdl,1);
% opt= struct('imgsz', [32 32], 'noise_figure', .5); imdl=mk_GREIT_model(img, 0.25, [], opt);
% 
% clf
% pth='/svn/lab_data/2017/UKE-Hamburg/EIT Brain/'; for fnum = 1:1; switch fnum;
%    case 1; fn='EIT_nativ_2_zvk_10.2.eit';
%       hidx=5032; iidx=[6000:10:6300];
%       cut = [310,6700];
%    otherwise; error 'huh?'
% end

vv= seq.eit.data;
FR = seq.eit.fs;

ll = size(vv,2);

subplot(421);
    plot(1e3*vv(:,1:100), 'k+'); hold on;
    plot(1e3*vv(msel,1:100), 'b+');
    badelecs = [1];
    kk=meas_icov_rm_elecs(imdl,badelecs);
    ee = find(diag(kk)~=1);
    mm = find(msel);
    ee = mm(ee);
    plot(1e3*vv(ee,1:100), 'r+');
    hold off; box off
    title 'IQ plot';

subplot(422);
    vk = 1e3*real(vv(:,1:100));
    plot(vk,'k');      hold on;
    vk(~msel,:) = NaN;
    plot(vk,'b');
    vn = NaN*vk; vn(ee,:) = vk(ee,:);
    plot(vn,'ro');
    hold off;
    box off; xlim([0,400]);
    title 'U shapes';

subplot(423);
    ei = abs(aux.elec_impedance);
    aei = median(ei,2);
    bar(1:32,aei);
    hold on;
    eb= errorbar(1:32,aei,max(ei,[],2)-aei,min(ei,[],2)-aei);
    set(eb,'Color',[0,0,0],'LineStyle','none');
    hold off;
    box off; xlim([1,32]);
    title 'Elec Z';

subplot(424);
    ft = fft(detrend(vv.').',[],2);
    fax = linspace(0,FR,ll+1); fax(end)=[];
    semilogy(fax,mean(abs(ft)));
    ylim([1e-4,1]);
    box off; xlim([0,12]);
    title 'Spectrum (Hz)'

vf = freq_filt(real(vv),@(f) f<6, FR, 2);

FR = 1e6/median(diff(aux.t_rel));
tt= (1:ll)/FR;
ss= sum(vf); mm = [min(ss);max(ss)];
[~,hidx]= min(ss);
subplot(413);
    plot([1;1]*tt([hidx,iidx]), ...
        mm*ones(1,length(iidx)+1),'Color',0.5*[1,1,1]);
hold on;
    plot(tt,sum(vf),'Linewidth',2,'Color',[0,0,0.5]);
hold off;
    xlim([0,max(tt)]);
    box off; title 'Z vs t (s)';

subplot(414);
imgr = inv_solve(imdl,vf(:,hidx),vf(:,iidx));
imgr.show_slices.img_cols = 11;
show_slices(imgr);
title(fn)
end


end % end function

