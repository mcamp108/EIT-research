dd= 

find(any((dd-0.0014057)<1e-7)')


find(any((dd'-0.0014057)<1e-7))


find(any((dd'-0.0014057)<1e-8))


find(any((dd'-0.0014057)<1e-10))

find(any(abs((dd'-0.0014057)<1e-10)))

find(any(abs((dd'-0.0014057))<1e-10))

find(any(abs((dd'-0.0014057))<1e-8))

kk=meas_icov_rm_elecs(imdl, [28,1]); ee = find(diag(kk)~=1);plot(1:4638,dd','k', 1:4638,dd(ee,:)','r')
kk=meas_icov_rm_elecs(imdl, [28,32,1]); ee = find(diag(kk)~=1);plot(1:4638,dd','k', 1:4638,dd(ee,:)','r')
plot(1:4638,dd(ee,:)','Color',[1,0.7,0.7]); hold on; plot(1:4638,dd(notee,:)','k'); hold off
size(ee)

notee=1:544; notee(ee)=[];
plot(1:4638,dd(ee,:)','Color',[1,0.7,0.7]); hold on; plot(1:4638,dd(notee,:)','k'); hold off
find(any(abs((dd'-0.0013386))<1e-8))

find(any(abs((dd'-0.0014682))<1e-8))

find(any(abs((dd'-0.0014577))<1e-8))

format compact

find(any(abs((dd'-0.0014577))<1e-8))

find(any(abs((dd'-0.0013386))<1e-8))

find(any(abs((dd'-0.0014577))<1e-8))

find(any(abs((dd'-0.0014681))<1e-8))

find(any(abs((dd'-0.0013463))<1e-8))

find(any(abs((dd'-0.0013542))<1e-8))

find(any(abs((dd'-0.0013491))<1e-8))


for be=1:32; kk=meas_icov_rm_elecs(imdl, be); disp([be, kk(263,263)]); end
  
for be=1:32; kk=meas_icov_rm_elecs(imdl, be); disp(full([be, kk(263,263)])); end
     

find(any(abs((dd'-0.0014577))<1e-8))

kk=meas_icov_rm_elecs(imdl, [28,32,1]); ee = find(diag(kk)~=1);plot(1:4638,dd','k', 1:4638,dd(ee,:)','r')
notee=1:544; notee(ee)=[];
kk=meas_icov_rm_elecs(imdl, [28,32,1]); ee = find(diag(kk)~=1);plot(1:4638,dd','k', 1:4638,dd(ee,:)','r');notee=1:544; notee(ee)=[];
plot(1:4638,dd(ee,:)','Color',[1,0.7,0.7]); hold on; plot(1:4638,dd(notee,:)','k'); hold off
kk=meas_icov_rm_elecs(imdl, [28,32,1,8]); ee = find(diag(kk)~=1);plot(1:4638,dd','k', 1:4638,dd(ee,:)','r');notee=1:544; notee(ee)=[];
plot(1:4638,dd(ee,:)','Color',[1,0.7,0.7]); hold on; plot(1:4638,dd(notee,:)','k'); hold off
kk=meas_icov_rm_elecs(imdl, [28,32,1,16]); ee = find(diag(kk)~=1);plot(1:4638,dd','k', 1:4638,dd(ee,:)','r');notee=1:544; notee(ee)=[];
plot(1:4638,dd(ee,:)','Color',[1,0.7,0.7]); hold on; plot(1:4638,dd(notee,:)','k'); hold off
kk=meas_icov_rm_elecs(imdl, [28,32,1,21]); ee = find(diag(kk)~=1);plot(1:4638,dd','k', 1:4638,dd(ee,:)','r');notee=1:544; notee(ee)=[];
plot(1:4638,dd(ee,:)','Color',[1,0.7,0.7]); hold on; plot(1:4638,dd(notee,:)','k'); hold off
kk=meas_icov_rm_elecs(imdl, [28,32,1,13]); ee = find(diag(kk)~=1);plot(1:4638,dd','k', 1:4638,dd(ee,:)','r');notee=1:544; notee(ee)=[];
plot(1:4638,dd(ee,:)','Color',[1,0.7,0.7]); hold on; plot(1:4638,dd(notee,:)','k'); hold off
plot(sum(dd(notee,:)))
plot(sum(dd(notee,:)))
                      ?
plot(sum(dd(notee,:)))
plot(1:4638,dd(ee,:)','Color',[1,0.7,0.7]); hold on; plot(1:4638,dd(notee,:)','k'); hold off
plot(1:4638,dd(ee,:)','Color',[1,0.7,0.7]); hold on; plot(1:4638,dd(notee,:)','k'); hold off
plot(sum(dd(notee,:)))

for be=1:32; kk=meas_icov_rm_elecs(imdl, be); ee = find(diag(kk)~=1); plot(dd(ee,'k')'); title(sprintf('bad=%d',be)); pause; end

% plot each electrode and look for worst ones.
for be=1:32; kk=meas_icov_rm_elecs(imdl, be); ee = find(diag(kk)~=1); plot(dd(ee,:)','k'); title(sprintf('bad=%d',be)); pause; end
plot(sum(dd(notee,:)))
