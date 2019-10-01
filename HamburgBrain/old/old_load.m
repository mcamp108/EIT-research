function D= load_HamburgBrain_data(pig)

% Each sequence in the data struct is an experimental trial. Each sequence
% has the following fields:

% name: Formatted file name used as title of figures for this sequence
% pig: Pig name
% imgr: Reconstructed images using eit.fdata and imdl for this pig
% eit: EIT part of data struct
%     data: EIT data output from eidors_readdata()
%     fs: EIT data framerate
%     sync1: First syncronization spike in EIT data
%     sync2: Second syncronization spike in EIT data
%     apn: Using the LabChart file, the time relative to the beginning of the first EIT synchronization spike that the apnoea flush occured in the arterial pressure data
%     inj: Same as above, but for the injection flush
%     vnt: Same as above, but for the end of apnoea flush
%     fdata: 9 Hz lowpass-filtered EIT data
% perf: Arterial pressure part of data struct
%     data: Perfusion data
%     tickrate: Arterial pressure framerate
%     apn: Apnoea flush spike in arterial pressure data
%     inj: Injection flush spike in arterial pressure data
%     vnt: End of apnoea flush spike in arterial pressure data
% fdata now calculated in align sequence to trim filter edge effect.
% 
% To change the lowpass filter stopband, go to lowpass_iir function.
% 
% Author:
%   Mark Campbell
%   Carleton University
%   markacampbell@cmail.carleton.ca
eitfs= 47.68;
maxsz= 0.2; maxh= 2; imgsize= [64 64];
[fmdl, imdl]= mk_pighead_fmdl(maxsz, maxh, imgsize);
% lpass= 9;

% 8.2
if pig== "8.2"
    files= {'EIT_8.2_nativ_1.eit',                  'EIT_8.2_nativ_2.eit',... 
            'EIT_8.2_30_Min_nach_Embolisation.eit', 'EIT_8.2_rechts_embolisiert.eit'};
    sync1= [681, 98, 941, 1027];
    sync2= [6370, 4990, 4633, 5264];
    e_apn= [42.5, 14.5, -9.6, 15.7];
    e_inj= [50.8, 31.1, 12.7, 26.1];
    e_vnt= [98.2, 83.9, 58, 78.8];
    p_apn= [53517, 14282, 7622, 30161];
    p_inj= [61910, 30896, 30011, 40589];
    p_vnt= [109097, 83697, 75311, 93234];
    
function D= load_data(files, sync1, sync2, e_apn, e_inj, e_vnt, p_apn, p_inj, p_vnt)

end % end function
A.seq1.eit.data= eidors_readdata('EIT_8.2_nativ_1.eit');
A.seq1.name= char(remove_underscores('EIT_8.2_nativ_1.eit'));
A.seq1.pig= "8.2";
A.seq1.eit.fs= eitfs;
A.seq1.eit.sync1= 681;
A.seq1.eit.sync2= 6370;
A.seq1.eit.apn= A.seq1.eit.sync1 + round(42.5* eitfs);
A.seq1.eit.inj= A.seq1.eit.sync1 + round(50.8* eitfs);
A.seq1.eit.vnt= A.seq1.eit.sync1 + round(98.2* eitfs);
A.seq1.perf= load('EIT_8.2_nativ_1.mat');
A.seq1.perf.apn= 53517;
A.seq1.perf.inj= 61910;
A.seq1.perf.vnt= 109097;
A.seq1= allign_eit_and_perf(A.seq1);
% A.seq1.eit.fdata= lowpass_iir(A.seq1.eit.data, lpass, eitfs);
A.seq1.imgr= get_imgr(A.seq1, imdl);

A.seq2.eit.data= eidors_readdata('EIT_8.2_nativ_2.eit');
A.seq2.name= char(remove_underscores('EIT_8.2_nativ_2.eit'));
A.seq2.pig= "8.2";
A.seq2.eit.fs= eitfs;
A.seq2.eit.sync1= 98; % num frames in time from second sync to frame 1 in data (-1:44.6) - same in matlab file (-1:42.6)
A.seq2.eit.sync2= 4990;
A.seq2.eit.apn= A.seq2.eit.sync1 + round(14.5* eitfs);
A.seq2.eit.inj= A.seq2.eit.sync1 + round(31.1* eitfs);
A.seq2.eit.vnt= A.seq2.eit.sync1 + round(83.9* eitfs);
A.seq2.perf= load('EIT_8.2_nativ_2.mat');
A.seq2.perf.apn= 14282;
A.seq2.perf.inj= 30896;
A.seq2.perf.vnt= 83697;
A.seq2= allign_eit_and_perf(A.seq2);
% A.seq2.eit.fdata= lowpass_iir(A.seq2.eit.data, lpass, eitfs);
A.seq2.imgr= get_imgr(A.seq2, imdl);

A.seq3.eit.data= eidors_readdata('EIT_8.2_30_Min_nach_Embolisation.eit');
A.seq3.name= char(remove_underscores('EIT_8.2_30_Min_nach_Embolisation.eit'));
A.seq3.pig= "8.2";
A.seq3.eit.fs= eitfs;
A.seq3.eit.sync1= 941;
A.seq3.eit.sync2= 4633;
A.seq3.eit.apn= A.seq3.eit.sync1 + round(-9.6* eitfs);
A.seq3.eit.inj= A.seq3.eit.sync1 + round(12.7* eitfs);
A.seq3.eit.vnt= A.seq3.eit.sync1 + round(58* eitfs);
A.seq3.perf= load('EIT_8.2_30_Min_nach_Embolisation.mat');
A.seq3.perf.apn= 7622;
A.seq3.perf.inj= 30011;
A.seq3.perf.vnt= 75311;
A.seq3= allign_eit_and_perf(A.seq3);
% A.seq3.eit.fdata= lowpass_iir(A.seq3.eit.data, lpass, eitfs);
A.seq3.imgr= get_imgr(A.seq3, imdl);

A.seq4.eit.data= eidors_readdata('EIT_8.2_rechts_embolisiert.eit');
A.seq4.name= char(remove_underscores('EIT_8.2_rechts_embolisiert.eit'));
A.seq4.pig= "8.2";
A.seq4.eit.fs= eitfs;
A.seq4.eit.sync1= 1027;
A.seq4.eit.sync2= 5264;
A.seq4.eit.apn= A.seq4.eit.sync1 + round(15.7* eitfs);
A.seq4.eit.inj= A.seq4.eit.sync1 + round(26.1* eitfs);
A.seq4.eit.vnt= A.seq4.eit.sync1 + round(78.8* eitfs);
A.seq4.perf= load('EIT_8.2_rechts_embolisiert.mat');
A.seq4.perf.apn= 30161;
A.seq4.perf.inj= 40589;
A.seq4.perf.vnt= 93234;
A.seq4= allign_eit_and_perf(A.seq4);
% A.seq4.eit.fdata= lowpass_iir(A.seq4.eit.data, lpass, eitfs);
A.seq4.imgr= get_imgr(A.seq4, imdl);
D= A;

% 9.2
elseif pig== "9.2"
B.seq1.eit.data= eidors_readdata('EIT_9.2_Nativ_1.eit');
B.seq1.name= char(remove_underscores('EIT_9.2_Nativ_1.eit'));
B.seq1.pig= "9.2";
B.seq1.eit.fs= eitfs;
B.seq1.eit.sync1= 246;
B.seq1.eit.sync2= 5239;
B.seq1.eit.apn= B.seq1.eit.sync1 + round(12.3* eitfs);
B.seq1.eit.inj= B.seq1.eit.sync1 + round(27* eitfs);
B.seq1.eit.vnt= B.seq1.eit.sync1 + round(82.1* eitfs);
B.seq1.perf= load('EIT_9.2_Nativ_1.mat');
B.seq1.perf.apn= 132012;
B.seq1.perf.inj= 146787;
B.seq1.perf.vnt= 201855;
B.seq1= allign_eit_and_perf(B.seq1);
% B.seq1.eit.fdata= lowpass_iir(B.seq1.eit.data, lpass, eitfs);
B.seq1.imgr= get_imgr(B.seq1, imdl);

B.seq2.eit.data= eidors_readdata('EIT_9.2_Nativ_2.eit');
B.seq2.name= char(remove_underscores('EIT_9.2_Nativ_2.eit'));
B.seq2.pig= "9.2";
B.seq2.eit.fs= eitfs;
B.seq2.eit.sync1= 291; % first syncronization spike in EIT data
B.seq2.eit.sync2= 5948; % second syncronization spike in EIT data
B.seq2.eit.apn= B.seq2.eit.sync1 + round(12.1* eitfs);
B.seq2.eit.inj= B.seq2.eit.sync1 + round(19.9* eitfs);
B.seq2.eit.vnt= B.seq2.eit.sync1 + round(98.7* eitfs);
B.seq2.perf= load('EIT_9.2_Nativ_2.mat');
B.seq2.perf.apn= 140143; % apnoea flush spike in arterial pressure data
B.seq2.perf.inj= 148520; % injection flush spike in arterial pressure data
B.seq2.perf.vnt= 226761; % end of apnoea flush spike in arterial pressure data
B.seq2= allign_eit_and_perf(B.seq2);
% B.seq2.eit.fdata= lowpass_iir(B.seq2.eit.data, lpass, eitfs);
B.seq2.imgr= get_imgr(B.seq2, imdl);

B.seq3.eit.data= eidors_readdata('EIT_nach_Embolisation_1_9.2.eit');
B.seq3.name= char(remove_underscores('EIT_nach_Embolisation_1_9.2.eit'));
B.seq3.pig= "9.2";
B.seq3.eit.fs= eitfs;
B.seq3.eit.sync1= 979;
B.seq3.eit.sync2= 5956;
B.seq3.eit.apn= B.seq3.eit.sync1 + round(14.7* eitfs);
B.seq3.eit.inj= B.seq3.eit.sync1 + round(20.7* eitfs);
B.seq3.eit.vnt= B.seq3.eit.sync1 + round(83.9* eitfs);
B.seq3.perf= load('EIT_nach_Embolisation_1_9.2.mat');
B.seq3.perf.apn= 173558;
B.seq3.perf.inj= 179538;
B.seq3.perf.vnt= 242739;
B.seq3= allign_eit_and_perf(B.seq3);
% B.seq3.eit.fdata= lowpass_iir(B.seq3.eit.data, lpass, eitfs);
B.seq3.imgr= get_imgr(B.seq3, imdl);

B.seq4.eit.data= eidors_readdata('EIT_nach_Perfusionsminderun_2_9.2.eit');
B.seq4.name= char(remove_underscores('EIT_nach_Perfusionsminderun_2_9.2.eit'));
B.seq4.pig= "9.2";
B.seq4.eit.fs= eitfs;
B.seq4.eit.sync1= 239;
B.seq4.eit.sync2= 6382;
B.seq4.eit.apn= B.seq4.eit.sync1 + round(33.7* eitfs);
B.seq4.eit.inj= B.seq4.eit.sync1 + round(41* eitfs);
B.seq4.eit.vnt= B.seq4.eit.sync1 + round(112.1* eitfs);
B.seq4.perf= load('EIT_nach_Perfusionsminderun_2_9.2.mat');
B.seq4.perf.apn= 236714;
B.seq4.perf.inj= 243968;
B.seq4.perf.vnt= 315069;
B.seq4= allign_eit_and_perf(B.seq4);
% B.seq4.eit.fdata= lowpass_iir(B.seq4.eit.data, lpass, eitfs);
B.seq4.imgr= get_imgr(B.seq4, imdl);
D= B;

% 10.2
elseif pig== "10.2"
[vv, aux]= eidors_readdata('EIT_nativ_10.2._schleuse.eit');
C.seq1.eit.data= vv;
C.seq1.eit.elec_impedance= aux.elec_impedance;
C.seq1.name= char(remove_underscores('EIT_nativ_10.2._schleuse.eit'));
C.seq1.pig= "10.2";
C.seq1.eit.fs= eitfs;
C.seq1.eit.sync1= 399;
C.seq1.eit.sync2= 5232;
C.seq1.eit.apn= C.seq1.eit.sync1 + round(7.4* eitfs);
C.seq1.eit.inj= C.seq1.eit.sync1 + round(13.4* eitfs);
C.seq1.eit.vnt= C.seq1.eit.sync1 + round(79.1* eitfs);
C.seq1.perf= load('EIT_nativ_10.2._schleuse.mat');
C.seq1.perf.apn= 13720;
C.seq1.perf.inj= 19708;
C.seq1.perf.vnt= 85388;
C.seq1= allign_eit_and_perf(C.seq1);
% C.seq1.eit.fdata= lowpass_iir(C.seq1.eit.data, lpass, eitfs);
C.seq1.imgr= get_imgr(C.seq1, imdl);

[vv, aux]= eidors_readdata('EIT_nativ_2_zvk_10.2.eit');
C.seq2.eit.data= vv;
C.seq2.eit.elec_impedance= aux.elec_impedance;
C.seq2.name= char(remove_underscores('EIT_nativ_2_zvk_10.2.eit'));
C.seq2.pig= "10.2";
C.seq2.eit.fs= eitfs;
C.seq2.eit.sync1= 262;
C.seq2.eit.sync2= 6729;
C.seq2.eit.apn= C.seq2.eit.sync1 + round(9.3* eitfs);
C.seq2.eit.inj= C.seq2.eit.sync1 + round(15.6* eitfs);
C.seq2.eit.vnt= C.seq2.eit.sync1 + round(85.5* eitfs);
C.seq2.perf= load('EIT_nativ_2_zvk_10.2.mat');
C.seq2.perf.apn= 12336;
C.seq2.perf.inj= 18832;
C.seq2.perf.vnt= 88633;
C.seq2= allign_eit_and_perf(C.seq2);
% C.seq2.eit.fdata= lowpass_iir(C.seq2.eit.data, lpass, eitfs);
C.seq2.imgr= get_imgr(C.seq2, imdl);

[vv, aux]= eidors_readdata('EIT_nach_Perfusionsminderung_Sequenz_3_Schleuse_10.2.eit');
C.seq3.eit.data= vv;
C.seq3.eit.elec_impedance= aux.elec_impedance;
C.seq3.name= char(remove_underscores('EIT_nach_Perfusionsminderung_Sequenz_3_Schleuse_10.2.eit'));
C.seq3.pig= "10.2";
C.seq3.eit.fs= eitfs;
C.seq3.eit.sync1= 124;
C.seq3.eit.sync2= 4866;
C.seq3.eit.apn= C.seq3.eit.sync1 + round(12.8* eitfs);
C.seq3.eit.inj= C.seq3.eit.sync1 + round(18.6* eitfs);
C.seq3.eit.vnt= C.seq3.eit.sync1 + round(73* eitfs);
C.seq3.perf= load('EIT_nach_Perfusionsminderung_Sequenz_3_Schleuse_10.2.mat');
C.seq3.perf.apn= 17460;
C.seq3.perf.inj= 23357;
C.seq3.perf.vnt= 77865;
C.seq3= allign_eit_and_perf(C.seq3);
% C.seq3.eit.fdata= lowpass_iir(C.seq3.eit.data, lpass, eitfs);
C.seq3.imgr= get_imgr(C.seq3, imdl);

[vv, aux]= eidors_readdata('EIT_nach_Perfusionsminderung_Sequenz_4_ZVK.eit');
C.seq4.eit.data= vv;
C.seq4.eit.elec_impedance= aux.elec_impedance;
C.seq4.name= char(remove_underscores('EIT_nach_Perfusionsminderung_Sequenz_4_ZVK.eit'));
C.seq4.pig= "10.2";
C.seq4.eit.fs= eitfs;
C.seq4.eit.sync1= 149;
C.seq4.eit.sync2= 5390;
C.seq4.eit.apn= C.seq4.eit.sync1 + round(11.9* eitfs);
C.seq4.eit.inj= C.seq4.eit.sync1 + round(21.9* eitfs);
C.seq4.eit.vnt= C.seq4.eit.sync1 + round(88.5* eitfs);
C.seq4.perf= load('EIT_nach_Perfusionsminderung_Sequenz_4_ZVK.mat');
C.seq4.perf.apn= 13854;
C.seq4.perf.inj= 23835;
C.seq4.perf.vnt= 90452;
C.seq4= allign_eit_and_perf(C.seq4);
% C.seq4.eit.fdata= lowpass_iir(C.seq4.eit.data, lpass, eitfs);
C.seq4.imgr= get_imgr(C.seq4, imdl);

[vv, aux]= eidors_readdata('EIT_10.2._nach_6h_Perfusionsmin_Schleuse Sequ5.eit');
C.seq5.eit.data= vv;
C.seq5.eit.elec_impedance= aux.elec_impedance;
C.seq5.name= char(remove_underscores('EIT_10.2._nach_6h_Perfusionsmin_Schleuse Sequ5.eit'));
C.seq5.pig= "10.2";
C.seq5.eit.fs= eitfs;
C.seq5.eit.sync1= 177;
C.seq5.eit.sync2= 5558;
C.seq5.eit.apn= C.seq5.eit.sync1 + round(14* eitfs);
C.seq5.eit.inj= C.seq5.eit.sync1 + round(20.1* eitfs);
C.seq5.eit.vnt= C.seq5.eit.sync1 + round(87.7* eitfs);
C.seq5.perf= load('EIT_10.2._nach_6h_Perfusionsmin_Schleuse Sequ5.mat');
C.seq5.perf.apn= 30897;
C.seq5.perf.inj= 37773;
C.seq5.perf.vnt= 104515;
C.seq5= allign_eit_and_perf(C.seq5);
% C.seq5.eit.fdata= lowpass_iir(C.seq5.eit.data, lpass, eitfs);
C.seq5.imgr= get_imgr(C.seq5, imdl);

[vv, aux]= eidors_readdata('EIT_nach_6h_Perfusionsminderung_10.02 ZVK Sequ6.eit');
C.seq6.eit.data= vv;
C.seq6.eit.elec_impedance= aux.elec_impedance;
C.seq6.name= char(remove_underscores('EIT_nach_6h_Perfusionsminderung_10.02 ZVK Sequ6.eit'));
C.seq6.pig= "10.2";
C.seq6.eit.fs= eitfs;
C.seq6.eit.sync1= 191;
C.seq6.eit.sync2= 4995;
C.seq6.eit.apn= C.seq6.eit.sync1 + round(14.5* eitfs);
C.seq6.eit.inj= C.seq6.eit.sync1 + round(20.4* eitfs);
C.seq6.eit.vnt= C.seq6.eit.sync1 + round(81.2* eitfs);
C.seq6.perf= load('EIT_nach_6h_Perfusionsminderung_10.02 ZVK Sequ6.mat');
C.seq6.perf.apn= 15381;
C.seq6.perf.inj= 21370;
C.seq6.perf.vnt= 82125;
C.seq6= allign_eit_and_perf(C.seq6);
% C.seq6.eit.fdata= lowpass_iir(C.seq6.eit.data, lpass, eitfs);
C.seq6.imgr= get_imgr(C.seq6, imdl);
D= C;

% 11.2
elseif pig== "11.2"
D.seq1.eit.data= eidors_readdata('EIT_11.2_Nativ_1_Schleuse.eit');
D.seq1.name= char(remove_underscores('EIT_11.2_Nativ_1_Schleuse.eit'));
D.seq1.pig= "11.2";
D.seq1.eit.fs= eitfs;
D.seq1.eit.sync1= 1225;
D.seq1.eit.sync2= 6369;
D.seq1.eit.apn= D.seq1.eit.sync1 + round(-10.2* eitfs);
D.seq1.eit.inj= D.seq1.eit.sync1 + round(11.5* eitfs);
D.seq1.eit.vnt= D.seq1.eit.sync1 + round(89.1* eitfs);
D.seq1.perf= load('EIT_11.2_Nativ_1_Schleuse.mat');
D.seq1.perf.apn= 11690;
D.seq1.perf.inj= 33376;
D.seq1.perf.vnt= 111004;
D.seq1= allign_eit_and_perf(D.seq1);
% D.seq1.eit.fdata= lowpass_iir(D.seq1.eit.data, lpass, eitfs);
D.seq1.imgr= get_imgr(D.seq1, imdl);

D.seq2.eit.data= eidors_readdata('EIT_11.2_Nativ_2_ZVK.eit');
D.seq2.name= char(remove_underscores('EIT_11.2_Nativ_2_ZVK.eit'));
D.seq2.pig= "11.2";
D.seq2.eit.fs= eitfs;
D.seq2.eit.sync1= 211;
D.seq2.eit.sync2= 5694;
D.seq2.eit.apn= D.seq2.eit.sync1 + round(10.7* eitfs);
D.seq2.eit.inj= D.seq2.eit.sync1 + round(15.3* eitfs);
D.seq2.eit.vnt= D.seq2.eit.sync1 + round(91.9* eitfs);
D.seq2.perf= load('EIT_11.2_Nativ_2_ZVK.mat');
D.seq2.perf.apn= 10714;
D.seq2.perf.inj= 15269;
D.seq2.perf.vnt= 91873;
D.seq2= allign_eit_and_perf(D.seq2);
% D.seq2.eit.fdata= lowpass_iir(D.seq2.eit.data, lpass, eitfs);
D.seq2.imgr= get_imgr(D.seq2, imdl);

D.seq3.eit.data= eidors_readdata('EIT_Sequenz_3_Schleuse.eit');
D.seq3.name= char(remove_underscores('EIT_Sequenz_3_Schleuse.eit'));
D.seq3.pig= "11.2";
D.seq3.eit.fs= eitfs;
D.seq3.eit.sync1= 261;
D.seq3.eit.sync2= 7008;
D.seq3.eit.apn= D.seq3.eit.sync1 + round(25.1* eitfs);
D.seq3.eit.inj= D.seq3.eit.sync1 + round(32.2* eitfs);
D.seq3.eit.vnt= D.seq3.eit.sync1 + round(113.7* eitfs);
D.seq3.perf= load('EIT_Sequenz_3_Schleuse.mat');
D.seq3.perf.apn= 36930;
D.seq3.perf.inj= 43994;
D.seq3.perf.vnt= 125489; % inj + 81495
D.seq3= allign_eit_and_perf(D.seq3);
% D.seq3.eit.fdata= lowpass_iir(D.seq3.eit.data, lpass, eitfs);
D.seq3.imgr= get_imgr(D.seq3, imdl);

D.seq4.eit.data= eidors_readdata('EIT_Seqeunz4_zvk.eit');
D.seq4.name= char(remove_underscores('EIT_Seqeunz4_zvk.eit'));
D.seq4.pig= "11.2";
D.seq4.eit.fs= eitfs;
D.seq4.eit.sync1= 55;
D.seq4.eit.sync2= 5814;
D.seq4.eit.apn= D.seq4.eit.sync1 + round(22.3* eitfs);
D.seq4.eit.inj= D.seq4.eit.sync1 + round(66* eitfs);
D.seq4.eit.vnt= D.seq4.eit.sync1 + round(104.4* eitfs);
D.seq4.perf= load('EIT_Seqeunz4_zvk.mat');
D.seq4.perf.apn= 31092;
D.seq4.perf.inj= 74955;
D.seq4.perf.vnt= 113239;
D.seq4= allign_eit_and_perf(D.seq4);
% D.seq4.eit.fdata= lowpass_iir(D.seq4.eit.data, lpass, eitfs);
D.seq4.imgr= get_imgr(D.seq4, imdl);

D.seq5.eit.data= eidors_readdata('EIT_11.2_4h_nach_Stroke_Schleuse_Sequ5.eit');
D.seq5.name= char(remove_underscores('EIT_11.2_4h_nach_Stroke_Schleuse_Sequ5.eit'));
D.seq5.pig= "11.2";
D.seq5.eit.fs= eitfs;
D.seq5.eit.sync1= 495;
D.seq5.eit.sync2= 5322;
D.seq5.eit.apn= D.seq5.eit.sync1 + round(8.1* eitfs);
D.seq5.eit.inj= D.seq5.eit.sync1 + round(18.2* eitfs);
D.seq5.eit.vnt= D.seq5.eit.sync1 + round(80.4* eitfs);
D.seq5.perf= load('EIT_11.2_4h_nach_Stroke_Schleuse_Sequ5.mat');
D.seq5.perf.apn= 16000;
D.seq5.perf.inj= 26149;
D.seq5.perf.vnt= 88339;
D.seq5= allign_eit_and_perf(D.seq5);
% D.seq5.eit.fdata= lowpass_iir(D.seq5.eit.data, lpass, eitfs);
D.seq5.imgr= get_imgr(D.seq5, imdl);

D.seq6.eit.data= eidors_readdata('EIT_4h_nach_Stroke_ZVk_Sequ6.eit');
D.seq6.name= char(remove_underscores('EIT_4h_nach_Stroke_ZVk_Sequ6.eit'));
D.seq6.pig= "11.2";
D.seq6.eit.fs= eitfs;
D.seq6.eit.sync1= 169;
D.seq6.eit.sync2= 6666;
D.seq6.eit.apn= D.seq6.eit.sync1 + round(12.1* eitfs);
D.seq6.eit.inj= D.seq6.eit.sync1 + round(16.2* eitfs);
D.seq6.eit.vnt= D.seq6.eit.sync1 + round(85.8* eitfs);
D.seq6.perf= load('EIT_4h_nach_Stroke_ZVk_Sequ6.mat');
D.seq6.perf.apn= 19419;
D.seq6.perf.inj= 23516;
D.seq6.perf.vnt= 93187;
D.seq6= allign_eit_and_perf(D.seq6);
% D.seq6.eit.fdata= lowpass_iir(D.seq6.eit.data, lpass, eitfs);
D.seq6.imgr= get_imgr(D.seq6, imdl);


% 12.2
elseif pig== "12.2"
E.seq1.eit.data= eidors_readdata('EIT_12.2. Seqeunz 1 nativ Schleuse.eit');
E.seq1.name= char(remove_underscores('EIT_12.2. Seqeunz 1 nativ Schleuse.eit'));
E.seq1.pig= "12.2";
E.seq1.eit.fs= eitfs;
E.seq1.eit.sync1= 251;
E.seq1.eit.sync2= 5754;
E.seq1.eit.apn= E.seq1.eit.sync1 + round(8.5* eitfs);
E.seq1.eit.inj= E.seq1.eit.sync1 + round(15* eitfs);
E.seq1.eit.vnt= E.seq1.eit.sync1 + round(94.7* eitfs);
E.seq1.perf= load('EIT 12.2 sequence 3.mat');
E.seq1.perf.apn= 10415;
E.seq1.perf.inj= 16873;
E.seq1.perf.vnt= 96616;
E.seq1= allign_eit_and_perf(E.seq1);
% E.seq1.eit.fdata= lowpass_iir(E.seq1.eit.data, lpass, eitfs);
E.seq1.imgr= get_imgr(E.seq1, imdl);

E.seq2.eit.data= eidors_readdata('EIT_12.2_NATIV_ZVK_Sequ_2.eit');
E.seq2.name= char(remove_underscores('EIT_12.2_NATIV_ZVK_Sequ_2.eit'));
E.seq2.pig= "12.2";
E.seq2.eit.fs= eitfs;
E.seq2.eit.sync1= 156;
E.seq2.eit.sync2= 7179;
E.seq2.eit.apn= E.seq2.eit.sync1 + round(14.5* eitfs);
E.seq2.eit.inj= E.seq2.eit.sync1 + round(23.4* eitfs);
E.seq2.eit.vnt= E.seq2.eit.sync1 + round(117* eitfs);
E.seq2.perf= load('EIT 12.2 sequence 4.mat');
E.seq2.perf.apn= 15490;
E.seq2.perf.inj= 24463;
E.seq2.perf.vnt= 116350; % 24463+ 91887
E.seq2= allign_eit_and_perf(E.seq2);
% E.seq2.eit.fdata= lowpass_iir(E.seq2.eit.data, lpass, eitfs);
E.seq2.imgr= get_imgr(E.seq2, imdl);

E.seq3.eit.data= eidors_readdata('EIT_direkt_nach_Stroke_12.2._Sequ_3_Schleuse.eit');
E.seq3.name= char(remove_underscores('EIT_direkt_nach_Stroke_12.2._Sequ_3_Schleuse.eit'));
E.seq3.pig= "12.2";
E.seq3.eit.fs= eitfs;
E.seq3.eit.sync1= 126;
E.seq3.eit.sync2= 7600;
E.seq3.eit.apn= E.seq3.eit.sync1 + round(12.7* eitfs);
E.seq3.eit.inj= E.seq3.eit.sync1 + round(50* eitfs);
E.seq3.eit.vnt= E.seq3.eit.sync1 + round(128* eitfs);
E.seq3.perf= load('EIT 12.2 sequence 5.mat');
E.seq3.perf.apn= 15807;
E.seq3.perf.inj= 53138;
E.seq3.perf.vnt= 131095;
E.seq3= allign_eit_and_perf(E.seq3);
% E.seq3.eit.fdata= lowpass_iir(E.seq3.eit.data, lpass, eitfs);
E.seq3.imgr= get_imgr(E.seq3, imdl);

E.seq4.eit.data= eidors_readdata('EIT_12.2_Sequ_4_nach_Stroke_ZVK.eit');
E.seq4.name= char(remove_underscores('EIT_12.2_Sequ_4_nach_Stroke_ZVK.eit'));
E.seq4.pig= "12.2";
E.seq4.eit.fs= eitfs;
E.seq4.eit.sync1= 518;
E.seq4.eit.sync2= 7420;
E.seq4.eit.apn= E.seq4.eit.sync1 + round(13* eitfs);
E.seq4.eit.inj= E.seq4.eit.sync1 + round(36.2* eitfs);
E.seq4.eit.vnt= E.seq4.eit.sync1 + round(126.6* eitfs);
E.seq4.perf= load('EIT 12.2 sequence 6.mat');
E.seq4.perf.apn= 22384;
E.seq4.perf.inj= 45609;
E.seq4.perf.vnt= 136047;
E.seq4= allign_eit_and_perf(E.seq4);
% E.seq4.eit.fdata= lowpass_iir(E.seq4.eit.data, lpass, eitfs);
E.seq4.imgr= get_imgr(E.seq4, imdl);

E.seq5.eit.data= eidors_readdata('EIT_12.2._3,5_nach_Stroke_Sequ_5_Schleuse.eit');
E.seq5.name= char(remove_underscores('EIT_12.2._3,5_nach_Stroke_Sequ_5_Schleuse.eit'));
E.seq5.pig= "12.2";
E.seq5.eit.fs= eitfs;
E.seq5.eit.sync1= 132;
E.seq5.eit.sync2= 6586;
E.seq5.eit.apn= E.seq5.eit.sync1 + round(17.4* eitfs);
E.seq5.eit.inj= E.seq5.eit.sync1 + round(37.5* eitfs);
E.seq5.eit.vnt= E.seq5.eit.sync1 + round(107.6* eitfs);
E.seq5.perf= load('EIT 12.2 sequence 7.mat');
E.seq5.perf.apn= 24336;
E.seq5.perf.inj= 44446;
E.seq5.perf.vnt= 114530;
E.seq5= allign_eit_and_perf(E.seq5);
% E.seq5.eit.fdata= lowpass_iir(E.seq5.eit.data, lpass, eitfs);
E.seq5.imgr= get_imgr(E.seq5, imdl);

E.seq6.eit.data= eidors_readdata('EIT_12.2._Sequ_6_3,5h_nach_Stroke_ZVK.eit');
E.seq6.name= char(remove_underscores('EIT_12.2._Sequ_6_3,5h_nach_Stroke_ZVK.eit'));
E.seq6.pig= "12.2";
E.seq6.eit.fs= eitfs;
E.seq6.eit.sync1= 245;
E.seq6.eit.sync2= 7715;
E.seq6.eit.apn= E.seq6.eit.sync1 + round(12.6* eitfs);
E.seq6.eit.inj= E.seq6.eit.sync1 + round(19.9* eitfs);
E.seq6.eit.vnt= E.seq6.eit.sync1 + round(117* eitfs);
E.seq6.perf= load('EIT 12.2 sequence 8.mat');
E.seq6.perf.apn= 12550;
E.seq6.perf.inj= 19840;
E.seq6.perf.vnt= 116881;
E.seq6= allign_eit_and_perf(E.seq6);
% E.seq6.eit.fdata= lowpass_iir(E.seq6.eit.data, lpass, eitfs);
E.seq6.imgr= get_imgr(E.seq6, imdl);
D= E;
end % end if

end % end function

function seq= allign_eit_and_perf(seq)
% Find optimal alignment between eit and perfusion data
% Annotate peaks and troughs of perfusion signal and transfer annotations
% to EIT data
valDiffThresh= 400;

estart= seq.eit.sync1+ round(seq.eit.fs); % desired starting frame is 1 s after first sync
eend= seq.eit.sync2- round(seq.eit.fs); % desired end frame is 1 s before last sync
pstart= round((estart/ seq.eit.fs)* seq.perf.tickrate);
pend= round((eend/ seq.eit.fs)* seq.perf.tickrate);

eapn= seq.eit.apn/ seq.eit.fs;
einj= seq.eit.inj/ seq.eit.fs;
evnt= seq.eit.vnt/ seq.eit.fs;
papn= seq.perf.apn/ seq.perf.tickrate;
pinj= seq.perf.inj/ seq.perf.tickrate;
pvnt= seq.perf.vnt/ seq.perf.tickrate;

t1= eapn- papn;
t2= einj- pinj;
t3= evnt- pvnt;
shift= mean( [t1, t2, t3] );

if shift> 0 % EIT sequence is ahead of perfusion sequence
    ecut= estart+ round( shift*  seq.eit.fs);
    pcut= pstart;   
elseif shift< 0 % Perfusion sequence is ahead of EIT sequence
    ecut= estart;
    pcut= pstart+ abs(round( shift*  seq.perf.tickrate));
end % end if
% filtered data
lowpass_cutoff= 8;

% new filter (hilbert, moving mean)
seq.eit.fdata= lowpass_iir(seq.eit.data, lowpass_cutoff, seq.eit.fs);
seq.eit.fdata= seq.eit.fdata(:, ecut+1: eend);

seq.eit.data= seq.eit.data(:, ecut+1: eend);

seq.eit.sync1= seq.eit.sync1- ecut;
seq.eit.sync2= seq.eit.sync2- ecut;
seq.eit.apn= seq.eit.apn- ecut;
seq.eit.inj= seq.eit.inj- ecut;
seq.eit.vnt= seq.eit.vnt- ecut;

seq.perf.data= movmean(seq.perf.data(:, pcut+1: pend+ pcut), 10);
seq.perf.apn= seq.perf.apn- pcut;
seq.perf.inj= seq.perf.inj- pcut;
seq.perf.vnt= seq.perf.vnt- pcut;

% Find peaks and valleys in perfusion Signal
seq.perf.vals= peakfinder(seq.perf.data, 8, [], -1);

% Remove false valleys that are too close to true valleys
for i= 1:length(seq.perf.vals)- 1
    if seq.perf.vals(i+1)- seq.perf.vals(i) < valDiffThresh 
        seq.perf.vals(i+1)= seq.perf.vals(i);
    end % end if
end % end for

% Remove valleys within the synchronization flush
invalids= [seq.perf.apn- 100, seq.perf.inj- 100, seq.perf.vnt- 100];
for invalid= [seq.perf.apn, seq.perf.inj, seq.perf.vnt]
    rm= find( (seq.perf.vals> invalid) + (seq.perf.vals< invalid+ 1100)==2 );
    if ~isempty(rm)
        seq.perf.vals(rm)= 0;
    end % end if
end % end if

seq.perf.vals= seq.perf.vals(seq.perf.vals~= 0);
seq.perf.vals= unique(seq.perf.vals);
seq.perf.peaks= zeros(1, length(seq.perf.vals)- 1);

for i= 1:length(seq.perf.peaks)
    start= seq.perf.vals(i);
    stop= seq.perf.vals(i+1);
    loc= find( seq.perf.data== max(seq.perf.data(start: stop)));
    seq.perf.peaks(i)= loc(find((loc> start) + (loc< stop)==2, 1));
end % end for

% Remove peaks within the synchronization flush
for invalid= invalids
    rm= find( (seq.perf.peaks> invalid) + (seq.perf.peaks< invalid+ 1100)==2 );
    if ~isempty(rm)
        seq.perf.peaks(rm)= 0;
    end % end if
end % end if

seq.perf.peaks= seq.perf.peaks(seq.perf.peaks~= 0);

% Transfer these annotations to the EIT data
seq.eit.peaks= round( (seq.perf.peaks./ seq.perf.tickrate).*seq.eit.fs );
seq.eit.vals= round( (seq.perf.vals./ seq.perf.tickrate).*seq.eit.fs );

end % end function

function imgr= get_imgr(seq, imdl)

vv= real(seq.eit.fdata);

% Remove bad electrodes
vv_cleaned= rm_bad_elec(vv, imdl);

% meas_std= std(vv,{},2);
% lim= mean(meas_std)+ std(meas_std); % remove measurements more than 1 sd of mean
% noisy= std(vv,{},2) > lim;
% vv(noisy,:)= 0;
% % code for testing noise removal
% sorted= sort(std(vv,{},2));
% histogram(sorted);
% title("Distribution of cleaned measurement pair SD for D.seq1 w new filter", 'FontSize', 18);
% xlabel("standard deviation of measurement pair", 'FontSize', 12);
% ylabel("frequency", 'FontSize', 12);
% length(find(noisy))

vh= real( vv_cleaned(:,seq.eit.inj)); % reference frame is injection frame
% vh = mean( real( seq.eit.fdata(:,:)),2); % reference frame is mean of
% whole time series

imgr= inv_solve(imdl, vh, vv);

% Set colour limits
clim= max(imgr.elem_data(:));
imgr.calc_colours.ref_level= 0;
imgr.calc_colours.lim= clim;

end % end function

function varargout = peakfinder(x0, sel, thresh, extrema, includeEndpoints, interpolate)
% https://www.mathworks.com/matlabcentral/fileexchange/25500-peakfinder-x0-sel-thresh-extrema-includeendpoints-interpolate
%PEAKFINDER Noise tolerant fast peak finding algorithm
%   INPUTS:
%       x0 - A real vector from the maxima will be found (required)
%       sel - The amount above surrounding data for a peak to be,
%           identified (default = (max(x0)-min(x0))/4). Larger values mean
%           the algorithm is more selective in finding peaks.
%       thresh - A threshold value which peaks must be larger than to be
%           maxima or smaller than to be minima.
%       extrema - 1 if maxima are desired, -1 if minima are desired
%           (default = maxima, 1)
%       includeEndpoints - If true the endpoints will be included as
%           possible extrema otherwise they will not be included
%           (default = true)
%       interpolate - If true quadratic interpolation will be performed
%           around each extrema to estimate the magnitude and the
%           position of the peak in terms of fractional indicies. Note that
%           unlike the rest of this function interpolation assumes the
%           input is equally spaced. To recover the x_values of the input
%           rather than the fractional indicies you can do:
%           peakX = x0 + (peakLoc - 1) * dx
%           where x0 is the first x value and dx is the spacing of the
%           vector. Output peakMag to recover interpolated magnitudes.
%           See example 2 for more information.
%           (default = false)
%
%   OUTPUTS:
%       peakLoc - The indicies of the identified peaks in x0
%       peakMag - The magnitude of the identified peaks
%
%   [peakLoc] = peakfinder(x0) returns the indicies of local maxima that
%       are at least 1/4 the range of the data above surrounding data.
%
%   [peakLoc] = peakfinder(x0,sel) returns the indicies of local maxima
%       that are at least sel above surrounding data.
%
%   [peakLoc] = peakfinder(x0,sel,thresh) returns the indicies of local
%       maxima that are at least sel above surrounding data and larger
%       (smaller) than thresh if you are finding maxima (minima).
%
%   [peakLoc] = peakfinder(x0,sel,thresh,extrema) returns the maxima of the
%       data if extrema > 0 and the minima of the data if extrema < 0
%
%   [peakLoc] = peakfinder(x0,sel,thresh,extrema, includeEndpoints)
%       returns the endpoints as possible extrema if includeEndpoints is
%       considered true in a boolean sense
%
%   [peakLoc, peakMag] = peakfinder(x0,sel,thresh,extrema,interpolate)
%       returns the results of results of quadratic interpolate around each
%       extrema if interpolate is considered to be true in a boolean sense
%
%   [peakLoc, peakMag] = peakfinder(x0,...) returns the indicies of the
%       local maxima as well as the magnitudes of those maxima
%
%   If called with no output the identified maxima will be plotted along
%       with the input data.
%
%   Note: If repeated values are found the first is identified as the peak
%
% Example 1:
% t = 0:.0001:10;
% x = 12*sin(10*2*pi*t)-3*sin(.1*2*pi*t)+randn(1,numel(t));
% x(1250:1255) = max(x);
% peakfinder(x)
%
% Example 2:
% ds = 100;  % Downsample factor
% dt = .001; % Time step
% ds_dt = ds*dt; % Time delta after downsampling
% t0 = 1;
% t = t0:dt:5 + t0;
% x = 0.2-sin(0.01*2*pi*t)+3*cos(7/13*2*pi*t+.1)-2*cos((1+pi/10)*2*pi*t+0.2)-0.2*t;
% x(end) = min(x);
% x_ds = x(1:ds:end); % Downsample to test interpolation
% [minLoc, minMag] = peakfinder(x_ds, .8, 0, -1, false, true);
% minT = t0 + (minLoc - 1) * ds_dt; % Take into account 1 based indexing
% p = plot(t,x,'-',t(1:ds:end),x_ds,'o',minT,minMag,'rv');
% set(p(2:end), 'linewidth', 2); % Show the markers more clearly
% legend('Actual Data', 'Input Data', 'Estimated Peaks');
% Copyright Nathanael C. Yoder 2015 (nyoder@gmail.com)
% Perform error checking and set defaults if not passed in
narginchk(1, 6);
nargoutchk(0, 2);
s = size(x0);
flipData =  s(1) < s(2);
len0 = numel(x0);
if len0 ~= s(1) && len0 ~= s(2)
    error('PEAKFINDER:Input','The input data must be a vector')
elseif isempty(x0)
    varargout = {[],[]};
    return;
end % end if
if ~isreal(x0)
    warning('PEAKFINDER:NotReal','Absolute value of data will be used')
    x0 = abs(x0);
end % end if
if nargin < 2 || isempty(sel)
    sel = (max(x0)-min(x0))/4;
elseif ~isnumeric(sel) || ~isreal(sel)
    sel = (max(x0)-min(x0))/4;
    warning('PEAKFINDER:InvalidSel',...
        'The selectivity must be a real scalar.  A selectivity of %.4g will be used',sel)
elseif numel(sel) > 1
    warning('PEAKFINDER:InvalidSel',...
        'The selectivity must be a scalar.  The first selectivity value in the vector will be used.')
    sel = sel(1);
end % end if
if nargin < 3 || isempty(thresh)
    thresh = [];
elseif ~isnumeric(thresh) || ~isreal(thresh)
    thresh = [];
    warning('PEAKFINDER:InvalidThreshold',...
        'The threshold must be a real scalar. No threshold will be used.')
elseif numel(thresh) > 1
    thresh = thresh(1);
    warning('PEAKFINDER:InvalidThreshold',...
        'The threshold must be a scalar.  The first threshold value in the vector will be used.')
end % end if
if nargin < 4 || isempty(extrema)
    extrema = 1;
else
    extrema = sign(extrema(1)); % Should only be 1 or -1 but make sure
    if extrema == 0
        error('PEAKFINDER:ZeroMaxima','Either 1 (for maxima) or -1 (for minima) must be input for extrema');
    end % end if
end % end if
if nargin < 5 || isempty(includeEndpoints)
    includeEndpoints = true;
end
if nargin < 6 || isempty(interpolate)
    interpolate = false;
end % end if
x0 = extrema*x0(:); % Make it so we are finding maxima regardless
thresh = thresh*extrema; % Adjust threshold according to extrema.
dx0 = diff(x0); % Find derivative
dx0(dx0 == 0) = -eps; % This is so we find the first of repeated values
ind = find(dx0(1:end-1).*dx0(2:end) < 0)+1; % Find where the derivative changes sign
% Include endpoints in potential peaks and valleys as desired
if includeEndpoints
    x = [x0(1);x0(ind);x0(end)];
    ind = [1;ind;len0];
    minMag = min(x);
    leftMin = minMag;
else
    x = x0(ind);
    minMag = min(x);
    leftMin = min(x(1), x0(1));
end % end if
% x only has the peaks, valleys, and possibly endpoints
len = numel(x);
if len > 2 % Function with peaks and valleys
    % Set initial parameters for loop
    tempMag = minMag;
    foundPeak = false;
    if includeEndpoints
        % Deal with first point a little differently since tacked it on
        % Calculate the sign of the derivative since we tacked the first
        %  point on it does not neccessarily alternate like the rest.
        signDx = sign(diff(x(1:3)));
        if signDx(1) <= 0 % The first point is larger or equal to the second
            if signDx(1) == signDx(2) % Want alternating signs
                x(2) = [];
                ind(2) = [];
                len = len-1;
            end
        else % First point is smaller than the second
            if signDx(1) == signDx(2) % Want alternating signs
                x(1) = [];
                ind(1) = [];
                len = len-1;
            end % end if
        end % end if
    end % end if
    % Skip the first point if it is smaller so we always start on a
    %   maxima
    if x(1) >= x(2)
        ii = 0;
    else
        ii = 1;
    end % end if
    % Preallocate max number of maxima
    maxPeaks = ceil(len/2);
    peakLoc = zeros(maxPeaks,1);
    peakMag = zeros(maxPeaks,1);
    cInd = 1;
    % Loop through extrema which should be peaks and then valleys
    while ii < len
        ii = ii+1; % This is a peak
        % Reset peak finding if we had a peak and the next peak is bigger
        %   than the last or the left min was small enough to reset.
        if foundPeak
            tempMag = minMag;
            foundPeak = false;
        end % end if
        % Found new peak that was lager than temp mag and selectivity larger
        %   than the minimum to its left.
        if x(ii) > tempMag && x(ii) > leftMin + sel
            tempLoc = ii;
            tempMag = x(ii);
        end % end if
        % Make sure we don't iterate past the length of our vector
        if ii == len
            break; % We assign the last point differently out of the loop
        end % end if
        ii = ii+1; % Move onto the valley
        % Come down at least sel from peak
        if ~foundPeak && tempMag > sel + x(ii)
            foundPeak = true; % We have found a peak
            leftMin = x(ii);
            peakLoc(cInd) = tempLoc; % Add peak to index
            peakMag(cInd) = tempMag;
            cInd = cInd+1;
        elseif x(ii) < leftMin % New left minima
            leftMin = x(ii);
        end % end if
    end % end while
    % Check end point
    if includeEndpoints
        if x(end) > tempMag && x(end) > leftMin + sel
            peakLoc(cInd) = len;
            peakMag(cInd) = x(end);
            cInd = cInd + 1;
        elseif ~foundPeak && tempMag > minMag % Check if we still need to add the last point
            peakLoc(cInd) = tempLoc;
            peakMag(cInd) = tempMag;
            cInd = cInd + 1;
        end % end if
    elseif ~foundPeak
        if x(end) > tempMag && x(end) > leftMin + sel
            peakLoc(cInd) = len;
            peakMag(cInd) = x(end);
            cInd = cInd + 1;
        elseif tempMag > min(x0(end), x(end)) + sel
            peakLoc(cInd) = tempLoc;
            peakMag(cInd) = tempMag;
            cInd = cInd + 1;
        end % end if
    end % end if
    % Create output
    if cInd > 1
        peakInds = ind(peakLoc(1:cInd-1));
        peakMags = peakMag(1:cInd-1);
    else
        peakInds = [];
        peakMags = [];
    end % end if
else % This is a monotone function where an endpoint is the only peak
    [peakMags,xInd] = max(x);
    if includeEndpoints && peakMags > minMag + sel
        peakInds = ind(xInd);
    else
        peakMags = [];
        peakInds = [];
    end % end if
end % end if
% Apply threshold value.  Since always finding maxima it will always be
%   larger than the thresh.
if ~isempty(thresh)
    m = peakMags>thresh;
    peakInds = peakInds(m);
    peakMags = peakMags(m);
end % end if
if interpolate && ~isempty(peakMags)
    middleMask = (peakInds > 1) & (peakInds < len0);
    noEnds = peakInds(middleMask);
    magDiff = x0(noEnds + 1) - x0(noEnds - 1);
    magSum = x0(noEnds - 1) + x0(noEnds + 1)  - 2 * x0(noEnds);
    magRatio = magDiff ./ magSum;
    peakInds(middleMask) = peakInds(middleMask) - magRatio/2;
    peakMags(middleMask) = peakMags(middleMask) - magRatio .* magDiff/8;
end % end if
% Rotate data if needed
if flipData
    peakMags = peakMags.';
    peakInds = peakInds.';
end % end if
% Change sign of data if was finding minima
if extrema < 0
    peakMags = -peakMags;
    x0 = -x0;
end % end if
% Plot if no output desired
if nargout == 0
    if isempty(peakInds)
        disp('No significant peaks found')
    else
        figure;
        plot(1:len0,x0,'.-',peakInds,peakMags,'ro','linewidth',2);
    end % end if
else
    varargout = {peakInds,peakMags};
end % end if
end % end function
