function [m_spec_deriv , m_AM, m_FM ,m_Entropy , m_amplitude ,gravity_center, m_PitchGoodness , m_Pitch , Pitch_chose , Pitch_weight ]=deriv(TS,fs);


% [m_spec_deriv , m_AM, m_FM ,m_Entropy , m_amplitude ,gravity_center, m_PitchGoodness , m_Pitch , Pitch_chose , Pitch_weight ]=deriv(TS,fs);
%	 uses the derivative estimates to calculate the spectrum's frequency
%	 and time derivatives
% 
%        S: estimated spectrum; S_f: estimated frequency derivative; 
%        S_t: estimated time derivative
%        NW: time bandwidth parameter (e.g. 3)
%        K : number of tapers kept, approx. 2*NW-1
%        pad: length to which data will be padded (preferably power of 2
%        window: time window size
%        winstep: distance between centers of adjacent time windows

%        Written by Sigal Saar August 08 2005







TS=filter_sound_sam(TS);

 load('parameters');

% Parameters;

E=taper_read();
N=length(TS);
%if N>300000
%    TS_all=TS;
    
TSM=running_windows(TS',param.window, param.winstep);
%TSM=running_windows(TS',param.window, 44.1);
S=0;
SF=0;
if floor(param.winstep)~=param.winstep
    E=[E ; [0 0]];
end
    J1=(fft(TSM(:,:).*(ones(size(TSM,1),1)*(E(:,1))'),param.pad,2));
    J1=J1(:,1:param.spectrum_range)* ( 27539);
    J2=(fft(TSM(:,:).*(ones(size(TSM,1),1)*(E(:,2))'),param.pad,2));
    J2=J2(:,1:param.spectrum_range)* ( 27539);

    %==============Power spectrum=============
    m_powSpec=real(J1).^2+real(J2).^2+imag(J1).^2+imag(J2).^2;
    m_time_deriv=-1*(real(J1).*real(J2)+imag(J1).*imag(J2));
    m_freq_deriv=((imag(J1).*real(J2)-real(J1).*imag(J2)));

m_time_deriv_max=max(m_time_deriv.^2,[],2);
m_freq_deriv_max=max(m_freq_deriv.^2,[],2);

%===

freq_winer_ampl_index=[param.min_freq_winer_ampl:param.max_freq_winer_ampl];

m_amplitude=sum(m_powSpec(:,freq_winer_ampl_index),2);

log_power=m_time_deriv(:,freq_winer_ampl_index).^2+m_freq_deriv(:,freq_winer_ampl_index).^2; 
m_SumLog=sum(log(m_powSpec(:,freq_winer_ampl_index)+eps),2);
m_LogSum=(sum(m_powSpec(:,freq_winer_ampl_index),2)); 

gravity_center=sum((ones(size(log_power,1),1)*(freq_winer_ampl_index)).*log_power,2);
gc_base=sum(log_power,2);
m_AM=sum(m_time_deriv(:,freq_winer_ampl_index),2);


gravity_center=gravity_center./max(gc_base,1)*fs/param.pad;
m_AM=m_AM./(m_amplitude+eps);
m_amplitude=log10(m_amplitude+1)*10-70; %units in Db


%===========Wiener entropy==================

m_LogSum(find(m_LogSum==0))=length(freq_winer_ampl_index); 
m_LogSum=log(m_LogSum/length(freq_winer_ampl_index)); %divide by the number of frequencies
m_Entropy=(m_SumLog/length(freq_winer_ampl_index))-m_LogSum;
m_Entropy(find(m_LogSum==0))=0; 


%============FM===================

m_FM=atan(m_time_deriv_max./(m_freq_deriv_max+eps));
%m_FM(find(m_freq_deriv_max==0))=0;

%%%%%%%%%%%%%%%%%%%%%%%%

%==========Directional Spectral derivatives=================


cFM=cos(m_FM);
sFM=sin(m_FM);

%==The image==
m_spec_deriv=m_time_deriv(:,3:255).*(sFM*ones(1,255-3+1))+m_freq_deriv(:,3:255).*(cFM*ones(1,255-3+1));
Cepstrum=(fft(m_spec_deriv./(m_powSpec(:,3:255)+eps),512,2))*( 1/2);
x=(real(Cepstrum(:,param.up_pitch:param.low_pitch))).^2+(imag(Cepstrum(:,param.up_pitch:param.low_pitch))).^2;
[m_PitchGoodness,m_Pitch]=sort(x,2);
m_PitchGoodness=m_PitchGoodness(:,end);
m_Pitch=m_Pitch(:,end);

m_Pitch(find(m_PitchGoodness<1))=1;
m_PitchGoodness=max(m_PitchGoodness,1);

m_Pitch=m_Pitch+3;
Pitch_chose= 22050./m_Pitch ; %1./(m_Pitch/1024*fs*512);

index_m_freq=find(Pitch_chose>param.pitch_HoP & (m_PitchGoodness<param.gdn_HoP | m_Entropy>param.up_wiener));

Pitch_chose(index_m_freq)=gravity_center(index_m_freq);

Pitch_weight=Pitch_chose.*m_PitchGoodness./sum(m_PitchGoodness);

m_FM=m_FM*180/pi;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function E=taper_read()
E=[0.10082636	0.488103509
0.105322801	0.500162601
0.109902844	0.512267828
0.114566609	0.52441597
0.119314209	0.536603808
0.124145724	0.548828006
0.129061237	0.561085165
0.13406077	0.573371947
0.139144346	0.585684955
0.144311979	0.598020673
0.149563625	0.610375643
0.15489924	0.622746408
0.160318747	0.635129333
0.165822059	0.647520959
0.17140907	0.659917593
0.177079603	0.672315657
0.182833508	0.684711576
0.188670605	0.697101533
0.194590658	0.709481955
0.200593442	0.721849024
0.206678703	0.734199166
0.21284613	0.746528447
0.219095424	0.75883323
0.225426242	0.7711097
0.231838226	0.783353984
0.238331005	0.795562387
0.244904146	0.807730973
0.251557231	0.819855988
0.258289784	0.831933498
0.265101343	0.843959749
0.271991402	0.855930805
0.278959394	0.867842793
0.286004782	0.879691899
0.293127	0.891474187
0.300325394	0.903185785
0.307599366	0.914822817
0.314948261	0.926381409
0.322371393	0.937857628
0.329868048	0.949247718
0.337437481	0.960547745
0.345078975	0.971753776
0.352791727	0.982862055
0.360574961	0.993868709
0.368427813	1.004769802
0.376349419	1.0155617
0.384338975	1.026240349
0.392395496	1.036802173
0.400518119	1.047243237
0.40870589	1.057559848
0.416957796	1.067748189
0.425272882	1.077804685
0.433650136	1.087725401
0.442088515	1.097506762
0.450586915	1.107145071
0.459144324	1.116636872
0.467759579	1.125978231
0.476431549	1.135165811
0.485159129	1.144196033
0.493941098	1.153065205
0.502776325	1.161769986
0.511663496	1.170306802
0.520601451	1.178672433
0.529588878	1.186863184
0.538624585	1.194875956
0.5477072	1.20270741
0.556835413	1.21035409
0.566007853	1.217812896
0.575223267	1.225080609
0.584480166	1.232154131
0.593777239	1.239030242
0.603113055	1.245705962
0.612486124	1.252178192
0.621895015	1.258444071
0.631338298	1.264500499
0.640814483	1.270344853
0.65032202	1.275974274
0.659859419	1.281385779
0.66942513	1.286576867
0.679017663	1.291544795
0.688635349	1.296286941
0.698276699	1.3008008
0.707940042	1.30508399
0.71762383	1.309133887
0.727326393	1.312948227
0.737046123	1.316524744
0.746781349	1.319861054
0.756530404	1.322955132
0.766291559	1.32580471
0.776063263	1.328407884
0.78584367	1.330762625
0.79563117	1.332866907
0.805423975	1.334718943
0.815220356	1.336316943
0.825018644	1.33765924
0.834816992	1.338744044
0.844613671	1.339569926
0.854406953	1.340135336
0.864194989	1.340438843
0.873976052	1.340479016
0.883748353	1.340254664
0.893510044	1.339764476
0.903259337	1.339007378
0.912994385	1.337982178
0.922713459	1.336688161
0.932414711	1.335124135
0.942096293	1.333289266
0.951756358	1.331183076
0.961393118	1.328804493
0.971004725	1.326153278
0.98058933	1.323228598
0.990145147	1.320030212
0.999670267	1.316557646
1.009162784	1.312810659
1.018621087	1.308789015
1.028043151	1.304492474
1.037427068	1.299921274
1.046771288	1.295075178
1.056073666	1.289954305
1.065332532	1.284559011
1.074546099	1.278889418
1.083712459	1.272946
1.092829704	1.266728997
1.101896167	1.260239124
1.110909939	1.253476977
1.119869232	1.246443033
1.128772378	1.239138246
1.13761735	1.231563449
1.146402478	1.223719358
1.155125856	1.215607166
1.163785934	1.207227945
1.172380805	1.198582649
1.180908799	1.189672828
1.18936801	1.180499554
1.197756886	1.171064377
1.206073523	1.161368728
1.214316368	1.151414156
1.222483516	1.141202331
1.230573535	1.130734921
1.238584518	1.120013714
1.246514797	1.109040618
1.254362941	1.097817659
1.262127161	1.086346745
1.269805789	1.074629903
1.277397275	1.062669516
1.28489995	1.05046773
1.292312384	1.03802681
1.299632907	1.025349379
1.30685997	1.012437582
1.313992023	0.999294221
1.321027637	0.98592186
1.32796526	0.97232312
1.334803343	0.958500803
1.341540575	0.94445771
1.348175406	0.930196762
1.354706526	0.91572094
1.361132383	0.901033223
1.367451787	0.886136711
1.373663187	0.871034622
1.379765391	0.855730176
1.385756969	0.84022665
1.391636729	0.824527323
1.39740324	0.808635712
1.403055549	0.792555213
1.408592105	0.776289403
1.414011955	0.759841859
1.419313788	0.743216217
1.424496531	0.726416171
1.429558992	0.709445536
1.434500098	0.692308068
1.439318776	0.675007641
1.444013953	0.657548249
1.448584676	0.639933705
1.453029752	0.622168124
1.457348466	0.604255557
1.461539626	0.586200118
1.465602517	0.568005979
1.469536185	0.549677253
1.473339677	0.53121829
1.477012277	0.512633264
1.48055315	0.493926555
1.483961463	0.475102544
1.487236381	0.456165582
1.490377426	0.43712011
1.493383765	0.417970628
1.496254683	0.398721576
1.498989701	0.379377574
1.501588106	0.359943092
1.504049301	0.340422779
1.506372809	0.320821226
1.508558035	0.30114311
1.51060462	0.281393051
1.512511969	0.261575758
1.514279842	0.24169597
1.515907645	0.221758381
1.517395139	0.201767772
1.518741965	0.181728885
1.519947767	0.16164653
1.521012425	0.141525477
1.521935582	0.121370547
1.522717118	0.101186559
1.523356676	0.080978341
1.523854375	0.060750734
1.524209857	0.040508576
1.524423242	0.020256713
1.52449429	-2.11E-09
1.524423242	-0.020256717
1.524209857	-0.040508579
1.523854375	-0.060750738
1.523356676	-0.080978349
1.522717118	-0.101186566
1.521935582	-0.121370554
1.521012425	-0.141525477
1.519947767	-0.16164653
1.518741965	-0.181728899
1.517395139	-0.201767772
1.515907645	-0.221758395
1.514279842	-0.24169597
1.512511969	-0.261575758
1.51060462	-0.281393051
1.508558035	-0.30114311
1.506372809	-0.320821226
1.504049301	-0.340422779
1.501588106	-0.359943092
1.498989701	-0.379377574
1.496254683	-0.398721606
1.493383765	-0.417970628
1.490377426	-0.43712011
1.487236381	-0.456165582
1.483961463	-0.475102544
1.48055315	-0.493926555
1.477012277	-0.512633264
1.473339677	-0.53121829
1.469536185	-0.549677253
1.465602517	-0.568005979
1.461539626	-0.586200118
1.457348466	-0.604255617
1.453029752	-0.622168124
1.448584676	-0.639933705
1.444013953	-0.657548249
1.439318776	-0.675007701
1.434500098	-0.692308068
1.429558992	-0.709445536
1.424496531	-0.726416171
1.419313788	-0.743216217
1.414011955	-0.759841859
1.408592105	-0.776289403
1.403055549	-0.792555213
1.39740324	-0.808635712
1.391636729	-0.824527323
1.385756969	-0.84022665
1.379765391	-0.855730176
1.373663187	-0.871034682
1.367451787	-0.88613677
1.361132383	-0.901033223
1.354706526	-0.91572094
1.348175406	-0.930196762
1.341540575	-0.94445771
1.334803343	-0.958500803
1.32796526	-0.97232312
1.321027637	-0.98592186
1.313992023	-0.999294281
1.30685997	-1.012437582
1.299632907	-1.025349379
1.292312384	-1.03802681
1.28489995	-1.05046773
1.277397275	-1.062669516
1.269805789	-1.074629903
1.262127161	-1.086346745
1.254362941	-1.097817659
1.246514797	-1.109040618
1.238584518	-1.120013714
1.230573535	-1.130734921
1.222483516	-1.141202331
1.214316368	-1.151414156
1.206073523	-1.161368728
1.197756886	-1.171064377
1.18936801	-1.180499554
1.180908799	-1.189672828
1.172380805	-1.198582649
1.163785934	-1.207227945
1.155125856	-1.215607166
1.146402478	-1.223719358
1.13761735	-1.231563449
1.128772378	-1.239138246
1.119869232	-1.246443033
1.110909939	-1.253476977
1.101896167	-1.260239124
1.092829704	-1.266728997
1.083712459	-1.272946
1.074546099	-1.278889418
1.065332532	-1.284559011
1.056073666	-1.289954305
1.046771288	-1.295075178
1.037427068	-1.299921274
1.028043151	-1.304492474
1.018621087	-1.308789015
1.009162784	-1.312810659
0.999670208	-1.316557646
0.990145147	-1.320030212
0.98058933	-1.323228598
0.971004725	-1.326153278
0.961393118	-1.328804493
0.951756358	-1.331183076
0.942096293	-1.333289266
0.932414711	-1.335124135
0.922713459	-1.336688161
0.912994385	-1.337982178
0.903259337	-1.339007378
0.893510044	-1.339764476
0.883748353	-1.340254664
0.873976052	-1.340479016
0.864194989	-1.340438843
0.854406953	-1.340135336
0.844613671	-1.339569926
0.834816992	-1.338744044
0.825018644	-1.33765924
0.815220356	-1.336316943
0.805423975	-1.334718943
0.79563117	-1.332866907
0.78584367	-1.330762625
0.776063263	-1.328407884
0.766291559	-1.32580471
0.756530404	-1.322955132
0.746781349	-1.319861054
0.737046123	-1.316524744
0.727326393	-1.312948227
0.71762383	-1.309133887
0.707940042	-1.30508399
0.698276699	-1.3008008
0.688635349	-1.296286941
0.679017663	-1.291544795
0.66942513	-1.286576867
0.659859419	-1.281385779
0.65032202	-1.275974274
0.640814483	-1.270344853
0.631338298	-1.264500499
0.621895015	-1.258444071
0.612486124	-1.252178192
0.603113055	-1.245705962
0.593777239	-1.239030242
0.584480166	-1.232154131
0.575223267	-1.225080609
0.566007853	-1.217812896
0.556835413	-1.21035409
0.5477072	-1.20270741
0.538624585	-1.194875956
0.529588878	-1.186863184
0.520601451	-1.178672433
0.511663496	-1.170306802
0.502776325	-1.161769986
0.493941098	-1.153065205
0.485159129	-1.144196033
0.476431549	-1.135165811
0.467759579	-1.125978231
0.459144324	-1.116636872
0.450586915	-1.107145071
0.442088515	-1.097506762
0.433650136	-1.087725401
0.425272882	-1.077804685
0.416957796	-1.067748189
0.40870586	-1.057559848
0.400518119	-1.047243237
0.392395496	-1.036802173
0.384338945	-1.026240349
0.376349419	-1.0155617
0.368427813	-1.004769802
0.360574961	-0.993868709
0.352791727	-0.982862055
0.345078975	-0.971753776
0.337437481	-0.960547745
0.329868048	-0.949247718
0.322371393	-0.937857628
0.314948261	-0.926381409
0.307599366	-0.914822817
0.300325394	-0.903185785
0.293126971	-0.891474187
0.286004782	-0.879691899
0.278959394	-0.867842793
0.271991402	-0.855930805
0.265101343	-0.843959749
0.258289784	-0.831933498
0.251557231	-0.819855988
0.244904146	-0.807730973
0.238331005	-0.795562387
0.231838226	-0.783353984
0.225426242	-0.7711097
0.219095424	-0.75883323
0.21284613	-0.746528447
0.206678703	-0.734199166
0.200593442	-0.721849024
0.194590658	-0.709481955
0.188670605	-0.697101533
0.182833508	-0.684711576
0.177079603	-0.672315657
0.171409056	-0.659917593
0.165822059	-0.647520959
0.160318747	-0.635129333
0.15489924	-0.622746408
0.149563611	-0.610375643
0.144311965	-0.598020673
0.139144346	-0.585684955
0.134060755	-0.573371947
0.129061222	-0.561085165
0.124145724	-0.548828006
0.119314201	-0.536603808
0.114566602	-0.52441597
0.109902844	-0.512267828
0.105322801	-0.500162601
0.100826353	-0.488103509];
