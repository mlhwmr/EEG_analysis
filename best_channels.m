% Dla zadanej kombinacji kana³ów wylicza parametry stabilnoœci dla epochów
% z R (z baselinem i bez) oraz z S (z baselinem i bez)
% kombinacje kana³ów trzeba si³¹ rzeczy wpisaæ samemu o ile nie dobieramy
% ich na razie losowo, np.channels = {'Fp1', 'Fp2'}

%Dostêpna pula kana³ów
% C3, C4, CP1, CP2, CP5, CP6, F3, F4, F7, F8, Fz, FC1, FC2, FC5, FC6, Fp1,
% Fp2, FT10, FT9, O1, O2, Oz, P3, P4, P7, P8, Pz,T7, T8, TP10, TP9

%-------------------------------------------------------------------------
% Tylko to zmieniamy podaj¹c wybrane kane³y
%-------------------------------------------------------------------------
% Przyk³ad
%{
channels = "temporal_mix"
chans =["FT10", "FT9", "T7", "T8", "TP10", "TP9"];
chans_cell = {'FT10', 'FT9', 'T7', 'T8', 'TP10', 'TP9'}
%}
%---------------------------R_basline-------------------------------------
function best_channels(channels,chans,chans_cell)
    %eeglab 

    %EEG = eeg_checkset( EEG );
    %EEG.etc.eeglabvers = '2019.1'; % this tracks which version of EEGLAB is being used, you may ignore it
    EEG = pop_loadset('filename','sb23_1-1245Hz.set','filepath','C:\\Users\\Asus\\Desktop\\licencjat\\zadanie2_lic\\sb23_flanker\\');
    EEG = eeg_checkset( EEG );
    EEG = pop_select( EEG, 'channel',{'C3' 'C4' 'CP1' 'CP2' 'CP5' 'CP6' 'F3' 'F4' 'F7' 'F8' 'Fz' 'FC1' 'FC2' 'FC5' 'FC6' 'Fp1' 'Fp2' 'FT10' 'FT9' 'O1' 'O2' 'Oz' 'P3' 'P4' 'P7' 'P8' 'Pz' 'T7' 'T8' 'TP10' 'TP9'});
    EEG = pop_loadset('filename','flanker_31_baseline.set','filepath','C:\\Users\\Asus\\Desktop\\licencjat\\zadanie2_lic\\sb23_flanker\\modified_flanker\\');
    EEG = eeg_checkset( EEG );
    EEG = pop_select( EEG, 'channel',chans_cell);
    EEG = eeg_checkset( EEG );
    EEG = pop_epoch( EEG, {  'R  1'  'R  2'  'R  8'  }, [-0.2         0.5], 'newname', 'sb_fl_raw epochs', 'epochinfo', 'yes');
    EEG = eeg_checkset( EEG );
    EEG = pop_rmbase( EEG, [-200 0] ,[]);
    EEG = eeg_checkset( EEG );

    stabil_params = {};
    As = {};
    len = 384;
    ch = length(chans);

    for i = 1:len
        [w,A,C,SBC,FPE,th]=arfit(EEG.data(1:ch,:,i)',1,1); % transponowaæ
        [S,Serr,per,tau,exctn]=armode(A,C,th); %- ok
        stabil_params{i} = abs(eig(S));
        As{i} = A;
    end

    nameA = "A_" + channels;
    nameSP = "SP_" + channels;
    fileA = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\with_baseline\R_kombinacje_rozne\A\" + nameA;
    fileSP = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\with_baseline\R_kombinacje_rozne\SP\" + nameSP;
    save(fileA, 'As')
    save(fileSP, 'stabil_params')

    %------------------------------R_no_baseline----------------------------------

    EEG = pop_rmbase( EEG, [-200 0] ,[]);
    EEG = eeg_checkset( EEG );

    stabil_params = {};
    As = {};

    for i = 1:len
        [w,A,C,SBC,FPE,th]=arfit(EEG.data(1:ch,:,i)',1,1); % transponowaæ
        [S,Serr,per,tau,exctn]=armode(A,C,th); %- ok
        stabil_params{i} = abs(eig(S));
        As{i} = A;
    end

    nameA = "A_" + channels;
    nameSP = "SP_" + channels;
    fileA = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\no_baseline\R_kombinacje_rozne\A\"+ nameA;
    fileSP = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\no_baseline\R_kombinacje_rozne\SP\"+ nameSP;
    save(fileA, 'As')
    save(fileSP, 'stabil_params')

    %-------------------------R1,R2,R8_baseline&no_baseline--------------------

    %---------------------------------R1_baseline------------------------------
    EEG_R1 = pop_selectevent( EEG, 'type',{'R  1'},'deleteevents','off','deleteepochs','on','invertepochs','off');
    EEG = eeg_checkset( EEG );

    stabil_params = {};
    As = {};
    size_n = size(EEG_R1.data);
    len = size_n(3);
    for i = 1:len 
        [w,A,C,SBC,FPE,th]=arfit(EEG_R1.data(1:ch,:,i)',1,1); % transponowaæ
        [S,Serr,per,tau,exctn]=armode(A,C,th); %- ok
        stabil_params{i} = abs(eig(S));
        As{i} = A;
    end

    nameA = "A_" + channels;
    nameSP = "SP_" + channels;
    fileA = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\with_baseline\R1_kombinacje_rozne\A\" + nameA;
    fileSP = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\with_baseline\R1_kombinacje_rozne\SP\" + nameSP;
    save(fileA, 'As')
    save(fileSP, 'stabil_params')

    %--------------------------------R1_no_baseline----------------------------
    EEG_R1 = pop_rmbase( EEG_R1, [-200 0] ,[]);
    EEG = eeg_checkset( EEG );

    stabil_params = {};
    As = {};
    size_n = size(EEG_R1.data);
    len = size_n(3);
    for i = 1:len 
        [w,A,C,SBC,FPE,th]=arfit(EEG_R1.data(1:ch,:,i)',1,1); % transponowaæ
        [S,Serr,per,tau,exctn]=armode(A,C,th); %- ok
        stabil_params{i} = abs(eig(S));
        As{i} = A;
    end

    nameA = "A_" + channels;
    nameSP = "SP_" + channels;
    fileA = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\no_baseline\R1_kombinacje_rozne\A\" + nameA;
    fileSP = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\no_baseline\R1_kombinacje_rozne\SP\" + nameSP;
    save(fileA, 'As')
    save(fileSP, 'stabil_params')

    %----------------------------------R2_baseline-----------------------------

    EEG_R2 = pop_selectevent( EEG, 'type',{'R  2'},'deleteevents','off','deleteepochs','on','invertepochs','off');
    EEG = eeg_checkset( EEG );

    size_R2 = size(EEG_R2.data);
    len_R2 = size_R2(3);

    stabil_params = {};
    As = {};
    size_n = size(EEG_R2.data);
    len = size_n(3);
    for i = 1:len 
        [w,A,C,SBC,FPE,th]=arfit(EEG_R2.data(1:ch,:,i)',1,1); % transponowaæ
        [S,Serr,per,tau,exctn]=armode(A,C,th); %- ok
        stabil_params{i} = abs(eig(S));
        As{i} = A;
    end

    nameA = "A_" + channels;
    nameSP = "SP_" + channels;
    fileA = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\with_baseline\R2_kombinacje_rozne\A\" + nameA;
    fileSP = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\with_baseline\R2_kombinacje_rozne\SP\" + nameSP;
    save(fileA, 'As')
    save(fileSP, 'stabil_params')

    %--------------------------------R2_no_baseline----------------------------
    EEG_R2 = pop_rmbase( EEG_R2, [-200 0] ,[]);
    EEG = eeg_checkset( EEG );

    stabil_params = {};
    As = {};
    size_n = size(EEG_R2.data);
    len = size_n(3);
    for i = 1:len 
        [w,A,C,SBC,FPE,th]=arfit(EEG_R2.data(1:ch,:,i)',1,1); % transponowaæ
        [S,Serr,per,tau,exctn]=armode(A,C,th); %- ok
        stabil_params{i} = abs(eig(S));
        As{i} = A;
    end

    nameA = "A_" + channels;
    nameSP = "SP_" + channels;
    fileA = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\no_baseline\R2_kombinacje_rozne\A\" + nameA;
    fileSP = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\no_baseline\R2_kombinacje_rozne\SP\" + nameSP;
    save(fileA, 'As')
    save(fileSP, 'stabil_params')


    %--------------------------------R8_baseline-------------------------------
    EEG_R8 = pop_selectevent( EEG, 'type',{'R  8'},'deleteevents','off','deleteepochs','on','invertepochs','off');
    EEG = eeg_checkset( EEG );

    size_R8 = size(EEG_R8.data);
    len_R8 = size_R8(3);

    stabil_params = {};
    As = {};
    size_n = size(EEG_R8.data);
    len = size_n(3);
    for i = 1:len 
        [w,A,C,SBC,FPE,th]=arfit(EEG_R8.data(1:ch,:,i)',1,1); % transponowaæ
        [S,Serr,per,tau,exctn]=armode(A,C,th); %- ok
        stabil_params{i} = abs(eig(S));
        As{i} = A;
    end

    nameA = "A_" + channels;
    nameSP = "SP_" + channels;
    fileA = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\with_baseline\R8_kombinacje_rozne\A\" + nameA;
    fileSP = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\with_baseline\R8_kombinacje_rozne\SP\" + nameSP;
    save(fileA, 'As')
    save(fileSP, 'stabil_params')

    %--------------------------------R8_no_baseline----------------------------
    EEG_R8 = pop_rmbase( EEG_R8, [-200 0] ,[]);
    EEG = eeg_checkset( EEG );

    stabil_params = {};
    As = {};
    size_n = size(EEG_R8.data);
    len = size_n(3);
    for i = 1:len 
        [w,A,C,SBC,FPE,th]=arfit(EEG_R8.data(1:ch,:,i)',1,1); % transponowaæ
        [S,Serr,per,tau,exctn]=armode(A,C,th); %- ok
        stabil_params{i} = abs(eig(S));
        As{i} = A;
    end

    nameA = "A_" + channels;
    nameSP = "SP_" + channels;
    fileA = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\no_baseline\R8_kombinacje_rozne\A\" + nameA;
    fileSP = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\no_baseline\R8_kombinacje_rozne\SP\" + nameSP;
    save(fileA, 'As')
    save(fileSP, 'stabil_params')


    %------------------------------S_baseline----------------------------------
    EEG = eeg_checkset( EEG );
    EEG.etc.eeglabvers = '2019.1'; % this tracks which version of EEGLAB is being used, you may ignore it
    EEG = pop_loadset('filename','sb23_1-1245Hz.set','filepath','C:\\Users\\Asus\\Desktop\\licencjat\\zadanie2_lic\\sb23_flanker\\');
    EEG = eeg_checkset( EEG );
    EEG = pop_select( EEG, 'channel',{'C3' 'C4' 'CP1' 'CP2' 'CP5' 'CP6' 'F3' 'F4' 'F7' 'F8' 'Fz' 'FC1' 'FC2' 'FC5' 'FC6' 'Fp1' 'Fp2' 'FT10' 'FT9' 'O1' 'O2' 'Oz' 'P3' 'P4' 'P7' 'P8' 'Pz' 'T7' 'T8' 'TP10' 'TP9'});
    EEG = pop_loadset('filename','flanker_31_baseline.set','filepath','C:\\Users\\Asus\\Desktop\\licencjat\\zadanie2_lic\\sb23_flanker\\modified_flanker\\');
    EEG = eeg_checkset( EEG );
    EEG = pop_select( EEG, 'channel',chans_cell);
    EEG = eeg_checkset( EEG );
    EEG = pop_epoch( EEG, {  'S  3'  'S  4'  'S  5'  'S  6'  }, [-0.2         0.5], 'newname', 'sb_fl_raw epochs', 'epochinfo', 'yes');
    EEG = eeg_checkset( EEG );

    stabil_params = {};
    As = {};

    size_n = size(EEG.data);
    len = size_n(3);

    for i = 1:len
        [w,A,C,SBC,FPE,th]=arfit(EEG.data(1:ch,:,i)',1,1); % transponowaæ
        [S,Serr,per,tau,exctn]=armode(A,C,th); %- ok
        stabil_params{i} = abs(eig(S));
        As{i} = A;
    end

    nameA = "A_" + channels;
    nameSP = "SP_" + channels;
    fileA = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\with_baseline\S_kombinacje_rozne\A\"+ nameA;
    fileSP = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\with_baseline\S_kombinacje_rozne\SP\"+ nameSP;
    save(fileA, 'As')
    save(fileSP, 'stabil_params')
    %--------------------------------S_no baseline-----------------------------

    EEG = pop_rmbase( EEG, [-200 0] ,[]);
    EEG = eeg_checkset( EEG );

    stabil_params = {};
    As = {};

    for i = 1:len
        [w,A,C,SBC,FPE,th]=arfit(EEG.data(1:ch,:,i)',1,1); % transponowaæ
        [S,Serr,per,tau,exctn]=armode(A,C,th); %- ok
        stabil_params{i} = abs(eig(S));
        As{i} = A;
    end

    nameA = "A_" + channels;
    nameSP = "SP_" + channels;
    fileA = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\no_baseline\S_kombinacje_rozne\A\"+ nameA;
    fileSP = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\no_baseline\S_kombinacje_rozne\SP\"+ nameSP;
    save(fileA, 'As')
    save(fileSP, 'stabil_params')
    %--------------------------------------------------------------------------
    %--------------------S3,S4,S5,S6_baseline&no_baseline----------------------

    %---------------------------------S3_baseline------------------------------
    EEG_S3 = pop_selectevent( EEG, 'type',{'S  3'},'deleteevents','off','deleteepochs','on','invertepochs','off');
    EEG = eeg_checkset( EEG );

    size_S3456 = size(EEG_S3.data);
    len_S3456 = size_S3456(3);

    stabil_params = {};
    As = {};
    size_n = size(EEG_S3.data);
    len = size_n(3);

    for i = 1:len 
        [w,A,C,SBC,FPE,th]=arfit(EEG_S3.data(1:ch,:,i)',1,1); % transponowaæ
        [S,Serr,per,tau,exctn]=armode(A,C,th); %- tu siê wywala
        stabil_params{i} = abs(eig(S));
        As{i} = A;
    end

    nameA = "A_" + channels;
    nameSP = "SP_" + channels;
    fileA = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\with_baseline\S3_kombinacje_rozne\A\" + nameA;
    fileSP = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\with_baseline\S3_kombinacje_rozne\SP\" + nameSP;
    save(fileA, 'As')
    save(fileSP, 'stabil_params')

    %--------------------------------S3_no_baseline----------------------------
    EEG_S3 = pop_rmbase( EEG_S3, [-200 0] ,[]);
    EEG = eeg_checkset( EEG );

    stabil_params = {};
    As = {};
    size_n = size(EEG_S3.data);
    len = size_n(3);
    for i = 1:len 
        [w,A,C,SBC,FPE,th]=arfit(EEG_S3.data(1:ch,:,i)',1,1); % transponowaæ
        [S,Serr,per,tau,exctn]=armode(A,C,th); %- tu siê wywala
        stabil_params{i} = abs(eig(S));
        As{i} = A;
    end

    nameA = "A_" + channels;
    nameSP = "SP_" + channels;
    fileA = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\no_baseline\S3_kombinacje_rozne\A\" + nameA;
    fileSP = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\no_baseline\S3_kombinacje_rozne\SP\" + nameSP;
    save(fileA, 'As')
    save(fileSP, 'stabil_params')

    %----------------------------------S4_baseline-----------------------------

    EEG_S4 = pop_selectevent( EEG, 'type',{'S  4'},'deleteevents','off','deleteepochs','on','invertepochs','off');
    EEG = eeg_checkset( EEG );

    stabil_params = {};
    As = {};
    size_n = size(EEG_S4.data);
    len = size_n(3);
    for i = 1:len 
        [w,A,C,SBC,FPE,th]=arfit(EEG_S4.data(1:ch,:,i)',1,1); % transponowaæ
        [S,Serr,per,tau,exctn]=armode(A,C,th); %- tu siê wywala
        stabil_params{i} = abs(eig(S));
        As{i} = A;
    end

    nameA = "A_" + channels;
    nameSP = "SP_" + channels;
    fileA = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\with_baseline\S4_kombinacje_rozne\A\" + nameA;
    fileSP = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\with_baseline\S4_kombinacje_rozne\SP\" + nameSP;
    save(fileA, 'As')
    save(fileSP, 'stabil_params')

    %--------------------------------S4_no_baseline----------------------------
    EEG_S4 = pop_rmbase( EEG_S4, [-200 0] ,[]);
    EEG = eeg_checkset( EEG );

    stabil_params = {};
    As = {};
    size_n = size(EEG_S4.data);
    len = size_n(3);
    for i = 1:len 
        [w,A,C,SBC,FPE,th]=arfit(EEG_S4.data(1:ch,:,i)',1,1); % transponowaæ
        [S,Serr,per,tau,exctn]=armode(A,C,th); %- tu siê wywala
        stabil_params{i} = abs(eig(S));
        As{i} = A;
    end

    nameA = "A_" + channels;
    nameSP = "SP_" + channels;
    fileA = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\no_baseline\S4_kombinacje_rozne\A\" + nameA;
    fileSP = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\no_baseline\S4_kombinacje_rozne\SP\" + nameSP;
    save(fileA, 'As')
    save(fileSP, 'stabil_params')


    %--------------------------------S5_baseline-------------------------------
    EEG_S5 = pop_selectevent( EEG, 'type',{'S  5'},'deleteevents','off','deleteepochs','on','invertepochs','off');
    EEG = eeg_checkset( EEG );

    stabil_params = {};
    As = {};
    size_n = size(EEG_S5.data);
    len = size_n(3);
    for i = 1:len 
        [w,A,C,SBC,FPE,th]=arfit(EEG_S5.data(1:ch,:,i)',1,1); % transponowaæ
        [S,Serr,per,tau,exctn]=armode(A,C,th); %- tu siê wywala
        stabil_params{i} = abs(eig(S));
        As{i} = A;
    end

    nameA = "A_" + channels;
    nameSP = "SP_" + channels;
    fileA = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\with_baseline\S5_kombinacje_rozne\A\" + nameA;
    fileSP = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\with_baseline\S5_kombinacje_rozne\SP\" + nameSP;
    save(fileA, 'As')
    save(fileSP, 'stabil_params')

    %--------------------------------S5_no_baseline----------------------------
    EEG_S5 = pop_rmbase( EEG_S5, [-200 0] ,[]);
    EEG = eeg_checkset( EEG );

    stabil_params = {};
    As = {};
    size_n = size(EEG_S5.data);
    len = size_n(3);
    for i = 1:len 
        [w,A,C,SBC,FPE,th]=arfit(EEG_S5.data(1:ch,:,i)',1,1); % transponowaæ
        [S,Serr,per,tau,exctn]=armode(A,C,th); %- tu siê wywala
        stabil_params{i} = abs(eig(S));
        As{i} = A;
    end

    nameA = "A_" + channels;
    nameSP = "SP_" + channels;
    fileA = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\no_baseline\S5_kombinacje_rozne\A\" + nameA;
    fileSP = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\no_baseline\S5_kombinacje_rozne\SP\" + nameSP;
    save(fileA, 'As')
    save(fileSP, 'stabil_params')


    %--------------------------------S6_baseline-------------------------------
    EEG_S6 = pop_selectevent( EEG, 'type',{'S  6'},'deleteevents','off','deleteepochs','on','invertepochs','off');
    EEG = eeg_checkset( EEG );

    stabil_params = {};
    As = {};
    size_n = size(EEG_S6.data);
    len = size_n(3);
    for i = 1:len 
        [w,A,C,SBC,FPE,th]=arfit(EEG_S6.data(1:ch,:,i)',1,1); % transponowaæ
        [S,Serr,per,tau,exctn]=armode(A,C,th); %- tu siê wywala
        stabil_params{i} = abs(eig(S));
        As{i} = A;
    end

    nameA = "A_" + channels;
    nameSP = "SP_" + channels;
    fileA = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\with_baseline\S6_kombinacje_rozne\A\" + nameA;
    fileSP = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\with_baseline\S6_kombinacje_rozne\SP\" + nameSP;
    save(fileA, 'As')
    save(fileSP, 'stabil_params')

    %--------------------------------S6_no_baseline----------------------------
    EEG_S6 = pop_rmbase( EEG_S6, [-200 0] ,[]);
    EEG = eeg_checkset( EEG );

    stabil_params = {};
    As = {};
    size_n = size(EEG_S6.data);
    len = size_n(3);
    for i = 1:len 
        [w,A,C,SBC,FPE,th]=arfit(EEG_S6.data(1:ch,:,i)',1,1); % transponowaæ
        [S,Serr,per,tau,exctn]=armode(A,C,th); %- tu siê wywala
        stabil_params{i} = abs(eig(S));
        As{i} = A;
    end

    nameA = "A_" + channels;
    nameSP = "SP_" + channels;
    fileA = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\no_baseline\S6_kombinacje_rozne\A\" + nameA;
    fileSP = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\no_baseline\S6_kombinacje_rozne\SP\" + nameSP;
    save(fileA, 'As')
    save(fileSP, 'stabil_params')

    %============================REO_REC=======================================
    k = ["AF3","AF4","AF7","AF8","AFz","C1","C2","C3","C4","C5","C6","CP1","CP2","CP3","CP4","CP5","CP6","CPz","F1","F2","F3","F4","F5","F6","F7","F8","Fz","FC1","FC2","FC3","FC4","FC5","FC6","FCz","Fp1","Fp2","FT10","FT7","FT8","FT9","O1","O2","Oz","P1","P2","P3","P4","P5","P6","P7","P8","Pz","PO3","PO4","PO7","PO8","POz","T7","T8","TP10","TP7","TP8","TP9"];         

    ch_ind = zeros(1,length(chans));

    for i=1:length(chans)
        ch_ind(i) = find(k==chans(i));
    end
    %--------------------------------reo_baseline------------------------------
    EEG = pop_loadset('filename','sb23_1-1245Hz.set','filepath','C:\\Users\\Asus\\Desktop\\licencjat\\zadanie2_lic\\sb23_reo\\');
    EEG = eeg_checkset( EEG );
    EEG = eeg_regepochs(EEG, 0.7, [-0.2 0.5], NaN)
    EEG = eeg_checkset( EEG );

    stabil_params = {};
    As = {};

    size_n = size(EEG.data);
    len = size_n(3);

    for i = 1:len 
        [w,A,C,SBC,FPE,th]=arfit(EEG.data(ch_ind,:,i)',1,1); % transponowaæ
        [S,Serr,per,tau,exctn]=armode(A,C,th); %- ok
        stabil_params{i} = abs(eig(S));
        As{i} = A;
    end

    nameA = "A_" + channels;
    nameSP = "SP_" + channels;
    fileA = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\with_baseline\reo_kombinacje_rozne\A\" + nameA;
    fileSP = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\with_baseline\reo_kombinacje_rozne\SP\" + nameSP;
    save(fileA, 'As')
    save(fileSP, 'stabil_params')
    %----------------------------reoR1_baseline--------------------------------
    stabil_params = {};
    As = {};

    size_n = size(EEG.data);
    len = size_n(3);

    for i = 1:len 
        [w,A,C,SBC,FPE,th]=arfit(EEG.data(ch_ind,:,i)',1,1); % transponowaæ
        [S,Serr,per,tau,exctn]=armode(A,C,th); %- ok
        stabil_params{i} = abs(eig(S));
        As{i} = A;
    end

    nameA = "A_" + channels;
    nameSP = "SP_" + channels;
    fileA = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\with_baseline\reoR1_kombinacje_rozne\A\" + nameA;
    fileSP = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\with_baseline\reoR1_kombinacje_rozne\SP\" + nameSP;
    save(fileA, 'As')
    save(fileSP, 'stabil_params')
    %--------------------------------------------------------------------------
    %----------------------------reoR2_baseline--------------------------------
    stabil_params = {};
    As = {};

    len = len_R2;

    for i = 50:(len+50) 
        [w,A,C,SBC,FPE,th]=arfit(EEG.data(ch_ind,:,i)',1,1); % transponowaæ
        [S,Serr,per,tau,exctn]=armode(A,C,th); %- ok
        stabil_params{i} = abs(eig(S));
        As{i} = A;
    end

    nameA = "A_" + channels;
    nameSP = "SP_" + channels;
    fileA = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\with_baseline\reoR2_kombinacje_rozne\A\" + nameA;
    fileSP = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\with_baseline\reoR2_kombinacje_rozne\SP\" + nameSP;
    save(fileA, 'As')
    save(fileSP, 'stabil_params')
    %--------------------------------------------------------------------------
    %----------------------------reoR8_baseline--------------------------------
    stabil_params = {};
    As = {};

    len = len_R8;

    for i = 50:(len+50) 
        [w,A,C,SBC,FPE,th]=arfit(EEG.data(ch_ind,:,i)',1,1); % transponowaæ
        [S,Serr,per,tau,exctn]=armode(A,C,th); %- ok
        stabil_params{i} = abs(eig(S));
        As{i} = A;
    end

    nameA = "A_" + channels;
    nameSP = "SP_" + channels;
    fileA = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\with_baseline\reoR8_kombinacje_rozne\A\" + nameA;
    fileSP = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\with_baseline\reoR8_kombinacje_rozne\SP\" + nameSP;
    save(fileA, 'As')
    save(fileSP, 'stabil_params')
    %--------------------------------------------------------------------------
    %----------------------------reoS3456_baseline--------------------------------
    stabil_params = {};
    As = {};

    len = len_S3456;

    for i = 50:(len+50) 
        [w,A,C,SBC,FPE,th]=arfit(EEG.data(ch_ind,:,i)',1,1); % transponowaæ
        [S,Serr,per,tau,exctn]=armode(A,C,th); %- ok
        stabil_params{i} = abs(eig(S));
        As{i} = A;
    end

    nameA = "A_" + channels;
    nameSP = "SP_" + channels;
    fileA = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\with_baseline\reoS3456_kombinacje_rozne\A\" + nameA;
    fileSP = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\with_baseline\reoS3456_kombinacje_rozne\SP\" + nameSP;
    save(fileA, 'As')
    save(fileSP, 'stabil_params')
    %--------------------------------------------------------------------------
    %------------------------------reo_no_baseline-----------------------------
    EEG = pop_rmbase(EEG, [-200 0]);
    EEG = eeg_checkset( EEG );

    stabil_params = {};
    As = {};

    for i = 1:len 
        [w,A,C,SBC,FPE,th]=arfit(EEG.data(ch_ind,:,i)',1,1); % transponowaæ
        [S,Serr,per,tau,exctn]=armode(A,C,th); %- ok
        stabil_params{i} = abs(eig(S));
        As{i} = A;
    end

    nameA = "A_" + channels;
    nameSP = "SP_" + channels;
    fileA = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\no_baseline\reo_kombinacje_rozne\A\" + nameA;
    fileSP = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\no_baseline\reo_kombinacje_rozne\SP\" + nameSP;
    save(fileA, 'As');
    save(fileSP, 'stabil_params');
    %--------------------------------------------------------------------------
    %----------------------------reoR1_no_baseline--------------------------------
    stabil_params = {};
    As = {};

    size_n = size(EEG.data);
    len = size_n(3);

    for i = 1:len 
        [w,A,C,SBC,FPE,th]=arfit(EEG.data(ch_ind,:,i)',1,1); % transponowaæ
        [S,Serr,per,tau,exctn]=armode(A,C,th); %- ok
        stabil_params{i} = abs(eig(S));
        As{i} = A;
    end

    nameA = "A_" + channels;
    nameSP = "SP_" + channels;
    fileA = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\no_baseline\reoR1_kombinacje_rozne\A\" + nameA;
    fileSP = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\no_baseline\reoR1_kombinacje_rozne\SP\" + nameSP;
    save(fileA, 'As')
    save(fileSP, 'stabil_params')
    %--------------------------------------------------------------------------
    %----------------------------reoR2_no_baseline--------------------------------
    stabil_params = {};
    As = {};

    len = len_R2;

    for i = 50:(len+50) 
        [w,A,C,SBC,FPE,th]=arfit(EEG.data(ch_ind,:,i)',1,1); % transponowaæ
        [S,Serr,per,tau,exctn]=armode(A,C,th); %- ok
        stabil_params{i} = abs(eig(S));
        As{i} = A;
    end

    nameA = "A_" + channels;
    nameSP = "SP_" + channels;
    fileA = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\no_baseline\reoR2_kombinacje_rozne\A\" + nameA;
    fileSP = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\no_baseline\reoR2_kombinacje_rozne\SP\" + nameSP;
    save(fileA, 'As')
    save(fileSP, 'stabil_params')
    %--------------------------------------------------------------------------
    %----------------------------reoR8_no_baseline--------------------------------
    stabil_params = {};
    As = {};

    len = len_R8;

    for i = 50:(len+50) 
        [w,A,C,SBC,FPE,th]=arfit(EEG.data(ch_ind,:,i)',1,1); % transponowaæ
        [S,Serr,per,tau,exctn]=armode(A,C,th); %- ok
        stabil_params{i} = abs(eig(S));
        As{i} = A;
    end

    nameA = "A_" + channels;
    nameSP = "SP_" + channels;
    fileA = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\no_baseline\reoR8_kombinacje_rozne\A\" + nameA;
    fileSP = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\no_baseline\reoR8_kombinacje_rozne\SP\" + nameSP;
    save(fileA, 'As')
    save(fileSP, 'stabil_params')
    %--------------------------------------------------------------------------
    %----------------------------reoS3456_no_baseline--------------------------------
    stabil_params = {};
    As = {};

    len = len_S3456;

    for i = 50:(len+50) 
        [w,A,C,SBC,FPE,th]=arfit(EEG.data(ch_ind,:,i)',1,1); % transponowaæ
        [S,Serr,per,tau,exctn]=armode(A,C,th); %- ok
        stabil_params{i} = abs(eig(S));
        As{i} = A;
    end

    nameA = "A_" + channels;
    nameSP = "SP_" + channels;
    fileA = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\no_baseline\reoS3456_kombinacje_rozne\A\" + nameA;
    fileSP = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\no_baseline\reoS3456_kombinacje_rozne\SP\" + nameSP;
    save(fileA, 'As')
    save(fileSP, 'stabil_params')
    %--------------------------------------------------------------------------
    %--------------------------------rec_baseline------------------------------
    EEG = pop_loadset('filename','sb23_1-1245Hz.set','filepath','C:\\Users\\Asus\\Desktop\\licencjat\\zadanie2_lic\\sb23_rec\\');
    EEG = eeg_checkset( EEG );
    EEG = eeg_regepochs(EEG, 0.7, [-0.2 0.5], NaN)
    EEG = eeg_checkset( EEG );
    %EEG = pop_rmbase(EEG, [-200 0]);
    %EEG = eeg_checkset( EEG );

    size_n = size(EEG.data);
    len = size_n(3);


    stabil_params = {};
    As = {};

    for i = 1:len 
        [w,A,C,SBC,FPE,th]=arfit(EEG.data(ch_ind,:,i)',1,1); % transponowaæ
        [S,Serr,per,tau,exctn]=armode(A,C,th); %- ok
        stabil_params{i} = abs(eig(S));
        As{i} = A;
    end

    nameA = "A_" + channels;
    nameSP = "SP_" + channels;
    fileA = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\with_baseline\rec_kombinacje_rozne\A\" + nameA;
    fileSP = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\with_baseline\rec_kombinacje_rozne\SP\" + nameSP;
    save(fileA, 'As')
    save(fileSP, 'stabil_params')

    %-------------------------------------------------------------------------

    %----------------------------recR1_baseline--------------------------------
    stabil_params = {};
    As = {};

    size_n = size(EEG.data);
    len = size_n(3);

    for i = 1:len 
        [w,A,C,SBC,FPE,th]=arfit(EEG.data(ch_ind,:,i)',1,1); % transponowaæ
        [S,Serr,per,tau,exctn]=armode(A,C,th); %- ok
        stabil_params{i} = abs(eig(S));
        As{i} = A;
    end

    nameA = "A_" + channels;
    nameSP = "SP_" + channels;
    fileA = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\with_baseline\recR1_kombinacje_rozne\A\" + nameA;
    fileSP = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\with_baseline\recR1_kombinacje_rozne\SP\" + nameSP;
    save(fileA, 'As')
    save(fileSP, 'stabil_params')
    %--------------------------------------------------------------------------
    %----------------------------recR2_baseline--------------------------------
    stabil_params = {};
    As = {};

    len = len_R2;

    for i = 50:(len+50) 
        [w,A,C,SBC,FPE,th]=arfit(EEG.data(ch_ind,:,i)',1,1); % transponowaæ
        [S,Serr,per,tau,exctn]=armode(A,C,th); %- ok
        stabil_params{i} = abs(eig(S));
        As{i} = A;
    end

    nameA = "A_" + channels;
    nameSP = "SP_" + channels;
    fileA = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\with_baseline\recR2_kombinacje_rozne\A\" + nameA;
    fileSP = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\with_baseline\recR2_kombinacje_rozne\SP\" + nameSP;
    save(fileA, 'As')
    save(fileSP, 'stabil_params')
    %--------------------------------------------------------------------------
    %----------------------------recR8_baseline--------------------------------
    stabil_params = {};
    As = {};

    len = len_R8;

    for i = 50:(len+50) 
        [w,A,C,SBC,FPE,th]=arfit(EEG.data(ch_ind,:,i)',1,1); % transponowaæ
        [S,Serr,per,tau,exctn]=armode(A,C,th); %- ok
        stabil_params{i} = abs(eig(S));
        As{i} = A;
    end

    nameA = "A_" + channels;
    nameSP = "SP_" + channels;
    fileA = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\with_baseline\recR8_kombinacje_rozne\A\" + nameA;
    fileSP = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\with_baseline\recR8_kombinacje_rozne\SP\" + nameSP;
    save(fileA, 'As');
    save(fileSP, 'stabil_params');
    %--------------------------------------------------------------------------
    %----------------------------recS3456_baseline--------------------------------
    stabil_params = {};
    As = {};

    len = len_S3456;

    for i = 50:(len+50) 
        [w,A,C,SBC,FPE,th]=arfit(EEG.data(ch_ind,:,i)',1,1); % transponowaæ
        [S,Serr,per,tau,exctn]=armode(A,C,th); %- ok
        stabil_params{i} = abs(eig(S));
        As{i} = A;
    end

    nameA = "A_" + channels;
    nameSP = "SP_" + channels;
    fileA = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\with_baseline\recS3456_kombinacje_rozne\A\" + nameA;
    fileSP = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\with_baseline\recS3456_kombinacje_rozne\SP\" + nameSP;
    save(fileA, 'As')
    save(fileSP, 'stabil_params')
    %--------------------------------------------------------------------------


    %------------------------------rec_no_baseline-----------------------------
    EEG = pop_rmbase(EEG, [-200 0]);
    EEG = eeg_checkset( EEG );
    stabil_params = {};
    As = {};

    for i = 1:len 
        [w,A,C,SBC,FPE,th]=arfit(EEG.data(ch_ind,:,i)',1,1); % transponowaæ
        [S,Serr,per,tau,exctn]=armode(A,C,th); %- ok
        stabil_params{i} = abs(eig(S));
        As{i} = A;
    end

    nameA = "A_" + channels;
    nameSP = "SP_" + channels;
    fileA = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\no_baseline\rec_kombinacje_rozne\A\" + nameA;
    fileSP = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\no_baseline\rec_kombinacje_rozne\SP\" + nameSP;
    save(fileA, 'As')
    save(fileSP, 'stabil_params')

    %--------------------------------------------------------------------------
    %----------------------------recR1_no_baseline--------------------------------
    stabil_params = {};
    As = {};

    size_n = size(EEG.data);
    len = size_n(3);

    for i = 1:len 
        [w,A,C,SBC,FPE,th]=arfit(EEG.data(ch_ind,:,i)',1,1); % transponowaæ
        [S,Serr,per,tau,exctn]=armode(A,C,th); %- ok
        stabil_params{i} = abs(eig(S));
        As{i} = A;
    end

    nameA = "A_" + channels;
    nameSP = "SP_" + channels;
    fileA = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\no_baseline\recR1_kombinacje_rozne\A\" + nameA;
    fileSP = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\no_baseline\recR1_kombinacje_rozne\SP\" + nameSP;
    save(fileA, 'As')
    save(fileSP, 'stabil_params')
    %--------------------------------------------------------------------------
    %----------------------------recR2_no_baseline--------------------------------
    stabil_params = {};
    As = {};

    len = len_R2;

    for i = 50:(len+50) 
        [w,A,C,SBC,FPE,th]=arfit(EEG.data(ch_ind,:,i)',1,1); % transponowaæ
        [S,Serr,per,tau,exctn]=armode(A,C,th); %- ok
        stabil_params{i} = abs(eig(S));
        As{i} = A;
    end

    nameA = "A_" + channels;
    nameSP = "SP_" + channels;
    fileA = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\no_baseline\recR2_kombinacje_rozne\A\" + nameA;
    fileSP = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\no_baseline\recR2_kombinacje_rozne\SP\" + nameSP;
    save(fileA, 'As')
    save(fileSP, 'stabil_params')
    %--------------------------------------------------------------------------
    %----------------------------recR8_no_baseline--------------------------------
    stabil_params = {};
    As = {};

    len = len_R8;

    for i = 50:(len+50) 
        [w,A,C,SBC,FPE,th]=arfit(EEG.data(ch_ind,:,i)',1,1); % transponowaæ
        [S,Serr,per,tau,exctn]=armode(A,C,th); %- ok
        stabil_params{i} = abs(eig(S));
        As{i} = A;
    end

    nameA = "A_" + channels;
    nameSP = "SP_" + channels;
    fileA = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\no_baseline\recR8_kombinacje_rozne\A\" + nameA;
    fileSP = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\no_baseline\recR8_kombinacje_rozne\SP\" + nameSP;
    save(fileA, 'As')
    save(fileSP, 'stabil_params')
    %--------------------------------------------------------------------------
    %----------------------------recS3456_no_baseline--------------------------------
    stabil_params = {};
    As = {};

    len = len_S3456;

    for i = 50:(len+50) 
        [w,A,C,SBC,FPE,th]=arfit(EEG.data(ch_ind,:,i)',1,1); % transponowaæ
        [S,Serr,per,tau,exctn]=armode(A,C,th); %- ok
        stabil_params{i} = abs(eig(S));
        As{i} = A;
    end

    nameA = "A_" + channels;
    nameSP = "SP_" + channels;
    fileA = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\no_baseline\recS3456_kombinacje_rozne\A\" + nameA;
    fileSP = "C:\Users\Asus\Desktop\licencjat\zadanie2_lic\uzyskane_dane\spotkanie220420\zmienne_flanker\no_baseline\recS3456_kombinacje_rozne\SP\" + nameSP;
    save(fileA, 'As')
    save(fileSP, 'stabil_params')
    %--------------------------------------------------------------------------

    sprintf("Policzono i zapisano dane dla kana³ów: " + channels)

end