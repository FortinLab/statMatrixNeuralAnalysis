source_dir = uigetdir([]); 
d = dir([source_dir, '\*.continuous']);
n=length(d);

if n>=78 && n<=86
    
for channelNumber = 1:35
    tic
    CreateSpectrogram(channelNumber,recStart,fileDir)
    
    if channelNumber <=4
        tet = 22;
    elseif channelNumber > 4 && channelNumber <=8
        tet = 21;
    elseif channelNumber > 8 && channelNumber <=12
        tet = 17;
    elseif channelNumber > 12 && channelNumber <=16
        tet = 18;
    elseif channelNumber > 16 && channelNumber <=20
        tet = 19;
    elseif channelNumber > 20 && channelNumber <=24
        tet = 20;
    elseif channelNumber > 24 && channelNumber <=28
        tet = 24;
    elseif channelNumber > 28 && channelNumber <=32
        tet = 23;
    end
    
    folder = strcat('C:\Users\Gabo\Documents\MATLAB\GlennSpec_NCP04_2017-10-16\');
    
    pngFileName = sprintf('Tetrode%d_Channel%d.png', tet,channelNumber);
	fullFileName = fullfile(folder, pngFileName);
    
    saveas(gcf,fullFileName)
    toc
end
    
    
elseif n >= 100
    
for channelNumber = 1:64
    tic
    CreateSpectrogram(channelNumber,recStart,fileDir)
    
    if channelNumber <= 4
        tet = 14;
    elseif channelNumber > 4 && channelNumber <=8
        tet = 13;
    elseif channelNumber > 8 && channelNumber <=12
        tet = 9;
    elseif channelNumber > 12 && channelNumber <=16
        tet = 10;
    elseif channelNumber > 16 && channelNumber <=20
        tet = 11;
    elseif channelNumber > 20 && channelNumber <=24
        tet = 12;
    elseif channelNumber > 24 && channelNumber <=28
        tet = 16;
    elseif channelNumber > 28 && channelNumber <=32
        tet = 15;
    elseif channelNumber > 32 && channelNumber <=36
        tet = 22;
    elseif channelNumber > 36 && channelNumber <=40
        tet = 21;
    elseif channelNumber > 40 && channelNumber <=44
        tet = 17;
    elseif channelNumber > 44 && channelNumber <=48
        tet = 18;
    elseif channelNumber > 48 && channelNumber <=52
        tet = 19;
    elseif channelNumber > 52 && channelNumber <=56
        tet = 20;
    elseif channelNumber > 56 && channelNumber <=60
        tet = 24;
    elseif channelNumber > 60 && channelNumber <=64
        tet = 23;
    end
    
    folder = strcat('C:\Users\Gabo\Documents\MATLAB\GlennSpec_NCP04_2017-10-16\');
    
    pngFileName = sprintf('Tetrode%d_Channel%d.png', tet,channelNumber);
	fullFileName = fullfile(folder, pngFileName);
    
    saveas(gcf,fullFileName)
    toc
end
    
end
    