clear all, close all, clc

fpath = './';
files = dir('*.mat');

for fi = 1:length(files),
    
    load([ fpath '/' files(fi).name ])
    
    %% Burst
    
    if isfield(Data,'Burst_VelBeam1'),
        
        figure(1), clf
        
        for i = 1:4,
            
            ax(i) = subplot(5,1,i);
            
            pcolor(Data.Burst_MatlabTimeStamp-datenum(2014,0,0), double(Data.Burst_Range), double(eval(['Data.Burst_VelBeam' num2str(i)]))' ),
            shading flat,
            datetick
            set(gca,'YDir','reverse')
            ylabel(['Vel Beam ' num2str(i)])
            caxis([-2 2])
            colorbar,
            
            
        end
        
        if isfield(Data,'IBurst_VelBeam5'),
            
            ax(5) = subplot(5,1,5);
            pcolor(Data.IBurst_MatlabTimeStamp-datenum(2014,0,0), double(Data.IBurst_Range), double(Data.IBurst_VelBeam5)' ),
            shading flat,
            datetick
            set(gca,'YDir','reverse')
            ylabel(['Vel Beam 5'])
            caxis([-2 2])
            colorbar,
            
        else
        end
        
        linkaxes(ax,'x')
        
        print('-dpng',[files(fi).name(1:end-4) '_Burst.png'])
        
    else
    end
    
    
    %% Average
    
    
    if isfield(Data,'Average_VelX'),
        
        figure(2), clf
        
        ax(i) = subplot(3,1,1);
        
        speed = sqrt( double(Data.Average_VelX).^2 + double(Data.Average_VelY).^2 );
        pcolor(Data.Average_MatlabTimeStamp, double(Data.Average_Range), speed' ),
        shading flat,
        datetick
        set(gca,'YDir','reverse')
        ylabel(['Avg Vel Mag '])
        caxis([-2 2])
        colorbar,
        
        if isfield(Data,'Average_AltimeterIceAST'),
            hold on,
            plot( Data.Average_MatlabTimeStamp, Data.Average_AltimeterIceLE,'k.')
            plot( Data.AverageIce_MatlabTimeStamp, -Data.AverageIce_VelZ1,'r.')
            legend('Vel','IceLE','-IceVelZ1')
            
            ax(2) = subplot(3,1,2);
            plot( Data.AverageIce_MatlabTimeStamp, Data.AverageIce_VelX,'gd'),
            hold on
            plot( Data.AverageIce_MatlabTimeStamp, Data.AverageIce_VelY,'r+'),
            ylabel('Ice vel [m/s]')
            datetick
            legend('VelX','VelY','Location','NorthEastOutside')
            
            ax(3) = subplot(3,1,3);
            plot( Data.AverageIce_MatlabTimeStamp, Data.AverageIce_DistanceBeam1,'c*'),
            hold on
            plot( Data.AverageIce_MatlabTimeStamp, Data.AverageIce_DistanceBeam2,'ys'),
            plot( Data.AverageIce_MatlabTimeStamp, Data.AverageIce_DistanceBeam3,'go'),
            plot( Data.AverageIce_MatlabTimeStamp, Data.AverageIce_DistanceBeam4,'mx'),
            ylabel('Ice dist [m]')
            datetick
            legend('beam 1','beam 2','beam 3','beam 4','Location','NorthEastOutside')

        else
        end
        
        %linkaxes(ax,'x')
        
        print('-dpng',[files(fi).name(1:end-4) '_Average.png'])
        
    else
    end
    
    
    
end

