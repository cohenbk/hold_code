linewidth=3;
if(strcmp(date_save,'20110524'))
    if(strcmp(tornado,'Goldsby'))
        if(radar_choice(2)==1)
            for ldx=1:max(size(xf0))
                xf0{ldx}=xf0{ldx}+KOUN_adjust(1); yf0{ldx}=yf0{ldx}+KOUN_adjust(2);
            end
            if(plot_all_flag)
                for ldx=1:max(size(xf1))
                    xf1{ldx}=xf1{ldx}+KOUN_adjust(1); yf1{ldx}=yf1{ldx}+KOUN_adjust(2);
                end
                for ldx=1:max(size(xf2))
                    xf2{ldx}=xf2{ldx}+KOUN_adjust(1); yf2{ldx}=yf2{ldx}+KOUN_adjust(2);
                end
                for ldx=1:max(size(xf3))
                    xf3{ldx}=xf3{ldx}+KOUN_adjust(1); yf3{ldx}=yf3{ldx}+KOUN_adjust(2);
                end
                for ldx=1:max(size(xf4))
                    xf4{ldx}=xf4{ldx}+KOUN_adjust(1); yf4{ldx}=yf4{ldx}+KOUN_adjust(2);
                end
            end
        end
        hold on
        for idx=1:5
            if(idx==1)
                for jdx=1:max(size(xf0))
                    if(white_flag)
                        plot(xf0{jdx},yf0{jdx},'-w','LineWidth',linewidth);
                    else
                        plot(xf0{jdx},yf0{jdx},'-k','LineWidth',linewidth);
                    end
                end
            elseif(idx==2 && plot_all_flag)
                for jdx=1:max(size(xf1))
                    plot(xf1{jdx},yf1{jdx},'-g','LineWidth',linewidth);
                end
            elseif(idx==3 && plot_all_flag)
                for jdx=1:max(size(xf2))
                    plot(xf2{jdx},yf2{jdx},'-y','LineWidth',linewidth);
                end
            elseif(idx==4 && plot_all_flag)
                for jdx=1:max(size(xf3))
                    plot(xf3{jdx},yf3{jdx},'-r','LineWidth',linewidth);
                end
            elseif(idx==5 && plot_all_flag)
                for jdx=1:max(size(xf4))
                    plot(xf4{jdx},yf4{jdx},'-m','LineWidth',linewidth+1);
                end
            end
        end
        hold off
    elseif(strcmp(tornado,'Chickasha'))
        if(radar_choice(2)==1)
            for ldx=1:max(size(xf0))
                xf0{ldx}=xf0{ldx}+KOUN_adjust(1); yf0{ldx}=yf0{ldx}+KOUN_adjust(2);
            end
            for ldx=1:max(size(xf1))
                xf1{ldx}=xf1{ldx}+KOUN_adjust(1); yf1{ldx}=yf1{ldx}+KOUN_adjust(2);
            end
            for ldx=1:max(size(xf2))
                xf2{ldx}=xf2{ldx}+KOUN_adjust(1); yf2{ldx}=yf2{ldx}+KOUN_adjust(2);
            end
            for ldx=1:max(size(xf3))
                xf3{ldx}=xf3{ldx}+KOUN_adjust(1); yf3{ldx}=yf3{ldx}+KOUN_adjust(2);
            end
            for ldx=1:max(size(xf4))
                xf4{ldx}=xf4{ldx}+KOUN_adjust(1); yf4{ldx}=yf4{ldx}+KOUN_adjust(2);
            end
            for ldx=1:max(size(xf4conf))
                xf4conf{ldx}=xf4conf{ldx}+KOUN_adjust(1); yf4conf{ldx}=yf4conf{ldx}+KOUN_adjust(2);
            end
        end
        hold on
        for idx=1:6
            if(idx==1)
                for jdx=1:max(size(xf0))
                    if(white_flag)
                        plot(xf0{jdx},yf0{jdx},'-w','LineWidth',linewidth);
                    else
                        plot(xf0{jdx},yf0{jdx},'-k','LineWidth',linewidth);
                    end
                end
            elseif(idx==2 && plot_all_flag)
                for jdx=1:max(size(xf1))
                    plot(xf1{jdx},yf1{jdx},'-g','LineWidth',linewidth);
                end
            elseif(idx==3 && plot_all_flag)
                for jdx=1:max(size(xf2))
                    plot(xf2{jdx},yf2{jdx},'-y','LineWidth',linewidth);
                end
            elseif(idx==4 && plot_all_flag)
                for jdx=1:max(size(xf3))
                    plot(xf3{jdx},yf3{jdx},'-r','LineWidth',linewidth);
                end
            elseif(idx==5 && plot_all_flag)
                for jdx=1:max(size(xf4))
                    plot(xf4{jdx},yf4{jdx},'-m','LineWidth',linewidth);
                end
            elseif(idx==6 && plot_all_flag)
                for jdx=1:max(size(xf4conf))
                    plot(xf4conf{jdx},yf4conf{jdx},'-b','LineWidth',linewidth);
                end
            end
        end
        hold off
    else
%         if(~exist('storm_num','var'))
%             storm_num=input(['Enter tornado number: ']);
%         end
%         eval('cd tornado_path')
%             if(storm_num==1)
%                 load(['Elreno.mat']);
%             elseif(storm_num==7)
%                 load('Goldsby_sat.mat');
%             elseif(storm_num==4)
%                 load('Mcloud.mat');
%             end
%         eval('cd ..')
%         if(radar_choice(2)==1)
%             damage_x=damage_x+KOUN_adjust(1); damage_y=damage_y+KOUN_adjust(2);
%         end
        [s,a]=shaperead('110524_tornadoes.shp');
        for idx=1:max(size(s))
            [x,y]=latlon2cart(s(idx).Y,s(idx).X,lat,lon);
            damage_x{idx}=x; damage_y{idx}=y;
            hold on
        %     subplot(3,2,[1 2 3 4])
            plot(damage_x{idx},damage_y{idx},'-k','LineWidth',linewidth)
            hold off
        end
%         hold on
%             plot(damage_x,damage_y,'-k','LineWidth',linewidth_damage)
%         hold off
    end
elseif(strcmp(date_save,'20100510'))
    if(~exist('storm_num','var'))
        storm_num=input(['Enter tornado number: ']);
    end
%     eval('cd tornado_path')
%         if(storm_num<7 && storm_num~=5)
%             load(['tornado' num2str(storm_num) '.mat']);
%             damage_x=x; damage_y=y; clear x y;
%         else
%             damage_x=NaN; damage_y=NaN;
%             eval('cd ..')
%     if(~ouprime_flag)
%             load KOUN_tornado.mat;
%             hold on
%             for qdx=1:max(size(x))
%                 plot(x{qdx},y{qdx},'-k','LineWidth',linewidth_damage)
% %                 text(x{qdx}(1),y{qdx}(1),num2str(qdx))
%             end
%             hold off
%     else
        eval('cd tornado_path')
            if(ouprime_flag)
                load(['OUP_tornado' num2str(storm_num) '.mat']);
            else
                load(['KOUN_tornado' num2str(storm_num) '.mat']);
            end
            damage_x=x; damage_y=y; clear x y;
            if(KOUN_adjust_flag && radar_choice==1)
               damage_x=damage_x+2.6; damage_y=damage_y-6.2; 
            end
            hold on
            plot(damage_x,damage_y,'-k','LineWidth',linewidth_damage)
            hold off
        eval('cd ..')
%     end
%             eval('cd tornado_path')
%         end
%     eval('cd ..')
%     hold on
%         if(radar_choice(2)==1)
%             damage_x=damage_x+KOUN_adjust(1); damage_y=damage_y+KOUN_adjust(2);
%         end
%         plot(damage_x,damage_y,'-k','LineWidth',linewidth_damage)
%     hold off
end
