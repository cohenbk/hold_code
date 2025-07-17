%options==13
if(~exist('input_fl','var'))
    input_fl=input('High Figure quality (1: Yes, 0: No)');
end
%Saves figures
title_str={'Z','vr','rhohv','ZDR','thres','Phidp','','','','Phidp_var','Zdr_var'};
if(~exist('date_save','var'))
    date_save=file_tmp(1:8);
end
figure(1)
picname=[num2str(date_save) '_' time_choice '_' num2str(ele_save(scan_num)) '_' title_str{plot_choice(1)} ];
eval('cd figs/')
if(input_fl)
%                 set(gca,'TitleFontSizeMultiplier',1.5)
    xl=get(gca,'XLabel');
    yl=get(gca,'YLabel');
    tl=get(gca,'Title');
    set(gca,'FontSize',14)
    set(xl,'FontSize',16);
    set(yl,'FontSize',16);
    set(tl,'FontSize',18);
    set(txt1(1),'FontSize',16)
    set(mark1(1),'MarkerSize',marker_size)
    set(gcf,'PaperPositionMode','auto')
    print(picname,'-dpng','-r300');
else
    print(picname,'-dpng');
end
eval('cd ..')
figure(2)
picname=[num2str(date_save) '_' time_choice '_' num2str(ele_save(scan_num)) '_' title_str{plot_choice(2)} ];
eval('cd figs')
if(input_fl)
    xl=get(gca,'XLabel');
    yl=get(gca,'YLabel');
    tl=get(gca,'Title');
    set(gca,'FontSize',14)
    set(xl,'FontSize',16);
    set(yl,'FontSize',16);
    set(tl,'FontSize',18);
    set(txt1(2),'FontSize',16)
    set(mark1(2),'MarkerSize',marker_size)
    set(gcf,'PaperPositionMode','auto')
    print(picname,'-dpng','-r300');
else
    print(picname,'-dpng');
end
eval('cd ..')
if(big_monitor)
    figure(3)
    picname=[num2str(date_save) '_' time_choice '_' num2str(ele_save(scan_num)) '_' title_str{plot_choice(3)} ];
    eval('cd figs')
    if(input_fl)
        xl=get(gca,'XLabel');
        yl=get(gca,'YLabel');
        tl=get(gca,'Title');
        set(gca,'FontSize',14)
        set(xl,'FontSize',16);
        set(yl,'FontSize',16);
        set(tl,'FontSize',18);
        set(txt1(3),'FontSize',16)
        set(mark1(3),'MarkerSize',marker_size)
        set(gcf,'PaperPositionMode','auto')
        print(picname,'-dpng','-r300');
    else
        print(picname,'-dpng');
    end
    eval('cd ..')
    figure(4)
    picname=[num2str(date_save) '_' time_choice '_' num2str(ele_save(scan_num)) '_' title_str{plot_choice(4)} ];
    eval('cd figs')
    if(input_fl)
         xl=get(gca,'XLabel');
        yl=get(gca,'YLabel');
        tl=get(gca,'Title');
        set(gca,'FontSize',14)
        set(xl,'FontSize',16);
        set(yl,'FontSize',16);
        set(tl,'FontSize',18);
        set(txt1(4),'FontSize',16)
        set(mark1(4),'MarkerSize',marker_size)
        set(gcf,'PaperPositionMode','auto')
        print(picname,'-dpng','-r300');
    else
        print(picname,'-dpng');
    end
    eval('cd ..')
    figure(5)
    picname=[num2str(date_save) '_' time_choice '_' num2str(ele_save(scan_num)) '_' title_str{plot_choice(5)} ];
    eval('cd figs')
    if(input_fl)
        xl=get(gca,'XLabel');
        yl=get(gca,'YLabel');
        tl=get(gca,'Title');
        set(gca,'FontSize',14)
        set(xl,'FontSize',16);
        set(yl,'FontSize',16);
        set(tl,'FontSize',18);
        set(txt1(5),'FontSize',16)
        set(mark1(5),'MarkerSize',marker_size)
        set(gcf,'PaperPositionMode','auto')
        print(picname,'-dpng','-r300');
    else
        print(picname,'-dpng');
    end
    eval('cd ..')
end
