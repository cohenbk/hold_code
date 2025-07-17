%Captures movie frame
eval('cd figs/')
figure(1)
currFrame=getframe(1);
writeVideo(vidObj_Z,currFrame);
figure(2)
currFrame=getframe(2);
writeVideo(vidObj_vr,currFrame);
figure(3)
currFrame=getframe(3);
writeVideo(vidObj_rhohv,currFrame);
figure(4)
currFrame=getframe(4);
writeVideo(vidObj_ZDR,currFrame);
% figure(4)
% currFrame=getframe(5);
% writeVideo(vidObj_phidp,currFrame);

eval('cd ..')
