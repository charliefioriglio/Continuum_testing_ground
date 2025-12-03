function figplot(X,Y,Z,F, axval,iso,no,LR,type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%      Routine to Plot Dyson Orbitals at specified iso value      %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%              Last Modified:    Jan 29 2022                      %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if type==1
    hold on
subplot(1,no,LR);
isoval=iso;
 axis([-axval axval -axval axval -axval axval], 'off');
 h=gca;
 h.XAxis.TickLength=[0 0];
 h.YAxis.TickLength=[0 0];
 h.ZAxis.TickLength=[0 0];
 h.XAxis.TickLabel=[];
 h.YAxis.TickLabel=[];
 h.ZAxis.TickLabel=[];
 p=patch(isosurface(X,Y,Z,F,isoval));
 p.FaceColor = 'none'; 
 p.EdgeColor = 'red';
 p.FaceAlpha = 0.3;
 hold on
 q=patch(isosurface(X,Y,Z,F,-isoval));
 q.FaceColor = 'none';
 q.EdgeColor = 'blue';
 q.FaceAlpha = 0.3;
 grid on
 daspect([1 1 1])
 view([0.1,0.1,0.05])
 x1=[0,0];
 y1=[0,0];
 z1=[-8,8];
 plot3(x1,y1,z1,'.-black');
 x1=[-4,4];
 y1=[0,0];
 z1=[0,0];
 plot3(x1,y1,z1,'.-black');
 x1=[0,0];
 y1=[-4,4];
 z1=[0,0];
 plot3(x1,y1,z1,'.-black');
 
 end
end