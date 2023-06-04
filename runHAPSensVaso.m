function runHAPSensVaso(y_hap,y_sens,y_vaso,t_hapIn,t_sensIn,t_vasoIn)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This is a companion code for the paper "Exploring chronic and   %%
%% transient tumor hypoxia for predicting the efficacy of hypoxia- %%
%% activated pro-drugs" by S. Mathur, S. Chen, and K.A. Rejniak    %%
%% This code simulates treatment schedule with up to 3 compounds   %%
%%                                                                 %%
%% The following parameters need to be specified:                  %%
%%  y_hap : equal 1 if HAP is administered, otherwise 0            %% 
%%  y_sens: equal 1 if sensitizer is administered, otherwise 0     %% 
%%  y_vaso: equal 1 if vasodilator is administered, otherwise 0    %%  
%%  t_sens: equal to time of sensitizer administration in          %%
%%        respect to HAP injection; can be positive or negative    %%
%%  t_vaso: equal to time od vasodilator administration in         %%
%%        respect to HAP injection; can be positive or negative    %%
%%                                                                 %%
%% It requires the following data in the DataIn/ directory:        %%
%%   GridPintsX.txt  -- x-coordinates for all grid points          %%
%%   GridPintsY.txt  -- y-coordinates for all grid points          %% 
%%   VeloUX.txt      -- x-coordinates for velocities on a grid     %%
%%   VeloUY.txt      -- y-coordinates for velocities on a grid     %%
%%   CellPoints.txt  -- x & y coordinates for all cell points      %%
%%   CellCenters.txt -- x & y coordinates for all cell centers     %%
%%   oxygen.txt      -- values of oxygen concentration on a grid   %%
%%                                                                 %%
%% June 10, 2023                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% switched that can be changed
toSave = 1;        % 1 to save data;    0 if no saving
toDraw = 1;        % 1 to draw results; 0 if no drawing
toSaveLast = 1;    % 1 to save data at the end of simulation

% parameters that can be changed
Nsave=10000;       % frequency of data saving (10,000=15 min) 
Ndraw=10000;       % frequency of data drawing (10,000=15 min)


%-------------------------------------------
% the rest of the code should not be changed

hap_yes =y_hap;    % indicates if HAP drug will be injected
sens_yes=y_sens;   % indicates if sensitizer will be injected
vaso_yes=y_vaso;   % indicated if vasodilator will be injected

t_hapOut  = t_hapIn +10;   % time to stop HAP  injection [min]
t_sensOut = t_sensIn+25;   % time to stop Vaso injection [min]
t_vasoOut = t_vasoIn+25;   % time to stop Sens injection [min]
t_schedule= 3*60;          % time to stop simulation [min]

 
dt=0.0015;                 % time step [min]
  NhapIn  = t_hapIn/dt;    % iteration to start HAP injection
  NhapOut = t_hapOut/dt;   % iteration to end HAP injection
  NsensIn = t_sensIn/dt;   % iteration to start Sens injection
  NsensOut= t_sensOut/dt;  % iteration to end Sens injection
  NvasoIn = t_vasoIn/dt;   % iteration to start Vaso injection
  NvasoOut= t_vasoOut/dt;  % iteration to end Vaso injection
Nschedule = t_schedule/dt; % total number of iterations


% oxygen parameters
DiffOxy=1000;         % diffusion coefficient [um^2/min]
Oxy_in=60;            % influx [mmHg]
Oxy_hypo=10;          % hypoxia level [mmHg]
Oxy_up=0.85;          % uptake per cell pixel [mmHg/min]

% inactive drug parameters
DiffiHAP=DiffOxy/50;  % diffusion coefficient [um^2/min]  
pHAP_in=50;           % influx [ag/um^3]
pHAP_act=0.9/dt;      % activation rate  

% acivated drug parameters
DiffaHAP=DiffOxy/25;  % diffusion coefficient [um^2/min] 
aHAP_in=0;            % influx [ag/um^3]   
aHAP_up=0.5;          % uptake [1/min] 
aHAP_dec=log(2)/10;   % decay [1/min]  
aHAP_kill=1;          % lethal level [ag/um^3]

% sensitizer parameters
DiffSens=DiffOxy;     % diffusion coefficient [um^2/min] 
Sens_in=101;          % influx [ag/um^3] 
Sens_act=[176,44,18]/2;  % levels for uptake increase [um^2/min]          
Sens_Oxup=[12.5/3,7.5/3,5/3,1];  % increase rates in O2 uptake [mmHg] 

% vasodiltor parameters
vasoRate=0.5;         % rate of influx decrease during dilation


% loading intial data from external files 
% grid point coordinates and velocity field
xgr=load('DataIn/GridPointsX.txt'); ygr=load('DataIn/GridPointsY.txt'); 
xx=xgr(:,1); yy=ygr(1,:);              % grid for drawing   
ux=load('DataIn/VeloUX.txt'); uy=load('DataIn/VeloUY.txt'); % fluid

% domain boundaries and grid 
xmin=min(min(xgr(:,1))); xmax=max(max(xgr(:,1)));  % domain boundries
ymin=min(min(ygr(1,:))); ymax=max(max(ygr(1,:))); 
Ngx=size(xgr,1); Ngy=size(xgr,2);      % number of grid points
hg=(xmax-xmin)/Ngx;                    % grid width 


% tumor cell points, centers, and cell segmentation by index
pts=load('DataIn/CellPoints.txt');     % tumor cell points
ptsNum=size(pts,1);                    % number of points
ptsUpt=zeros(ptsNum,3);                % uptake of [oxy,drug,sens]
center=load('DataIn/CellCenters.txt'); % cell centers
Ncells=size(center,1);                 % number of cells
cells=zeros(size(pts,1),1);            % cell segmentation   
dist_mat=pdist2(center(:,1:2),pts(:,1:2));    % to the nearby center
[dist,index]=min(dist_mat); cells(:,1)=index; % by point index


% oxygen, inactive drug, active drug, and sensitiser concentrations
oxy=load('DataIn/oxygen.txt');         % oxygen
drgI=zeros(size(xgr));                 % inactive drug
drgA=zeros(size(xgr));                 % activated drug
sens=zeros(size(xgr));                 % activated drug
  
% variables for diffusion calculations
oxyM=zeros(Ngx+2,Ngy+2);
oxyM(2:Ngx+1,2:Ngy+1)=oxy(1:Ngx,1:Ngy);
oxyM(1,:)   =Oxy_in;       oxyM(Ngx+2,:) =0;              
oxyM(  :,1)   =oxyM(:,2);  oxyM(  :,Ngy+2)=oxyM(:,Ngy+1); 
drgIM=zeros(Ngx+2,Ngy+2);
drgAM=zeros(Ngx+2,Ngy+2);
sensM=zeros(Ngx+2,Ngy+2);


% grid occupation by the cell points - for interior and exterior
dataInd=zeros(size(oxy));     % grid index of cell occupation
  for ii=1:size(pts,1)                 
    NptsX=1+floor((pts(ii,1)-xmin)/hg); 
    NptsY=1+floor((pts(ii,2)-ymin)/hg);
    dataInd(NptsX,NptsY)=1;
  end
cellInter=zeros(size(oxyM));  %grid points inside the cells
  for ii=1:ptsNum
    NptsX=2+floor((pts(ii,1)-xmin)/hg);  
    NptsY=2+floor((pts(ii,2)-ymin)/hg);
    cellInter(NptsX,NptsY)=1;
  end


% create a directory to save all data and parameters
if (toSave==1)||(toSaveLast==1)
  pathData=MakeDir(hap_yes,sens_yes,vaso_yes,t_hapIn,t_sensIn,t_vasoIn);

  parNames=DefineParameterNames; % save order of parameters
  fID=fopen([pathData,'/parameterNames.txt'],'w');
    for ii=1:length(parNames) fprintf(fID,'%s\n',parNames(ii)); end
    fclose(fID);
  % save parameters
  parameters=[hap_yes,vaso_yes,sens_yes,t_hapIn,t_hapOut,t_vasoIn,...
    t_vasoOut,t_sensIn,t_sensOut,t_schedule,NhapIn,NhapOut,NvasoIn,...
    NvasoOut,NsensIn,NsensOut,Nschedule,dt,DiffOxy,Oxy_in,Oxy_hypo,...
    Oxy_up,DiffiHAP,pHAP_in,pHAP_act,DiffaHAP,aHAP_in,aHAP_up,...
    aHAP_dec,aHAP_kill,vasoRate,DiffSens,Sens_in,Sens_act,Sens_Oxup,...
    Nsave,Ndraw]';
  fname=[pathData,'/parameters.txt']; save(fname,'parameters','-ascii')
end


% create a figure
if (toDraw==1)
  fig_hdl=figure('position',[100,150,1200,600]);
end


tic
%--------------------------------------------------
%  main program loop
%--------------------------------------------------
vx=ux; vy=uy;     % save interstitial fluid velocities for diation rate
for iter=0:Nschedule
  if (mod(iter,5000)==0)   
    disp(['iteration: ',num2str(iter),' out of ',num2str(Nschedule)])  
  end
 
  
  % oxygen influx under vasodilator when scheduled
  if (vaso_yes==1)&&(iter>=NvasoIn)&&(iter<=NvasoOut) % oxygen influx 
       oxyM(1:5,:)=vasoRate*Oxy_in; ux=vasoRate*vx; uy=vasoRate*vy; 
  else oxyM(1:5,:)=Oxy_in; ux=vx; uy=vy; end
  
  % inactive drug influx under vasodilator when scheduled
  if (hap_yes==1)&&(iter>=NhapIn)&&(iter<=NhapOut) % inactive drug influx 
    if (vaso_yes==1)&&(iter>=NvasoIn)&&(iter<=NvasoOut)   
         drgIM(1:5,:)=vasoRate*pHAP_in;   
    else drgIM(1:5,:)=pHAP_in; end  
  else   drgIM(1:5,:)=0; end
 
  % sensitizer influx under vasodilator when scheduled
  if (sens_yes==1)&&(iter>=NsensIn)&&(iter<=NsensOut) % sensitizer influx 
    if (vaso_yes==1)&&(iter>=NvasoIn)&&(iter<=NvasoOut)   
         sensM(1:5,:)=vasoRate*Sens_in;   
    else sensM(1:5,:)=Sens_in; end  
  else   sensM(1:5,:)=0; end
 

  % oxygen and drug boundary conditions 
  oxyM(end,:)=0; oxyM(:,1)=oxyM(:,2); oxyM(:,end)=oxyM(:,end-1); 
  drgIM(end,:)=drgIM(end-1,:); drgIM(:,1)=drgIM(:,2);
  drgIM(:,end)=drgIM(:,end-1);  
  drgAM(1,:)=drgAM(2,:); drgAM(end,:)=drgAM(end-1,:); 
  drgAM(:,1)=drgAM(:,2); drgAM(:,end)=drgAM(:,end-1);
  sensM(1,:)=sensM(2,:); sensM(end,:)=sensM(end-1,:); 
  sensM(:,1)=sensM(:,2); sensM(:,end)=sensM(:,end-1);
 
  % compound kinetics calculations
  for ii=2:Ngx+1
    for jj=2:Ngy+1
   
      % DIFFUSION  
      if (cellInter(ii,jj)==0)  %diffusion outside the cell bodies  
        lf=(1-cellInter(ii-1,jj));  rg=(1-cellInter(ii+1,jj));
        bt=(1-cellInter(ii,jj-1));  tp=(1-cellInter(ii,jj+1));
        il=lf+rg+tp+bt;       
        
        %%OXYGEN
        cent=oxyM(ii,jj);
        left=(1-cellInter(ii-1,jj))*oxyM(ii-1,jj);
        righ=(1-cellInter(ii+1,jj))*oxyM(ii+1,jj);
         bot=(1-cellInter(ii,jj-1))*oxyM(ii,jj-1);
         top=(1-cellInter(ii,jj+1))*oxyM(ii,jj+1);       
        oxyM(ii,jj)=cent+(DiffOxy*dt/(hg*hg))*(left+righ+top+bot-il*cent);
      
        %%INACTIVE DRUG
        cent=drgIM(ii,jj);
        left=(1-cellInter(ii-1,jj))*drgIM(ii-1,jj);
        righ=(1-cellInter(ii+1,jj))*drgIM(ii+1,jj);
         bot=(1-cellInter(ii,jj-1))*drgIM(ii,jj-1);
         top=(1-cellInter(ii,jj+1))*drgIM(ii,jj+1);       
        drgIM(ii,jj)=cent+(DiffiHAP*dt/(hg*hg))*(left+righ+top+bot-il*cent);

        %%ACTIVE DRUG
        cent=drgAM(ii,jj);
        left=(1-cellInter(ii-1,jj))*drgAM(ii-1,jj);
        righ=(1-cellInter(ii+1,jj))*drgAM(ii+1,jj);
         bot=(1-cellInter(ii,jj-1))*drgAM(ii,jj-1);
         top=(1-cellInter(ii,jj+1))*drgAM(ii,jj+1);       
        drgAM(ii,jj)=cent+(DiffaHAP*dt/(hg*hg))*(left+righ+top+bot-il*cent);
      
        %%SENSITIZER
        cent=sensM(ii,jj);
        left=(1-cellInter(ii-1,jj))*sensM(ii-1,jj);
        righ=(1-cellInter(ii+1,jj))*sensM(ii+1,jj);
         bot=(1-cellInter(ii,jj-1))*sensM(ii,jj-1);
         top=(1-cellInter(ii,jj+1))*sensM(ii,jj+1);       
        sensM(ii,jj)=cent+(DiffSens*dt/(hg*hg))*(left+righ+top+bot-il*cent);   
      end  % end diffusion
      
      
      % ADVECTION
      if ux(ii-1,jj-1)>=0
        oxyM(ii,jj) =oxyM(ii,jj)+(dt/hg)*(oxyM(ii-1,jj)-oxyM(ii,jj))*ux(ii-1,jj-1);
        drgIM(ii,jj)=drgIM(ii,jj)+(dt/hg)*(drgIM(ii-1,jj)-drgIM(ii,jj))*ux(ii-1,jj-1);
        drgAM(ii,jj)=drgAM(ii,jj)+(dt/hg)*(drgAM(ii-1,jj)-drgAM(ii,jj))*ux(ii-1,jj-1);
        sensM(ii,jj)=sensM(ii,jj)+(dt/hg)*(sensM(ii-1,jj)-sensM(ii,jj))*ux(ii-1,jj-1);
      else
        oxyM(ii,jj) =oxyM(ii,jj)+(dt/hg)*(oxyM(ii,jj)-oxyM(ii+1,jj))*ux(ii-1,jj-1);
        drgIM(ii,jj)=drgIM(ii,jj)+(dt/hg)*(drgIM(ii,jj)-drgIM(ii+1,jj))*ux(ii-1,jj-1);
        drgAM(ii,jj)=drgAM(ii,jj)+(dt/hg)*(drgAM(ii,jj)-drgAM(ii+1,jj))*ux(ii-1,jj-1);
        sensM(ii,jj)=sensM(ii,jj)+(dt/hg)*(sensM(ii,jj)-sensM(ii+1,jj))*ux(ii-1,jj-1);         
      end
      if uy(ii-1,jj-1)>=0
        oxyM(ii,jj) =oxyM(ii,jj)+(dt/hg)*(oxyM(ii,jj-1)-oxyM(ii,jj))*uy(ii-1,jj-1);
        drgIM(ii,jj)=drgIM(ii,jj)+(dt/hg)*(drgIM(ii,jj-1)-drgIM(ii,jj))*uy(ii-1,jj-1);
        drgAM(ii,jj)=drgAM(ii,jj)+(dt/hg)*(drgAM(ii,jj-1)-drgAM(ii,jj))*uy(ii-1,jj-1);
        sensM(ii,jj)=sensM(ii,jj)+(dt/hg)*(sensM(ii,jj-1)-sensM(ii,jj))*uy(ii-1,jj-1);
      else
        oxyM(ii,jj) =oxyM(ii,jj)+(dt/hg)*(oxyM(ii,jj)-oxyM(ii,jj+1))*uy(ii-1,jj-1);
        drgIM(ii,jj)=drgIM(ii,jj)+(dt/hg)*(drgIM(ii,jj)-drgIM(ii,jj+1))*uy(ii-1,jj-1);
        drgAM(ii,jj)=drgAM(ii,jj)+(dt/hg)*(drgAM(ii,jj)-drgAM(ii,jj+1))*uy(ii-1,jj-1);
        sensM(ii,jj)=sensM(ii,jj)+(dt/hg)*(sensM(ii,jj)-sensM(ii,jj+1))*uy(ii-1,jj-1);
      end  % end advection
      

      % ACTIVATION
      if (cellInter(ii,jj)==0)&&(oxyM(ii,jj)<=Oxy_hypo)
        drgAM(ii,jj)=drgAM(ii,jj)+dt*pHAP_act*drgIM(ii,jj);  
        drgIM(ii,jj)=max(0,drgIM(ii,jj)*(1-dt*pHAP_act));  
      end  % end activation

      
      % DECAY
      drgAM(ii,jj)=max(0,drgAM(ii,jj)*(1-dt*aHAP_dec));
       
    end 
  end  % end of the double loop for compound kinetics 
  
  
  % CELLULAR UPTAKE -- oxygen and active drug  ptsUpt=[oxy,drug,sens]
  for ii=1:ptsNum          % check all points that form all cells
    if (ptsUpt(ii,1)>=0)   % not dead
      ptsUpt(ii,1)=0;
      
      % uptake from 4 grid point around the cell point
      NptsX=2+floor((pts(ii,1)-xmin)/hg);  
      NptsY=2+floor((pts(ii,2)-ymin)/hg);
      if (cellInter(NptsX-1,NptsY-1)==0)  % if inside the cell
        if (sensM(NptsX-1,NptsY-1)>Sens_act(1)) Sscale=Sens_Oxup(1);
        elseif (sensM(NptsX-1,NptsY-1)>Sens_act(2)) Sscale=Sens_Oxup(2);
        elseif (sensM(NptsX-1,NptsY-1)>Sens_act(3)) Sscale=Sens_Oxup(3);
        else Sscale=Sens_Oxup(4); end     % calculate O2 uptake rate
        
        if (oxyM(NptsX-1,NptsY-1)<=5)     % calculate O2 uptake
          oxyTmp=oxyM(NptsX-1,NptsY-1)*(1-Sscale*Oxy_up*dt);
        else       
          oxyTmp=max(0,oxyM(NptsX-1,NptsY-1)-Sscale*Oxy_up*dt);
        end
        oxyM(NptsX-1,NptsY-1)=oxyTmp;     % on the domain grid 
        ptsUpt(ii,1)=ptsUpt(ii,1)+oxyTmp; % absorbed by cell     

        drgTmp=drgAM(NptsX-1,NptsY-1)*aHAP_up*dt;  % calculate drug uptake
        drgAM(NptsX-1,NptsY-1)=max(0,drgAM(NptsX-1,NptsY-1)-drgTmp); % grid
        ptsUpt(ii,2)=ptsUpt(ii,2)+drgTmp; % absorbed by cell
      end
    
      if (cellInter(NptsX-1,NptsY+1)==0)  % if inside the cell
        if (sensM(NptsX-1,NptsY+1)>Sens_act(1)) Sscale=Sens_Oxup(1);
        elseif (sensM(NptsX-1,NptsY+1)>Sens_act(2)) Sscale=Sens_Oxup(2);
        elseif (sensM(NptsX-1,NptsY+1)>Sens_act(3)) Sscale=Sens_Oxup(3);
        else Sscale=Sens_Oxup(4); end     % calculate O2 uptake rate
        
        if (oxyM(NptsX-1,NptsY+1)<=5)     % calculate O2 uptake 
          oxyTmp=oxyM(NptsX-1,NptsY+1)*(1-Sscale*Oxy_up*dt);
        else       
          oxyTmp=max(0,oxyM(NptsX-1,NptsY+1)-Sscale*Oxy_up*dt);
        end
        oxyM(NptsX-1,NptsY+1)=oxyTmp;     % on the domain grid 
        ptsUpt(ii,1)=ptsUpt(ii,1)+oxyTmp; % absorbed by cell
      
        drgTmp=drgAM(NptsX-1,NptsY+1)*aHAP_up*dt;  % calculate drug uptake
        drgAM(NptsX-1,NptsY+1)=max(0,drgAM(NptsX-1,NptsY+1)-drgTmp); % grid
        ptsUpt(ii,2)=ptsUpt(ii,2)+drgTmp; % absorbed by cell
      end
    
      if (cellInter(NptsX+1,NptsY-1)==0)  % if inside the cell
        if (sensM(NptsX+1,NptsY-1)>Sens_act(1)) Sscale=Sens_Oxup(1);
        elseif (sensM(NptsX+1,NptsY-1)>Sens_act(2)) Sscale=Sens_Oxup(2);
        elseif (sensM(NptsX+1,NptsY-1)>Sens_act(3)) Sscale=Sens_Oxup(3);
        else Sscale=Sens_Oxup(4); end     % calculate O2 uptake rate
        
        if (oxyM(NptsX+1,NptsY-1)<=5)     % calculate O2 uptake 
          oxyTmp=oxyM(NptsX+1,NptsY-1)*(1-Sscale*Oxy_up*dt);
        else         
          oxyTmp=max(0,oxyM(NptsX+1,NptsY-1)-Sscale*Oxy_up*dt);
        end
        oxyM(NptsX+1,NptsY-1)=oxyTmp;     % on the domain grid 
        ptsUpt(ii,1)=ptsUpt(ii,1)+oxyTmp; % absorbed by cell
      
        drgTmp=drgAM(NptsX+1,NptsY-1)*aHAP_up*dt; % calculate drug uptake
        drgAM(NptsX+1,NptsY-1)=max(0,drgAM(NptsX+1,NptsY-1)-drgTmp); % grid
        ptsUpt(ii,2)=ptsUpt(ii,2)+drgTmp; % absorbed by cell
      end
    
      if (cellInter(NptsX+1,NptsY+1)==0)  % if inside the cell
        if (sensM(NptsX+1,NptsY+1)>Sens_act(1)) Sscale=Sens_Oxup(1);
        elseif (sensM(NptsX+1,NptsY+1)>Sens_act(2)) Sscale=Sens_Oxup(2);
        elseif (sensM(NptsX+1,NptsY+1)>Sens_act(3)) Sscale=Sens_Oxup(3);
        else Sscale=Sens_Oxup(4); end     % calculate O2 uptake rate
        
        if (oxyM(NptsX+1,NptsY+1)<=5)     % calculate O2 uptake 
          oxyTmp=oxyM(NptsX+1,NptsY+1)*(1-Sscale*Oxy_up*dt);
        else       
          oxyTmp=max(0,oxyM(NptsX+1,NptsY+1)-Sscale*Oxy_up*dt);
        end
        oxyM(NptsX+1,NptsY+1)=oxyTmp;     % on the domain grid  
        ptsUpt(ii,1)=ptsUpt(ii,1)+oxyTmp; % absorbed by cell
       
        drgTmp=drgAM(NptsX+1,NptsY+1)*aHAP_up*dt; % calculate drug uptake
        drgAM(NptsX+1,NptsY+1)=max(0,drgAM(NptsX+1,NptsY+1)-drgTmp); % grid
        ptsUpt(ii,2)=ptsUpt(ii,2)+drgTmp; % absorbed by cell
      end 
    end 
  end  % end of the double loop for compound uptake 

  
  % CELL KILL
  mat_cell=zeros(Ncells,2);  % drug content per cell
  for kk=1:Ncells
    indCell=find(cells(:,1)==kk); % all points for a given cell
    if (length(indCell>0))
      drgUpt=sum(ptsUpt(indCell,2));   % all absorbed drug
      drgUptAv=drgUpt/length(indCell); % average drug
      mat_cell(kk,1:2)=[drgUpt,drgUptAv];
      if (drgUptAv>aHAP_kill)     % if drug above the threshold
        ptsUpt(indCell,1)=-100;   % the cell marked dead
      end
    end
  end
   
  
  % --------------------------
  % save the data periodically
  if (toSave==1)&&(mod(iter,Nsave)==0)
    SaveData(oxyM,drgIM,drgAM,sensM,ptsUpt,mat_cell,aHAP_kill,Ngx,Ngy,...
    iter,pathData);
  end
  % save the last data 
  if (toSaveLast==1)&&(iter==Nschedule)
    SaveData(oxyM,drgIM,drgAM,sensM,ptsUpt,mat_cell,aHAP_kill,Ngx,Ngy,...
    iter,pathData);
  end

  
  % ---------------------------------------------
  % draw the data periodically and last iteration
  if (toDraw==1)&&((mod(iter,Ndraw)==0)||(iter==Nschedule))              
    clf
    set(fig_hdl,'PaperPositionMode','auto')  
    set(fig_hdl,'InvertHardcopy','off')
    DrawFigurePanel(oxyM,drgIM,drgAM,sensM,dataInd,pts,mat_cell,cells,...
        Ncells,iter,Nschedule,xx,yy,Ngx,Ngy,dt,Oxy_hypo,aHAP_kill,hg,...
        xmin,xmax,ymin,ymax)        
    if (toSave==1)&&(mod(iter,Nsave)==0)
      fname=[pathData,'/Figs/Fig_',num2str(iter)];
      print('-djpeg',fname)
    end
    pause(0.01)        
  end
  
end  % end of the main loop
toc

clear all
end  % end of the whole program


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function creates a directory in which all data will be saved        % 
% directory name structure: HAPSensVaso_HtimeStimeVtime               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pathOut=MakeDir(hap_y,sens_y,vaso_y,t_hapI,t_sensI,t_vasoI)  
  pathOut=['HAPSensVaso_',num2str(hap_y),num2str(sens_y),num2str(vaso_y),...
           '_H',num2str(t_hapI),'S',num2str(t_sensI),'V',num2str(t_vasoI)];
  mkdir(pathOut)
  mkdir([pathOut,'/Figs'])
  mkdir([pathOut,'/Data'])
  disp(['your data will be saved in a directory: ',pathOut])  
end % end function MakeDir

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function draws four panels of the figure with (1) oxygen,           %
% (2) inactive drug, (3) acivated drug, and (4) sensitizer            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DrawFigurePanel(oxyM,drgIM,drgAM,sensM,dataInd,pts,mat_cell,...
cells,Ncells,iter,Nschedule,xx,yy,Ngx,Ngy,dt,Oxy_hypo,aHAP_kill,hg,...
xmin,xmax,ymin,ymax)
    
  oxymap=[ 0,255,255;  0,255,255;  0,215,255;  0,174,255;  0,134,255;...
           0, 94,255;  0, 40,255;  0, 13,255;  5, 11,226; 12, 27,183;...
          47, 29, 93;102,  0,  0;168,  0,  0;233,  0,  0;255, 26,  0;...
         255, 64,  0;255,102,  0;255,153,  0;255,178,  0;255,217,  0;...
         255,255,  0;255,255,225]/255;
  cellcol=[139,69,19]/255;  
  
  
  %% oxygen
  % draw oxy average profiles
  ax11=subplot('position',[0.08,0.9,0.3,0.1]); 
    col=colormap(ax11,oxymap);   colSz=size(col,1);
    datamax=max(max(oxyM)); rate=datamax/colSz;
    dataProfile=zeros(Ngx,2);
    for ii=1:Ngx
      for jj=1:Ngy  
        if (dataInd(ii,jj)==0)     
          dataProfile(ii,1)=dataProfile(ii,1)+oxyM(ii+1,jj+1);
          dataProfile(ii,2)=dataProfile(ii,2)+1;
        end
      end
    end
    for ii=1:Ngx
      dataProfile(ii,1)=dataProfile(ii,1)/dataProfile(ii,2);  
    end
    hold on
     
    hypI=find(dataProfile(:,1)>=Oxy_hypo);
    xpl=hg*hypI(end); ypl=dataProfile(hypI(end),1);
    ratpl=min(colSz,1+floor(dataProfile(hypI(end),1)/rate));
    line([xpl,xpl],[0,60],'Color','k','LineWidth',2)
    for ii=1:Ngx
      xpl=hg*ii; ypl=dataProfile(ii,1);
      ratpl=min(colSz,1+floor(dataProfile(ii,1)/rate));
      colpl=[col(ratpl,1),col(ratpl,2),col(ratpl,3)]; 
      line([xpl,xpl],[0,ypl],'Color',colpl,'LineWidth',7)
    end  
    axis([xmin-5,xmax-7,0,max(2*datamax,0.1)])
    caxis([0,60])
    axis off
    
  % draw oxy tissue distribution 
  ax1=subplot('position',[0.05,0.55,0.4,0.3]);
    axis([xmin,xmax,ymin,ymax])
    axis equal
    axis([xmin,xmax,ymin,ymax])
    ylabel('blood vessel','Color','r','FontSize',15)
    hold on
    oxy(1:Ngx,1:Ngy)=oxyM(2:Ngx+1,2:Ngy+1);   
    imagesc(xx,yy,oxy')       
    hold on 
    colormap(ax1,oxymap)
    colorbar
    caxis([0,60])
    shading interp
    plot(pts(:,1),pts(:,2),'o','MarkerSize',3,'markerfaceColor',...
        [0.5,0.5,0.5],'markeredgeColor',[0.5,0.5,0.5]) 

    title(['OXYGEN;    iter=',num2str(iter),' of ',num2str(Nschedule),...
           ';       time=',num2str(iter*dt),'  min of ',...
           num2str(Nschedule*dt), ' min']) 
 
     
      
  %% inactive drug  
  % draw inactive drug average profiles
  ax21=subplot('position',[0.08,0.4,0.3,0.1]);
    col=colormap(ax21,autumn);   colSz=size(col,1);
    datamax=max(max(drgIM)); rate=datamax/colSz;
    dataProfile=zeros(Ngx,2);
    for ii=1:Ngx
      for jj=1:Ngy  
        if (dataInd(ii,jj)==0)     
          dataProfile(ii,1)=dataProfile(ii,1)+drgIM(ii+1,jj+1);
          dataProfile(ii,2)=dataProfile(ii,2)+1;
        end
      end
    end
    for ii=1:Ngx
      dataProfile(ii,1)=dataProfile(ii,1)/dataProfile(ii,2);  
    end
    hold on
      
    for ii=1:Ngx
      xpl=hg*ii; ypl=max(0.001,dataProfile(ii,1));
      ratpl=min(colSz,1+floor(dataProfile(ii,1)/rate));
      colpl=[col(ratpl,1),col(ratpl,2),col(ratpl,3)]; 
      line([xpl,xpl],[0,ypl],'Color',colpl,'LineWidth',7)
    end
    axis([xmin-5,xmax-7,0,max(2*datamax,0.1)])  
    axis off
 
  % draw inactive drug tissue distribution 
  ax2=subplot('position',[0.05,0.05,0.4,0.3]);
    axis([xmin,xmax,ymin,ymax])
    axis equal
    axis([xmin,xmax,ymin,ymax])
    ylabel('blood vessel','Color','r','FontSize',15)
    hold on  
    drgI(1:Ngx,1:Ngy)=drgIM(2:Ngx+1,2:Ngy+1);  
    imagesc(xx,yy,drgI')
    hold on 
    colormap(ax2,autumn)
    colorbar
    caxis([0,5])
    shading interp
    plot(pts(:,1),pts(:,2),'o','MarkerSize',3,'markerfaceColor',...
        [0.5,0.5,0.5],'markeredgeColor',[0.5,0.5,0.5])

    title(['INACTIVE DRUG;   max drug=',num2str(max(max(drgI(:,:)))),...
           '  total=',num2str(sum(sum(drgI)))])    
   
     
  %% active drug  
  % draw active drug average profiles
  ax31=subplot('position',[0.53,0.4,0.3,0.1]);
    col=colormap(ax31,spring);   colSz=size(col,1);
    datamax=max(max(drgAM)); rate=datamax/colSz;
    dataProfile=zeros(Ngx,2);
    for ii=1:Ngx
      for jj=1:Ngy  
        if (dataInd(ii,jj)==0)     
          dataProfile(ii,1)=dataProfile(ii,1)+drgAM(ii+1,jj+1);
          dataProfile(ii,2)=dataProfile(ii,2)+1;
        end
      end
    end
    for ii=1:Ngx
      dataProfile(ii,1)=dataProfile(ii,1)/dataProfile(ii,2);  
    end
    hold on
    for ii=1:Ngx
      xpl=hg*ii; ypl=max(0.001,3*dataProfile(ii,1));
      ratpl=min(colSz,5+floor(dataProfile(ii,1)/rate));
      colpl=[col(ratpl,1),col(ratpl,2),col(ratpl,3)]; 
      line([xpl,xpl],[0,ypl],'Color',colpl,'LineWidth',7)
    end
    axis([xmin-5,xmax-7,0,max(2*datamax,0.1)])
    axis off

  % draw active drug tissue distribution         
  ax3=subplot('position',[0.5,0.05,0.4,0.3]);
    axis([xmin,xmax,ymin,ymax])
    axis equal
    axis([xmin,xmax,ymin,ymax])
    ylabel('blood vessel','Color','r','FontSize',15)
    hold on
    drgA(1:Ngx,1:Ngy)=drgAM(2:Ngx+1,2:Ngy+1); 
    imagesc(xx,yy,drgA')
    hold on 
    colormap(ax3,spring)
    colorbar
    caxis([0,0.5])
    shading interp
    plot(pts(:,1),pts(:,2),'o','MarkerSize',3,'markerfaceColor',...
         [0.5,0.5,0.5],'markeredgeColor',[0.5,0.5,0.5])  
           
    indMat=find(mat_cell(:,2)>=aHAP_kill);   
    CellDead=length(indMat);
    for kk=1:length(indMat)
      indCell=find(cells(:,1)==indMat(kk));
      plot(pts(indCell,1),pts(indCell,2),'.','color',[0,0,0],'MarkerSize',10) 
    end     
   
    title(['ACTIVE DRUG;    max=',num2str(max(max(drgA))),...
           '  total=',num2str(sum(sum(drgA(:,:)))),...
           '  #dead cells: ',num2str(CellDead),' of ',num2str(Ncells)])  
   
        
  %% sensitizer  
  % draw sensitizer average profiles
  ax41=subplot('position',[0.53,0.9,0.3,0.1]);
    col=colormap(ax41,summer);   colSz=size(col,1);
    datamax=max(max(sensM)); rate=datamax/colSz;
    dataProfile=zeros(Ngx,2);
    for ii=1:Ngx
      for jj=1:Ngy  
        if (dataInd(ii,jj)==0)     
          dataProfile(ii,1)=dataProfile(ii,1)+sensM(ii+1,jj+1);
          dataProfile(ii,2)=dataProfile(ii,2)+1;
        end
      end
    end
    for ii=1:Ngx
      dataProfile(ii,1)=dataProfile(ii,1)/dataProfile(ii,2);  
    end
    hold on
    for ii=1:Ngx
      xpl=hg*ii; ypl=max(0.001,dataProfile(ii,1));
      ratpl=min(colSz,5+floor(dataProfile(ii,1)/rate));
      colpl=[col(ratpl,1),col(ratpl,2),col(ratpl,3)]; 
      line([xpl,xpl],[0,ypl],'Color',colpl,'LineWidth',7)
    end
    axis([xmin-5,xmax-7,0,max(2*datamax,0.1)])
    axis off
 
  % draw sensitizer tissue distribution 
  ax4=subplot('position',[0.5,0.55,0.4,0.3]);
    axis([xmin,xmax,ymin,ymax])
    axis equal
    axis([xmin,xmax,ymin,ymax])
    ylabel('blood vessel','Color','r','FontSize',15)
    hold on
    sens(1:Ngx,1:Ngy)=sensM(2:Ngx+1,2:Ngy+1); 
    imagesc(xx,yy,sens')
    hold on 
    colormap(ax4,summer)
    colorbar
    caxis([0,25])
    shading interp
    plot(pts(:,1),pts(:,2),'o','MarkerSize',3,'markerfaceColor',...
        [0.5,0.5,0.5],'markeredgeColor',[0.5,0.5,0.5])  
           
    title(['SENSITIZER;    max=',num2str(max(max(sens))),...
           '  total=',num2str(sum(sum(sens(:,:))))])
      
end % end function DrawFigurePanel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function defines names of parameters to save                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function parNames=DefineParameterNames
  parNames=[...
    string('hap_yes'),   string('vaso_yes'), string('sens_yes'),...
    string('t_hapIn'),   string('t_hapOut'), string('t_vasoIn'),...
    string('t_vasoOut'), string('t_sensIn'), string('t_sensOut'),...
    string('t_schedule'),string('NhapIn'),   string('NhapOut'),...
    string('NvasoIn'),   string('NvasoOut'), string('NsensIn'),...
    string('NsensOut'),  string('Nschedule'),string('dt'),...
    string('DiffOxy'),   string('Oxy_in'),   string('Oxy_hypo'),...
    string('Oxy_up'),    string('DiffiHAP'), string('pHAP_in'),...
    string('pHAP_act'),  string('DiffaHAP'), string('aHAP_in'),...
    string('aHAP_up'),   string('aHAP_dec'), string('aHAP_kill'),...
    string('vasoRate'),  string('DiffSens'), string('Sens_in'),...
    string('Sens_act'),  string('Sens_Oxup'),...
    string('Nsave'),string('Ndraw')]';
end % end function DefineParameterNames

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function saves relevant data                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SaveData(oxyM,drgIM,drgAM,sensM,ptsUpt,mat_cell,aHAP_kill,...
Ngx,Ngy,iter,pathData)
    
  oxy(1:Ngx,1:Ngy)=oxyM(2:Ngx+1,2:Ngy+1); 
  fname=[pathData,'/Data/oxy_',num2str(iter),'.txt'];
  save(fname,'oxy','-ASCII')
    
  drgI(1:Ngx,1:Ngy)=drgIM(2:Ngx+1,2:Ngy+1); 
  fname=[pathData,'/Data/drgI_',num2str(iter),'.txt'];
  save(fname,'drgI','-ASCII')
    
  drgA(1:Ngx,1:Ngy)=drgAM(2:Ngx+1,2:Ngy+1); 
  fname=[pathData,'/Data/drgA_',num2str(iter),'.txt'];
  save(fname,'drgA','-ASCII')

  sens(1:Ngx,1:Ngy)=sensM(2:Ngx+1,2:Ngy+1); 
  fname=[pathData,'/Data/sens_',num2str(iter),'.txt'];
  save(fname,'sens','-ASCII')

  fname=[pathData,'/Data/ptsUpt_',num2str(iter),'.txt'];
  save(fname,'ptsUpt','-ASCII')

  fname=[pathData,'/Data/matCell_',num2str(iter),'.txt'];
  save(fname,'mat_cell','-ASCII')
  
  indMat=find(mat_cell(:,2)>=aHAP_kill); 
  dead_cell=[indMat,mat_cell(indMat,:)];
  fname=[pathData,'/Data/deadCell_',num2str(iter),'.txt'];
  save(fname,'dead_cell','-ASCII')
end  % function SaveData

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%