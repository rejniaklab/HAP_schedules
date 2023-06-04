%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This is a companion code for the paper "Exploring chronic    %%
%% and transient tumor hypoxia for predicting the efficacy      %%
%% of hypoxia-activated pro-drugs" by S. Mathur, S. Chen, and   %%
%% K.A. Rejniak                                                 %%
%% This script code prepares parameters for the main function:  %%
%% runHAPSensVaso.m                                             %%
%%                                                              %%
%% The following parameters need to be specified:               %%
%%  y_hap : equal 1 if HAP is administered, otherwise 0         %% 
%%  y_sens: equal 1 if sensitizer is administered, otherwise 0  %% 
%%  y_vaso: equal 1 if vasodilator is administered, otherwise 0 %%  
%%  t_sens: equal to time od sensitizer administration in       %%
%%        respect to HAP injection; can be positive or negative %%
%%  t_vaso: equal to time of vasodilator administration in      %%
%%        respect to HAP injection; can be positive or negative %%
%%                                                              %%
%% for the example of 3-drug combination with a sensitizer 15   %%
%% minutes before HAP and a vasodilator 5 minutes after HAP:    %%
%%   y_hap=1;  y_sens=1;  y_vaso=1;  t_sens=-15;  t_vaso=5;     %%
%%                                                              %%
%% June 10, 2023                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 % five parameters to change:
 y_hap  =  1;     % use 1, if HAP is injected, otherwise 0
 y_sens =  1;     % use 1, if sensitizer is injected, otherwise 0
 y_vaso =  1;     % use 1, if vasodilator is injected, otherwise 0
 t_sens =-15;     % time of Sens injection, negative before HAP
 t_vaso =  5;     % time of Vaso injection, positive after HAP
 
 
 % automatic assessment of the time of HAP injection, so that the first 
 % injection of HAP, Sens, or Vaso is at time 0  
 t_hap=min(min(t_vaso,t_sens),0); 
 times_HSV=[y_hap,y_sens,y_vaso,0-t_hap,t_sens-t_hap,t_vaso-t_hap]
 
 
 % run the code for treatment administration
 runHAPSensVaso(y_hap,y_sens,y_vaso,0-t_hap,t_sens-t_hap,t_vaso-t_hap);
 
 %%%%%%%%%%%%%%%%
 
 
 