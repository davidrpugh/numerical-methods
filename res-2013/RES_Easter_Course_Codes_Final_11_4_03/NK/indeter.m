clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%MODEL SPECIFIC BLOCK%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%CHOOSE THE MOD FILE NAME%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model_string='NK_Rotemberg_free_parameters_RES_Course';%%%%%%%%%%%%%%%%%%
%model_string='NK_free_parameters_RES_Course';%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%CHOOSE THE NAME OF THE TWO PARAMETERS FOR THE LOOP%%%%%%%%%%%%
parameter1_string=['\theta_{\pi}'];%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameter2_string=['\Pi_{ss}'];%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%CHOOSE THE GRID FOR EACH PARAMETER IN THE LOOP%%%%%%%%%%%%%%%%
parameter1_min=0;%starting value%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameter1_step=0.5;%step for the grid%%%%%%%%%%%%%%%%%%%%%%%%%%
parameter1_max=2.5;%maximum value%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%CHOOSE THE GRID FOR EACH PARAMETER IN THE LOOP%%%%%%%%%%%%%%%%
%parameter1_min=0.9;%starting value%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameter1_step=0.1;%step for the grid%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameter1_max=1.5;%maximum value%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameter2_min=1.0;%starting value%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameter2_step=0.01;%step for the grid%%%%%%%%%%%%%%%%%%%%%%%%%%
parameter2_max=1.03;%maximum value%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%END OF THE MODEL SPECIFIC BLOCK%%%%%%%%%%%%%%%%%%%%%%%%%%


%% DO NOT CHANGE THIS PART.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%DEFINE PARAMETERS GRIDS
parameter1=[parameter1_min:parameter1_step:parameter1_max];
parameter2=[parameter2_min:parameter2_step:parameter2_max];
for j=1:length(parameter1);
      for i=1:length(parameter2);
        %%%CALL DYNARE FOR EACH PARAMETER VALUE
          model=['dynare ' model_string ' noclearall'];

             
    eval(model)

    
    if oo_.dr.nsfwrd==oo_.dr.rank            %%IF # OF NON-PREDETERMINED VARIABLES == RANK -->UNIQUE EQUILIBRIUM, SADDLE PATH STABLE
        ind(j,i)=0;
        
    elseif oo_.dr.nsfwrd<oo_.dr.rank         %%IF # OF NON-PREDETERMINED VARIABLES  < RANK -->INDETERMINACY
    
        ind(j,i)=3;
        
    elseif oo_.dr.nsfwrd>oo_.dr.rank         %%IF # OF NON-PREDETERMINED VARIABLES  > RANK -->INSTABILITY
        ind(j,i)=4;
        
        pause(2)
    end
    
      end
end
   
save indet.mat  %SAVES THE RESULTS SO THERE IS NO NEED TO RUN THE LOOP AGAIN

%%PLOT THE REGIONS OF SADDLE PATH STABILITY, INDETERMINACY AND INSTABILITY
figure(1)
markersize=1;
spy(ind(:,:)==0,'r.');
hold on 
spy(ind(:,:)==4,'k.');
hold on 
spy(ind(:,:)==3,'g.');
axis xy;
set(gca,'XTick',1:1:length(parameter2));
set(gca,'XTicklabel',{parameter2_min:parameter2_step:parameter2_max});
set(gca,'YTick',1:1:length(parameter1));
set(gca,'YTicklabel',{parameter1_min:parameter1_step:parameter1_max});
title('Indeterminacy NK','fontsize',9); xlabel(parameter2_string,'fontsize',9);
ylabel(parameter1_string,'fontsize',9);    
legend('Determinacy','Indeterminacy','Instability',-1);