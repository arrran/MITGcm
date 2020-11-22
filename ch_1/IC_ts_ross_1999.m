clear all

cd ('/mnt/data/malyaral/MITgcm-checkpoint66j/exp/ross_5_8/run')
%%

 dt=75;
startdate=datenum(1992,02,01);



clear time_char
for i=1:550    
   time_options(i,1)=8064*(i-1); 
   
   time_char(i,:)=[year(startdate+time_options(i,1)*dt/(3600*24)) , month(startdate+time_options(i,1)*dt/(3600*24)), day(startdate+time_options(i,1)*dt/(3600*24)) ];
end

nx=210;
ny=200;
nz=70;

%% 30 Jan 1999 = 366
NF=366;
time=time_options(NF);

t_cavity=rdmds('T',time);
salt_cavity=rdmds('S',time);





%%
 cd ('/mnt/data/malyaral/MITgcm-checkpoint66j/exp/ross_7/input');
 acc = 'real*8';
 fid=fopen('T_ross_1999.box','w','b'); fwrite(fid,t_cavity,acc);fclose(fid); 
 fid=fopen('S_ross_1999.box','w','b'); fwrite(fid,salt_cavity,acc);fclose(fid);

%%


