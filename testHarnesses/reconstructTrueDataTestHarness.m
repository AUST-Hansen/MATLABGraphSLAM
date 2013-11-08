 
% test reconstructTrueCloud

N = 1500;

PC = reconstructTrueCloud(x_obs_t_true(:,1:N),euler_obs_t(:,1:N),x_berg_cm_t(:,1:N),euler_berg_t(:,1:N),imagenex,imagenexData(:,1:N),'plot');

%PC = reconstructTrueCloud(x_obs_t_true(:,1:N),euler_obs_t(:,1:N),x_berg_cm_t(:,1:N),0*euler_berg_t(:,1:N),imagenex,imagenexData(:,1:N),'plot');