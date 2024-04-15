clear all
clc
close all
start
filename = 'mag_10t.nc';
%filename = 'coq_2hrs.nc';
%filename = 'croco_Coq_his_sml.nc';
nc = netcdf(filename);
w = nc{'w'}(:,:,:,:).*nc{'w'}.scale_factor(:) + nc{'w'}.add_offset(:);
rho = nc{'rho'}(:,:,:,:).*nc{'rho'}.scale_factor(:) + nc{'rho'}.add_offset(:);
time = nc{'time'}(:);
zeta = nc{'zeta'}(:,:,:);
theta_s = nc.theta_s(:);
theta_b = nc.theta_b(:);
hc = nc.hc(:);
N = length(nc('s_rho'));
type = 'r';
vtransform = 2;
nc2 = netcdf('mag_grd.nc');
h = nc2{'h'}(:,:);
gridz = zlevs(h,squeeze(zeta(1,:,:)).*nc{'zeta'}.scale_factor(:) + nc{'zeta'}.add_offset(:),theta_s,theta_b,hc,N,type,vtransform);

% u y v a rho-points
u = nc{'u'}(:,:,:,:).*nc{'u'}.scale_factor(:) + nc{'u'}.add_offset(:);
v = nc{'v'}(:,:,:,:).*nc{'v'}.scale_factor(:) + nc{'v'}.add_offset(:);

% Calcular u_rho y v_rho
u_rho = zeros(size(rho));
v_rho = zeros(size(rho));

for i=1:length(time) % Numero de pasos de tiempo
    for j=1:N % Número de niveles verticales
        u_aux=squeeze(u(i,j,:,:));
        v_aux=squeeze(v(i,j,:,:));
        u_rho(i,j,:,:)=u2rho_2d(u_aux);
        v_rho(i,j,:,:)=v2rho_2d(v_aux);
    end
end

%%
[time_length,s_rho,y,x]=size(u_rho);
%ncid = netcdf.create('coq_paraview_2t.nc', 'NETCDF4');
ncid = netcdf.create('mag_paraview_10t.nc', 'NETCDF4');

% Definir las dimensiones en el archivo NetCDF
time_dimid = netcdf.defDim(ncid, 'time', time_length);
lon_dimid = netcdf.defDim(ncid, 'longitude', x);
lat_dimid = netcdf.defDim(ncid, 'latitude', y);
s_rho_dimid = netcdf.defDim(ncid, 's_rho', s_rho);

% Definir las variables en el archivo NetCDF
time_varid = netcdf.defVar(ncid, 'time', 'NC_DOUBLE', [time_dimid]);
u_varid = netcdf.defVar(ncid, 'u', 'NC_FLOAT', [lon_dimid, lat_dimid, s_rho_dimid, time_dimid]);
v_varid = netcdf.defVar(ncid, 'v', 'NC_FLOAT', [lon_dimid, lat_dimid, s_rho_dimid, time_dimid]);
w_varid = netcdf.defVar(ncid, 'w', 'NC_FLOAT', [lon_dimid, lat_dimid, s_rho_dimid, time_dimid]);
rho_varid = netcdf.defVar(ncid, 'rho', 'NC_FLOAT', [lon_dimid, lat_dimid, s_rho_dimid, time_dimid]);
gridz_varid = netcdf.defVar(ncid, 'gridz', 'NC_FLOAT', [lon_dimid, lat_dimid, s_rho_dimid]);

% Finalizar la definición de variables
netcdf.endDef(ncid);

% Escribir datos en las variables
netcdf.putVar(ncid, time_varid, time);
netcdf.putVar(ncid, u_varid, permute(u_rho, [4,3,2,1]));
netcdf.putVar(ncid, v_varid, permute(v_rho, [4,3,2,1]));
netcdf.putVar(ncid, w_varid, permute(w, [4,3,2,1]));
netcdf.putVar(ncid, rho_varid, permute(rho, [4,3,2,1]));
netcdf.putVar(ncid, gridz_varid, permute(gridz, [3,2,1]));

% Re-enter define mode.
netcdf.reDef(ncid);

% Atributos para la variable 'time'
time_varid = netcdf.inqVarID(ncid, 'time');
netcdf.putAtt(ncid, time_varid, 'long_name', 'time since initialization');
netcdf.putAtt(ncid, time_varid, 'units', 'second');
netcdf.putAtt(ncid, time_varid, 'field', 'time, scalar, series');
netcdf.putAtt(ncid, time_varid, 'standard_name', 'time');
netcdf.putAtt(ncid, time_varid, 'axis', 'T');

% Atributos para las vairables
u_varid = netcdf.inqVarID(ncid, 'u');
netcdf.putAtt(ncid, u_varid, 'long_name', 'u-momentum component');
netcdf.putAtt(ncid, u_varid, 'units', 'meter second-1');
netcdf.putAtt(ncid, u_varid, 'field', 'u-velocity, scalar, series');
netcdf.putAtt(ncid, u_varid, 'standard_name', 'sea_water_x_velocity_at_u_location');
netcdf.putAtt(ncid, u_varid, 'coordinates', 'latitude longitude');

v_varid = netcdf.inqVarID(ncid, 'v');
netcdf.putAtt(ncid, v_varid, 'long_name', 'v-momentum component');
netcdf.putAtt(ncid, v_varid, 'units', 'meter second-1');
netcdf.putAtt(ncid, v_varid, 'field', 'v-velocity, scalar, series');
netcdf.putAtt(ncid, v_varid, 'standard_name', 'sea_water_y_velocity_at_v_location');
netcdf.putAtt(ncid, v_varid, 'coordinates', 'latitude longitude');

w_varid = netcdf.inqVarID(ncid, 'w');
netcdf.putAtt(ncid, w_varid, 'long_name', 'vertical momentum component');
netcdf.putAtt(ncid, w_varid, 'units', 'meter second-1');
netcdf.putAtt(ncid, w_varid, 'field', 'w-velocity, scalar, series');
netcdf.putAtt(ncid, w_varid, 'standard_name', 'upward_sea_water_velocity');
netcdf.putAtt(ncid, w_varid, 'coordinates', 'latitude longitude');

rho_varid = netcdf.inqVarID(ncid, 'rho');
netcdf.putAtt(ncid, rho_varid, 'long_name', 'densitu anomaly');
netcdf.putAtt(ncid, rho_varid, 'units', 'kilogram meter-3');
netcdf.putAtt(ncid, rho_varid, 'field', 'density, scalar, series');
netcdf.putAtt(ncid, rho_varid, 'standard_name', 'sea_water_sigma_t');

gridz_varid = netcdf.inqVarID(ncid, 'gridz');
netcdf.putAtt(ncid, gridz_varid, 'long_name', 'vertical coordinate');
netcdf.putAtt(ncid, gridz_varid, 'units', 'meter')


% Cerrar el archivo NetCDF
netcdf.close(ncid);

disp('Archivo NetCDF creado exitosamente ');