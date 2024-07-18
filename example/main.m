% Dieno Diba 2022
% 2D finite element magnetotelluric inversion
% Run the inversion

clear
clc
close all

input_m0 = 'input_m0.txt';
input_dat = 'input_data.txt';
input_stg = 'input_setting.txt';
input_topo = 'input_topo.txt';

Inversion(input_m0,input_dat,input_topo,input_stg)
