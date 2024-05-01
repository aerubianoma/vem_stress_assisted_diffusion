%% Add_local_path

% Include all folders under the current path and the path is no longer valid after restarting MATLAB.

clc;
clear;

close all
warning off

cd ..

addpath(genpath(pwd)); 

% Change interpreter to latex for the plots

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

cd bin
