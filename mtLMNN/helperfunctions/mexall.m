% compiles all .c files in the current directory
%
% Kilian Q. Weinberger, 2006

dd=dir('*.c');
for i=1:length(dd)
fprintf('mex %s -lblas \n',dd(i).name);
 eval(['try;mex -lblas  ' dd(i).name '; catch; disp(lasterr); end;']);
end;

dd=dir('*.cpp');
for i=1:length(dd)
fprintf('mex %s -lblas \n',dd(i).name);
 eval(['try;mex -lblas  ' dd(i).name '; catch; disp(lasterr); end;']);
end;
