function parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Start of Phi parameters part                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi10  = 0.3301;
phi30  = 0.33;
del    = 0;
eps    = 0;
KK     = 0;
phirct = 0;
JJ     = 0;
phires = 0;
AA     = -2;
gamma  = -4;
timeend = -1*ones(8,8);
%timeendP = -1*ones(8,8);
%timeendD = -1*ones(8,8);
i = 1;
j = 1;
while (i<9)
    phi10 = 0.1*i
    while (j<10-i)
        text = 'phi0'
        phi30 = 0.1*j
        simul = strcat('phi0-',sprintf('%1i',i),'-',sprintf('%1i',j));
        timeend(i,j) = NewGypsum4(phi10,phi30,del,eps,KK,phirct,JJ,phires,AA,gamma,simul);
%        timeendP(i,j) = NewGypsum4P(phi10,phi30,del,eps,KK,phirct,JJ,phires,AA,gamma,simul);
%        timeendD(i,j) = NewGypsum4D(phi10,phi30,del,eps,KK,phirct,JJ,phires,AA,gamma,simul);
        j = j+1;
    end
    j = 1;
    i = i+1;
end
timeend(8,8) = NewGypsum4(1/3,1/3,del,eps,KK,phirct,JJ,phires,AA,gamma,'phi0-3rd-3rd');
%timeendP(8,8) = NewGypsum4P(1/3,1/3,del,eps,KK,phirct,JJ,phires,AA,gamma,'phi0-3rd-3rd');
%timeendD(8,8) = NewGypsum4D(1/3,1/3,del,eps,KK,phirct,JJ,phires,AA,gamma,'phi0-3rd-3rd');
filepath = 'D:\Simulations\Existence-simulations\Standard\phi0-timeend.txt'; %specific directory
%filepathP = 'D:\Simulations\3rd-set\Pressure\phi0-timeend.txt'; %specific directory
%filepathD = 'D:\Simulations\3rd-set\Drag\phi0-timeend.txt'; %specific directory
fileName = fopen(char(filepath),'w'); %create txt file at specific location
%fileNameP = fopen(char(filepathP),'w'); %create txt file at specific location
%fileNameD = fopen(char(filepathD),'w'); %create txt file at specific location
fmt1 = ['%8s\t|\t',repmat('% 1.1f\t',1,size(timeend,1)),'\r\n\v\r\n'];
fmt3 = ['%1.1f\t|\t',repmat('% 3i\t',1,size(timeend,1)),'\r\n'];
fprintf(fileName,fmt1,char('phi10'),0.1*[1:8]'); %fill the txt file with the correct variable
fprintf(fileName,fmt3,[0.1*[1:8];timeend']); %fill the txt file with the correct variable
%fprintf(fileNameP,fmt1,char('phi10'),0.1*[1:8]'); %fill the txt file with the correct variable
%fprintf(fileNameP,fmt3,[0.1*[1:8];timeendP']); %fill the txt file with the correct variable
%fprintf(fileNameD,fmt1,char('phi10'),0.1*[1:8]'); %fill the txt file with the correct variable
%fprintf(fileNameD,fmt3,[0.1*[1:8];timeendD']); %fill the txt file with the correct variable
fclose('all');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End of Phi parameters part                                               %
%Start of delta parameters part                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% phi10  = 0.2;
% phi30  = 0.5;
% del    = 0;
% eps    = 0;
% KK     = 0;
% phirct = 0;
% JJ     = 0;
% phires = 0;
% AA     = 0;
% gamma  = 0;
% timeend = -1*ones(11,1);
% timeendP = -1*ones(11,1);
% timeendD = -1*ones(11,1);
% i = 5;
% while (i<6)
%     text = 'del'
%     simul = strcat('del-',sprintf('%2i',i));
%     %timeend(i+6,1) = NewGypsum4(phi10,phi30,i,eps,KK,phirct,JJ,phires,AA,gamma,simul);
%     %timeendP(i+6,1) = NewGypsum4P(phi10,phi30,i,eps,KK,phirct,JJ,phires,AA,gamma,simul);
%     %timeendD(i+6,1) = NewGypsum4D(phi10,phi30,i,eps,KK,phirct,JJ,phires,AA,gamma,simul);
%     timeendD(i+6,1) = NewGypsum4D(phi10,phi30,0,eps,KK,phirct,JJ,phires,AA,gamma,simul);
%     i = i+1
% end
% %filepath = 'D:\Simulations\3rd-set\Standard\del-timeend.txt'; %specific directory
% %filepathP = 'D:\Simulations\3rd-set\Pressure\del-timeend.txt'; %specific directory
% %filepathD = 'D:\Simulations\3rd-set\Drag\del-timeend.txt'; %specific directory
% filepathD = 'D:\Simulations\3rd-set\Kminus\del-timeend.txt'; %specific directory
% %fileName = fopen(char(filepath),'w'); %create txt file at specific location
% %fileNameP = fopen(char(filepathP),'w'); %create txt file at specific location
% fileNameD = fopen(char(filepathD),'w'); %create txt file at specific location
% fmt1 = ['%8s\t|\t','%8s\t','\r\n\v\r\n'];
% fmt3 = ['%2.2g\t|\t','% 3i\t','\r\n'];
% %fprintf(fileName,fmt1,char('del'),char('time')); %fill the txt file with the correct variable
% %fprintf(fileName,fmt3,[10.^([-5:5]);timeend']); %fill the txt file with the correct variable
% %fprintf(fileNameP,fmt1,char('del'),char('time')); %fill the txt file with the correct variable
% %fprintf(fileNameP,fmt3,[10.^([-5:5]);timeendP']); %fill the txt file with the correct variable
% fprintf(fileNameD,fmt1,char('del'),char('time')); %fill the txt file with the correct variable
% fprintf(fileNameD,fmt3,[10.^([-5:5]);timeendD']); %fill the txt file with the correct variable
% fclose('all');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %End of delta parameters part                                             %
% %Start of epsilon parameters part                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 phi10  = 0.3;
 phi30  = 0.4;
 del    = 0;
 eps    = 0;
 KK     = 0;
 phirct = 0;
 JJ     = 0;
 phires = 0;
 AA     = -2;
 gamma  = -4;
 timeend = -1*ones(11,1);
% timeendP = -1*ones(11,1);
% timeendD = -1*ones(11,1);
 i = -4;
 while (i<7)
     text = 'eps'
     simul = strcat('eps-',sprintf('%2i',i));
     timeend(i+5,1) = NewGypsum4(phi10,phi30,del,i,KK,phirct,JJ,phires,AA,gamma,simul);
%     timeendP(i+5,1) = NewGypsum4P(phi10,phi30,del,i,KK,phirct,JJ,phires,AA,gamma,simul);
%     timeendD(i+5,1) = NewGypsum4D(phi10,phi30,del,i,KK,phirct,JJ,phires,AA,gamma,simul);
     i = i+1
 end
 filepath = 'D:\Simulations\Existence-simulations\Standard\eps-timeend.txt'; %specific directory
% filepathP = 'D:\Simulations\3rd-set\Pressure\eps-timeend.txt'; %specific directory
% filepathD = 'D:\Simulations\3rd-set\Drag\eps-timeend.txt'; %specific directory
 fileName = fopen(char(filepath),'w'); %create txt file at specific location
% fileNameP = fopen(char(filepathP),'w'); %create txt file at specific location
% fileNameD = fopen(char(filepathD),'w'); %create txt file at specific location
 fmt1 = ['%8s\t|\t','%8s\t','\r\n\v\r\n'];
 fmt3 = ['%2.2g\t|\t','% 3i\t','\r\n'];
 fprintf(fileName,fmt1,char('eps'),char('time')); %fill the txt file with the correct variable
 fprintf(fileName,fmt3,[1.4*10.^(0.5*[-10:0]);timeend']); %fill the txt file with the correct variable
% fprintf(fileNameP,fmt1,char('eps'),char('time')); %fill the txt file with the correct variable
% fprintf(fileNameP,fmt3,[1.4*10.^(0.5*[-10:0]);timeendP']); %fill the txt file with the correct variable
% fprintf(fileNameD,fmt1,char('eps'),char('time')); %fill the txt file with the correct variable
% fprintf(fileNameD,fmt3,[1.4*10.^(0.5*[-10:0]);timeendD']); %fill the txt file with the correct variable
 fclose('all');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %End of epsilon parameters part                                           %
% %Start of Kphireact parameters part                                       %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% phi10  = 0.2;
% phi30  = 0.5;
% del    = 0;
% eps    = 0;
% KK     = 0;
% phirct = 0;
% JJ     = 0;
% phires = 0;
% AA     = 0;
% gamma  = 0;
% timeend = -1*ones(11,11);
% timeendP = -1*ones(11,11);
% timeendD = -1*ones(11,11);
% i = -5;
% j = 0;
% %i=1;
% %j=1;
% while (i<6)
%     while (j<11)
%         text = 'Kphireact'
%         simul = strcat('Kphireact-',sprintf('%1i',i),'-',sprintf('%1i',j));
%         timeend(i+6,j+1) = NewGypsum4(phi10,phi30,del,eps,i,j,JJ,phires,AA,gamma,simul);
%         timeendP(i+6,j+1) = NewGypsum4P(phi10,phi30,del,eps,i,j,JJ,phires,AA,gamma,simul);
%         timeendD(i+6,j+1) = NewGypsum4D(phi10,phi30,del,eps,i,j,JJ,phires,AA,gamma,simul);
%         j = j+1
%     end
%     j = 0;
%     i = i+1
% end
% filepath = 'D:\Simulations\3rd-set\Standard\Kphireact-timeend.txt'; %specific directory
% filepathP = 'D:\Simulations\3rd-set\Pressure\Kphireact-timeend.txt'; %specific directory
% filepathD = 'D:\Simulations\3rd-set\Drag\Kphireact-timeend.txt'; %specific directory
% fileName = fopen(char(filepath),'w'); %create txt file at specific location
% fileNameP = fopen(char(filepathP),'w'); %create txt file at specific location
% fileNameD = fopen(char(filepathD),'w'); %create txt file at specific location
% fmt1 = ['%8s\t|\t',repmat('% 2.2g\t',1,size(timeend,1)),'\r\n\v\r\n'];
% fmt3 = ['%2.2g\t|\t',repmat('% 3i\t',1,size(timeend,1)),'\r\n'];
% fprintf(fileName,fmt1,char('Kappa'),10.^(-0.5*[0:10]')); %fill the txt file with the correct variable
% fprintf(fileName,fmt3,[23*10.^([-5:5]);timeend']); %fill the txt file with the correct variable
% fprintf(fileNameP,fmt1,char('Kappa'),10.^(-0.5*[0:10]')); %fill the txt file with the correct variable
% fprintf(fileNameP,fmt3,[23*10.^([-5:5]);timeendP']); %fill the txt file with the correct variable
% fprintf(fileNameD,fmt1,char('Kappa'),10.^(-0.5*[0:10]')); %fill the txt file with the correct variable
% fprintf(fileNameD,fmt3,[23*10.^([-5:5]);timeendD']); %fill the txt file with the correct variable
% fclose('all');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %End of Kphireact parameters part                                         %
% %Start of Jphires parameters part                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% phi10  = 0.2;
% phi30  = 0.5;
% del    = 0;
% eps    = 0;
% KK     = 0;
% phirct = 0;
% JJ     = 0;
% phires = 0;
% AA     = 0;
% gamma  = 0;
% timeend = -1*ones(11,11);
% timeendP = -1*ones(11,11);
% timeendD = -1*ones(11,11);
% i = -5;
% j = 0;
% while (i<6)
%     while (j<11)
%         text = 'Jphires'
%         simul = strcat('Jphires-',sprintf('%1i',i),'-',sprintf('%1i',j));
%         timeend(i+6,j+1) = NewGypsum4(phi10,phi30,del,eps,KK,phirct,i,j,AA,gamma,simul);
%         timeendP(i+6,j+1) = NewGypsum4P(phi10,phi30,del,eps,KK,phirct,i,j,AA,gamma,simul);
%         timeendD(i+6,j+1) = NewGypsum4D(phi10,phi30,del,eps,KK,phirct,i,j,AA,gamma,simul);
%         j = j+1
%     end
%     j = 0;
%     i = i+1
% end
% filepath = 'D:\Simulations\3rd-set\Standard\Jphires-timeend.txt'; %specific directory
% filepathP = 'D:\Simulations\3rd-set\Pressure\Jphires-timeend.txt'; %specific directory
% filepathD = 'D:\Simulations\3rd-set\Drag\Jphires-timeend.txt'; %specific directory
% fileName = fopen(char(filepath),'w'); %create txt file at specific location
% fileNameP = fopen(char(filepathP),'w'); %create txt file at specific location
% fileNameD = fopen(char(filepathD),'w'); %create txt file at specific location
% fmt1 = ['%8s\t|\t',repmat('% 2.2g\t',1,size(timeend,1)),'\r\n\v\r\n'];
% fmt3 = ['%2.2g\t|\t',repmat('% 3i\t',1,size(timeend,1)),'\r\n'];
% fprintf(fileName,fmt1,char('J'),10.^(-0.5*[0:10]')); %fill the txt file with the correct variable
% fprintf(fileName,fmt3,[0.4*10.^([-5:5]);timeend']); %fill the txt file with the correct variable
% fprintf(fileNameP,fmt1,char('J'),10.^(-0.5*[0:10]')); %fill the txt file with the correct variable
% fprintf(fileNameP,fmt3,[0.4*10.^([-5:5]);timeendP']); %fill the txt file with the correct variable
% fprintf(fileNameD,fmt1,char('J'),10.^(-0.5*[0:10]')); %fill the txt file with the correct variable
% fprintf(fileNameD,fmt3,[0.4*10.^([-5:5]);timeendD']); %fill the txt file with the correct variable
% fclose('all');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %End of Jphires parameters part                                           %
% %Start of Agamma parameters part                                          %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% phi10  = 0.3;
% phi30  = 0.4;
% del    = 0;
% eps    = 0;
% KK     = 0;
% phirct = 0;
% JJ     = 0;
% phires = 0;
% AA     = 0;
% gamma  = -4;
% timeend = -1*ones(8,8);
% timeendP = -1*ones(8,8);
% timeendD = -1*ones(8,8);
% i = -5;
% j = -5;
% while (i<3)
%     while (j<3)
%         text = 'Agamma'
%         simul = strcat('Agamma-',sprintf('%1i',i),'-',sprintf('%1i',j));
%         timeend(i+6,j+6) = NewGypsum4(phi10,phi30,del,eps,KK,phirct,JJ,phires,i,j,simul);
%         timeendP(i+6,j+6) = NewGypsum4P(phi10,phi30,del,eps,KK,phirct,JJ,phires,i,j,simul);
%         timeendD(i+6,j+6) = NewGypsum4D(phi10,phi30,del,eps,KK,phirct,JJ,phires,i,j,simul);
%         j = j+1
%     end
%     j = -5;
%     i = i+1
% end
% filepath = 'D:\Simulations\Existence-simulations\Standard\Agamma-timeend.txt'; %specific directory
% filepathP = 'D:\Simulations\3rd-set\Pressure\Agamma-timeend.txt'; %specific directory
% filepathD = 'D:\Simulations\3rd-set\Drag\Agamma-timeend.txt'; %specific directory
% fileName = fopen(char(filepath),'w'); %create txt file at specific location
% fileNameP = fopen(char(filepathP),'w'); %create txt file at specific location
% fileNameD = fopen(char(filepathD),'w'); %create txt file at specific location
% fmt1 = ['%8s\t|\t',repmat('% 2.6g\t',1,size(timeend,1)),'\r\n\v\r\n'];
% fmt3 = ['%2.6g\t|\t',repmat('% 3i\t',1,size(timeend,1)),'\r\n'];
% fprintf(fileName,fmt1,char('A'),10.^(0.5*[9:16]')); %fill the txt file with the correct variable
% fprintf(fileName,fmt3,[0.379+0.003*([0:7]);timeend']); %fill the txt file with the correct variable
% fprintf(fileNameP,fmt1,char('A'),0.5*10.^([-5:2]')); %fill the txt file with the correct variable
% fprintf(fileNameP,fmt3,[0.5*10.^([-5:2]);timeendP']); %fill the txt file with the correct variable
% fprintf(fileNameD,fmt1,char('A'),0.5*10.^([-5:2]')); %fill the txt file with the correct variable
% fprintf(fileNameD,fmt3,[0.5*10.^([-5:2]);timeendD']); %fill the txt file with the correct variable
% fclose('all');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End of Agamma parameters part                                            %
%Start of delta t parameters part                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% phi10  = 0.2;
% phi30  = 0.5;
% del    = 0;
% eps    = 0;
% KK     = 0;
% phirct = 0;
% JJ     = 0;
% phires = 0;
% AA     = 0;
% gamma  = 0;
% timeend = -1*ones(11,1);
% timeendP = -1*ones(11,1);
% timeendD = -1*ones(11,1);
% i = -5;
% while (i<6)
%     text = 'deltime'
%     simul = strcat('deltime-',sprintf('%2i',i));
%     timeend(i+6,1) = NewGypsum4(phi10,phi30,i,eps,KK,phirct,JJ,phires,AA,gamma,simul);
%     timeendP(i+6,1) = NewGypsum4P(phi10,phi30,i,eps,KK,phirct,JJ,phires,AA,gamma,simul);
%     timeendD(i+6,1) = NewGypsum4D(phi10,phi30,i,eps,KK,phirct,JJ,phires,AA,gamma,simul);
%     timeendD(i+6,1) = NewGypsum4D(phi10,phi30,del,eps,KK,phirct,JJ,phires,AA,gamma,simul,i);
%     i = i+1
% end
% filepath = 'D:\Simulations\3rd-set\Standard\del-timeend.txt'; %specific directory
% filepathP = 'D:\Simulations\3rd-set\Pressure\del-timeend.txt'; %specific directory
% filepathD = 'D:\Simulations\3rd-set\Drag\del-timeend.txt'; %specific directory
% filepathD = 'D:\Simulations\3rd-set\Kminus\deltime-timeend.txt'; %specific directory
% fileName = fopen(char(filepath),'w'); %create txt file at specific location
% fileNameP = fopen(char(filepathP),'w'); %create txt file at specific location
% fileNameD = fopen(char(filepathD),'w'); %create txt file at specific location
% fmt1 = ['%8s\t|\t','%8s\t','\r\n\v\r\n'];
% fmt3 = ['%2.2g\t|\t','% 3i\t','\r\n'];
% fprintf(fileName,fmt1,char('del'),char('time')); %fill the txt file with the correct variable
% fprintf(fileName,fmt3,[10.^([-5:5]);timeend']); %fill the txt file with the correct variable
% fprintf(fileNameP,fmt1,char('del'),char('time')); %fill the txt file with the correct variable
% fprintf(fileNameP,fmt3,[10.^([-5:5]);timeendP']); %fill the txt file with the correct variable
% fprintf(fileNameD,fmt1,char('deltime'),char('time')); %fill the txt file with the correct variable
% fprintf(fileNameD,fmt3,[10.^([-5:5]);timeendD']); %fill the txt file with the correct variable
% fclose('all');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %End of delta t parameters part                                           %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return
end