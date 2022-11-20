



%Fait par : Abdelazyz Rkhiss 


function varargout = IPproj(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @IPproj_OpeningFcn, ...
                   'gui_OutputFcn',  @IPproj_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end



% --- Executes just before IPproj is made visible.
function IPproj_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to IPproj (see VARARGIN)
ah = axes('unit', 'normalized', 'position', [0 0 1 1]);
bg = imread('kali.jpg'); imagesc(bg);
% Choose default command line output for IPproj
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% UIWAIT makes IPproj wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = IPproj_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
AproposFig=figure('Name','About', ...
   'NumberTitle','off',...
   'tag', 'About', ...
   'BusyAction','Queue','Interruptible','off', ...
   'position',[400 300 300 220],...
   'IntegerHandle', 'off', ...
   'WindowStyle','modal',...
   'Colormap', gray(256));

Std.Interruptible = 'off';
Std.BusyAction = 'queue';
Ctl = Std;
Ctl.Units = 'Pixels';
Ctl.Parent = AproposFig;

texte=Ctl;
texte.Interruptible='off';
texte.Style='text';
texte.Horiz='left';
texte.Foreground='red';
texte.FontWeight='bold';
texte.Fontsize=14;
texte.Background=[235 235 235]/255;

% %%/******* text label A Propos **************/
chaine1='Ce projet réalisé par :';
hs.hTexteApropos=uicontrol(texte,...
   'Foreground','black',...
   'Position',[10 110 220 100], ...
   'String',chaine1);
texte.Background=[235 235 235]/255;

chaine2='Abdelazyz RKHISS';
texte.Fontsize=12;
hs.hTexteApropos=uicontrol(texte,...
   'Position',[60 85 180 100], ...
   'String',chaine2);
texte.Background=[235 235 235]/255;

chaine3='Abdelazyz.Rkhiss@usmba.ac.ma';
texte.Fontsize=12;
texte.Foreground='red';
texte.Background=[235 235 235]/255;
hs.hTexteApropos=uicontrol(texte,...
   'Position',[10 65 280 100], ...
   'String',chaine3);
texte.Background=[235 235 235]/255;

chaine4='Professeur : ';
texte.Fontsize=14;

hs.hTexteApropos=uicontrol(texte,...
   'Foreground','black',...
   'Position',[10 26 280 100], ...
   'String',chaine4);

chaine5='M.Hamid Tairi';
hs.hTexteApropos=uicontrol(texte,...
   'Position',[60 5 280 100], ...
   'String',chaine5);

% --- Executes on button press in aboutt.
function aboutt_Callback(hObject, eventdata, handles)
AproposFig=figure('Name','About', ...
   'NumberTitle','off',...
   'tag', 'About', ...
   'BusyAction','Queue','Interruptible','off', ...
   'position',[400 300 300 220],...
   'IntegerHandle', 'off', ...
   'WindowStyle','modal',...
   'Colormap', gray(256));

Std.Interruptible = 'off';
Std.BusyAction = 'queue';
Ctl = Std;
Ctl.Units = 'Pixels';
Ctl.Parent = AproposFig;

texte=Ctl;
texte.Interruptible='off';
texte.Style='text';
texte.Horiz='left';
texte.Foreground='red';
texte.FontWeight='bold';
texte.Fontsize=14;
texte.Background=[235 235 235]/255;

% %%/******* text label A Propos **************/
chaine1='Ce projet réalisé par :';
hs.hTexteApropos=uicontrol(texte,...
   'Foreground','black',...
   'Position',[10 110 220 100], ...
   'String',chaine1);
texte.Background=[235 235 235]/255;

chaine2='Abdelazyz RKHISS';
texte.Fontsize=12;
hs.hTexteApropos=uicontrol(texte,...
   'Position',[60 85 180 100], ...
   'String',chaine2);
texte.Background=[235 235 235]/255;

chaine3='Abdelazyz.Rkhiss@usmba.ac.ma';
texte.Fontsize=12;
texte.Foreground='red';
texte.Background=[235 235 235]/255;
hs.hTexteApropos=uicontrol(texte,...
   'Position',[10 65 280 100], ...
   'String',chaine3);
texte.Background=[235 235 235]/255;

chaine4='Professeur : ';
texte.Fontsize=14;

hs.hTexteApropos=uicontrol(texte,...
   'Foreground','black',...
   'Position',[10 26 280 100], ...
   'String',chaine4);

chaine5='M.Hamid Tairi';
hs.hTexteApropos=uicontrol(texte,...
   'Position',[60 5 280 100], ...
   'String',chaine5);


% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% --------------------------------------------------------------------
function save_Callback(hObject, eventdata, handles)
image = handles.courant_data;
[file,path] = uiputfile('*.png','Save your image ...');
imwrite(image, sprintf('%s',path,file),'png');


% --------------------------------------------------------------------
function exit_Callback(hObject, eventdata, handles)
delete(handles.figure1)

function Open_Callback(hObject, eventdata, handles)
[file,path] = uigetfile('*.*');
handles.ima = imread(sprintf('%s',path,file));
axes(handles.affichage);
handles.courant_data = handles.ima;
subimage(handles.courant_data);

img=handles.courant_data;

[l c d]=size(img);

if(d==3)
    r=img(:,:,1);
    g=img(:,:,2);
    b=img(:,:,3);

    hr=imhist(r,256);
    hv=imhist(g,256);
    hb=imhist(b,256);
    axes(handles.histo);
    plot([hr,hv,hb]);
else
    axes(handles.histo);
    imhist(img,256);
end

handles.output = hObject;
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on button press in convertgris.
function convertgris_Callback(hObject, eventdata, handles)
axes(handles.affichage);

img=handles.courant_data;

[l,c,d]=size(img);
v=img;

if(d==3)
    v=zeros(l,c);
    for i=1 : l
        for j=1 : c
            v(i,j)=img(i,j,1)*0.333+img(i,j,2)*0.333+img(i,j,3)*0.333;
        end
    end
end

v=uint8(v);
imshow(v);

axes(handles.histo);
imhist(v,256);

handles.output = hObject;
guidata(hObject, handles);

% --- Executes on button press in lumin2.
function lumin2_Callback(hObject, eventdata, handles)
axes(handles.imagetr);

img=handles.courant_data;

[l c d]=size(img);

if(d==3)
    img=rgb2gray(img);
end

moy=mean2(img);
[m n]=size(img);
v=img;

for i=1 : m
    for j=1 : n
        pix=img(i,j)-moy;
        if(pix<0)
            pix=0;
        end
        v(i,j)=pix;
    end
end
v=uint8(v);

imshow(v);

axes(handles.histotrt);
imhist(v,256);

handles.output = hObject;
guidata(hObject, handles);


% --- Executes on button press in lumin1.
function lumin1_Callback(hObject, eventdata, handles)
axes(handles.imagetr);

img=handles.courant_data;
[l c d]=size(img);

if(d==3)
    img=rgb2gray(img);
end

moy=mean2(img);

[m n]=size(img);
v=img;

for i=1 : m
    for j=1 : n
        pix=img(i,j)+moy;
        if(pix>255)
            pix=255;
        end
            v(i,j)=pix;
    end
end

v=uint8(v);

imshow(v);

axes(handles.histotrt);
imhist(v,256);

handles.output = hObject;
guidata(hObject, handles);


% --- Executes on button press in ctr1.
function ctr1_Callback(hObject, eventdata, handles)
axes(handles.imagetr);

img=handles.courant_data;
[l c d]=size(img);

if(d==3)
    img=rgb2gray(img);
end

img=double(img);
ctr=std2(img);

maxi=max(max(img));
mini=min(min(img));

[m n]=size(img);

v=img;

for i=1 : m
    for j=1 : n
        
        pix=(255/(maxi-mini))*(img(i,j)-mini);
        
        if(pix>255)
            pix=255;
        elseif(pix<0)
            pix=0;
        end
            v(i,j)=pix;
    end
end

v=uint8(v);

imshow(v);

axes(handles.histotrt);
imhist(v,256);

handles.output = hObject;
guidata(hObject, handles);


% --- Executes on button press in negative.
function negative_Callback(hObject, eventdata, handles)
axes(handles.imagetr);

img=handles.courant_data;
[l c d]=size(img);

if(d==3)
    img=rgb2gray(img);
end

ctr=std2(img);

[m n]=size(img);

v=img;
for i=1 : m
    for j=1 : n
        v(i,j)=255-img(i,j);
    end
end
v=uint8(v);
imshow(v);

axes(handles.histotrt);
imhist(v,256);

handles.output = hObject;
guidata(hObject, handles);


% --- Executes on button press in meroir.
function meroir_Callback(hObject, eventdata, handles)
axes(handles.imagetr);

img=handles.courant_data;
[l c d]=size(img);

if(d==3)
    img=rgb2gray(img);
end

moy=mean2(img);

[m n]=size(img);
v=img;

for i=1 : m
    for j=1 : n
        v(i,j)=img(i,n-j+1);
    end
end
v=uint8(v);
imshow(v);

axes(handles.histotrt);
imhist(v,256);

handles.output = hObject;
guidata(hObject, handles);


% --- Executes on button press in gauss.
function gauss_Callback(hObject, eventdata, handles)
axes(handles.imagetr);

img=handles.courant_data;
[l c d]=size(img);

if(d==3)
    img=rgb2gray(img);
end

v=imnoise(img,'gaussian',0.01);
imshow(v);

axes(handles.histotrt);
imhist(v,256);

handles.output = hObject;
guidata(hObject, handles);


% --- Executes on button press in speckle.
function speckle_Callback(hObject, eventdata, handles)
axes(handles.imagetr);

img=handles.courant_data;
[l c d]=size(img);

if(d==3)
    img=rgb2gray(img);
end

v=imnoise(img,'speckle',0.01);
imshow(v);

axes(handles.histotrt);
imhist(v,256);

handles.output = hObject;
guidata(hObject, handles);


% --- Executes on button press in poivre.
function poivre_Callback(hObject, eventdata, handles)
axes(handles.imagetr);

img=handles.courant_data;
[l c d]=size(img);

if(d==3)
    img=rgb2gray(img);
end

v=imnoise(img,'salt & pepper',0.03);
imshow(v);

axes(handles.histotrt);
imhist(v,256);

handles.output = hObject;
guidata(hObject, handles);


% --- Executes on button press in gaussien.
function gaussien_Callback(hObject, eventdata, handles)
axes(handles.imagetr);

img=handles.courant_data;
[l c d]=size(img);

if(d==3)
    img=rgb2gray(img);
end

[m n]=size(img);
img=double(img);
v=img;
F=1/16 * [1 2 1; 2 4 2; 1 2 1];

for i=2 : m-1
    for j=2 : n-1
        M=img(i-1:i+1 , j-1:j+1);
        R=M.*F;
        S=sum(R(:));
        
         if(S>255)
            v(i,j)=255;
        elseif(S<0)
            v(i,j)=0;
         end
         
         v(i,j)=S;
    end
end
v=uint8(v);
imshow(v);

axes(handles.histotrt);
imhist(v,256);

handles.output = hObject;
guidata(hObject, handles);



% --- Executes on button press in median.
function median_Callback(hObject, eventdata, handles)
axes(handles.imagetr);

img=handles.courant_data;
[l c d]=size(img);

if(d==3)
    img=rgb2gray(img);
end

[m n]=size(img);
img=double(img);
v=img;

for i=2 : m-1
    for j=2 : n-1
        voisin=img(i-1:i+1 , j-1:j+1);
        
        listv=[voisin(1,:) voisin(2,:) voisin(3,:)];
        sort(listv);
        m=median(listv);
        v(i,j)=m;
    end
end

v=uint8(v);

imshow(v);

axes(handles.histotrt);
imhist(v,256);

handles.output = hObject;
guidata(hObject, handles);



% --- Executes on button press in kirsch.
function kirsch_Callback(hObject, eventdata, handles)
axes(handles.imagetr);

img=handles.courant_data;
[l c d]=size(img);

if(d==3)
    img=rgb2gray(img);
end

[m n]=size(img);
img=double(img);
v=img;
H=[-3 -3 -3 ; 5 0 -3 ; 5 5 -3 ];
v=img;

for i=2 : m-1
    for j=2 : n-1
        M=img(i-1:i+1 , j-1:j+1);
        R=M.*H;
        S=sum(R(:));
        
        if(S>255)
            v(i,j)=255;
        elseif(S<0)
            v(i,j)=0;
        end
        v(i,j)=S;
    end
end

v=uint8(v);
imshow(v);

axes(handles.histotrt);
imhist(v,256);

handles.output = hObject;
guidata(hObject, handles);


% --- Executes on button press in laplacien.
function laplacien_Callback(hObject, eventdata, handles)
axes(handles.imagetr);

img=handles.courant_data;
[l c d]=size(img);

if(d==3)
    img=rgb2gray(img);
end

[m n]=size(img);
img=double(img);
v=img;

H=[-1 -1 -1 ;-1 8 -1 ; -1 -1 -1];

for i=2 : m-1
    for j=2 : n-1
        
        M=img(i-1:i+1 , j-1:j+1);
        R=M.*H;
        S=sum(R(:));
        
        if(S>255)
            v(i,j)=255;
        elseif(S<0)
            v(i,j)=0;
        end
        v(i,j)=S;
    end
end

v=uint8(v);
imshow(v);

axes(handles.histotrt);
imhist(v,256);

handles.output = hObject;
guidata(hObject, handles);


%-----------------------------------------------------------menu----------------------------------------------------%
%-------------------------------------------------------------------------------------------------------------------%


% --------------------------------------------------------------------
function freq_Callback(hObject, eventdata, handles)
%%voila pass bat et haut de transformer de fourier

% --------------------------------------------------------------------
function fftbas_Callback(hObject, eventdata, handles)
img=handles.courant_data;
[l c d]=size(img);

if(d==3)
    img=rgb2gray(img);
end

img=double(img);

img=fftshift(fft2(img));
[m , n]=size(img); 
H=zeros(m,n);

v=H;

D=20;
 
 m2=round(m/2); 
 n2=round(n/2);
 H(m2-D:m2+D , n2-D:n2+D)=1;
 for i=1 : m
     for j=1 : n
         v(i,j)=img(i,j)*H(i,j);
     end
 end
 v=ifft2(v);
 v=abs(v);
 
 for i=1 : m
     for j=1 : n
        if(v(i,j)<0)
            v(i,j)=0;
        elseif(v(i,j)>255)
            v(i,j)=255;
        end
     end
 end
 
 v=uint8(v);
 
axes(handles.imagetr);
imshow(v);

axes(handles.histotrt);
imhist(v,256);

handles.output = hObject;
guidata(hObject, handles);


% --------------------------------------------------------------------
function ffthaut_Callback(hObject, eventdata, handles)
img=handles.courant_data;
[l c d]=size(img);

if(d==3)
    img=rgb2gray(img);
end

img=double(img);

img=fftshift(fft2(img));

[m , n]=size(img);

H=ones(m,n);

v=H;
D=10;

m2=round(m/2); 
n2=round(n/2);

H(m2-D:m2+D , n2-D:n2+D)=0;
 
for i=1 : m
     for j=1 : n
         
         v(i,j)=img(i,j)*H(i,j);

     end
end

v=ifft2(v);
v=abs(v);

for i=1 : m
     for j=1 : n
        if(v(i,j)<0)
            v(i,j)=0;
        elseif(v(i,j)>255)
            v(i,j)=255;
        end
     end
end
 
v=uint8(v);
axes(handles.imagetr);
imshow(v);

axes(handles.histotrt);
imhist(v,256);

handles.output = hObject;
guidata(hObject, handles);
% --------------------------------------------------------------------


% --------------------------------------------------------------------
function morpho_Callback(hObject, eventdata, handles)
% --------------------------------------------------------------------
function eros_Callback(hObject, eventdata, handles)

lementS=[0 1 0; 1 1 1; 0 1 0]; 
Img=handles.courant_data;
[l c d]=size(Img);

if(d==3)
    Img=rgb2gray(Img);
end
result = double(Img);
for i=2:l-1
    for j=2:c-1 
        M=double(Img((i-1:i+1),(j-1:j+1)));
        Ma=M-lementS;% 
        result(i,j) = min(min(Ma));
    end
end
result=uint8(result);
axes(handles.imagetr);
imshow(result);

axes(handles.histotrt);
imhist(result,256);

handles.output = hObject;
guidata(hObject, handles);


% --------------------------------------------------------------------
function delat_Callback(hObject, eventdata, handles)
lementS=[0 1 0; 1 1 1; 0 1 0]; 
Img=handles.courant_data;
[l c d]=size(Img);

if(d==3)
    Img=rgb2gray(Img);
end

result =double(Img);
for i = 2  :l-1
    for j = 2  :c-1
      M=double(Img((i-1:i+1),(j-1:j+1)));
      M=M+lementS;
      result(i,j)=max(max(M));
    end
end
result=uint8(result);

axes(handles.imagetr);
imshow(result);

axes(handles.histotrt);
imhist(result,256);

handles.output = hObject;
guidata(hObject, handles);


% --------------------------------------------------------------------
function ouver_Callback(hObject, eventdata, handles)
lementS=[0 1 0; 1 1 1; 0 1 0]; 
Img=handles.courant_data;
[l c d]=size(Img);

if(d==3)
    Img=rgb2gray(Img);
end
result =double(Img);
%erosion
for i = 2  :l-1
    for j = 2  :c-1
      M=double(Img((i-1:i+1),(j-1:j+1)));
      M=M-lementS;
      result(i,j)=min(min(M));
    end
end
real=result; 
%dilatation
for i = 2  :l-1
    for j = 2  :c-1
      M=double(result((i-1:i+1),(j-1:j+1)));
      M=M+lementS;
      real(i,j)=max(max(M));
    end
end

result=uint8(real);

axes(handles.imagetr);
imshow(result);

axes(handles.histotrt);
imhist(result,256);

handles.output = hObject;
guidata(hObject, handles);


% --------------------------------------------------------------------
function ferm_Callback(hObject, eventdata, handles)
Img=handles.courant_data;
lementS=[0 1 0; 1 1 1; 0 1 0]; 
[l c d]=size(Img);

if(d==3)
    Img=rgb2gray(Img);
end

real =double(Img);
%dilatation
for i = 2  :l-1
    for j = 2  :c-1
      M=double(Img((i-1:i+1),(j-1:j+1)));
      M=M+lementS;
      real(i,j)=max(max(M));
    end
end

result=real;
%erosion
for i = 2  :l-1
    for j = 2  :c-1
      M=double(real((i-1:i+1),(j-1:j+1)));
      M=M-lementS;
      result(i,j)=min(min(M));
    end
end

result=uint8(result);

axes(handles.imagetr);
imshow(result);

axes(handles.histotrt);
imhist(result,256);

handles.output = hObject;
guidata(hObject, handles);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%--------------------------------------------------------------------

function ptinteret_Callback(hObject, eventdata, handles)
% --------------------------------------------------------------------
function susan_Callback(hObject, eventdata, handles)

im=handles.courant_data;
image=double(im);

[l c d]=size(im);

if(d==3)
    image=rgb2gray(im);
end

[n,m]=size(image);
rayon=2;
alpha=50;
r=2;
alpha=alpha/100;
mask=zeros(2*rayon+1);
b=ones(rayon+1);
for i=1:rayon+1
    for j=1:rayon+1
        if (rayon==1)
           if(j>i)
            b(i,j)=0;
           end
         else
             if(j>i+1)
            b(i,j)=0;
         end
        end
    end
end
mask(1:rayon+1,rayon+1:2*rayon+1)=b;
mask(1:rayon+1,1:rayon+1)=rot90(b);
mask0=mask;
mask0=flipdim(mask0,1);
mask=mask0+mask;
mask(rayon+1,:)=mask(rayon+1,:)-1;
max_reponse=sum(sum(mask));
f=zeros(n,m);
for i=(rayon+1):n-rayon
    for j=(rayon+1):m-rayon
  
          image_courant=double(image(i-rayon:i+rayon,j-rayon:j+rayon));

    image_courant_mask=image_courant.*mask;

         inteniste_cental= image_courant_mask(rayon+1,rayon+1);
         s=exp(-1*(((image_courant_mask-inteniste_cental)/max_reponse).^6));
       somme=sum(sum(s));
%   si le centre du mask est un 0 il faut soustraire les zeros des filtres
                if (inteniste_cental==0)
                    somme=somme-length((find(mask==0)));
                end       
         f(i,j)=somme;           
     end
end
ff=f(rayon+1:n-(rayon+1),rayon+1:m-(rayon+1));
minf=min(min(ff));
maxf=max(max(f));
fff=f;
d=2*r+1;
temp1=round(n/d);
if (temp1-n/d)<0.5 &(temp1-n/d)>0
temp1=temp1-1;
end
temp2=round(m/d);
if (temp2-m/d)<0.5 &(temp2-m/d)>0
temp2=temp2-1;
end
fff(n:temp1*d+d,m:temp2*d+d)=0;
for i=(r+1):d:temp1*d+d
for j=(r+1):d:temp2*d+d
window=fff(i-r:i+r,j-r:j+r);
window0=window;
[xx,yy]=find(window0==0);
for k=1:length(xx)
window0(xx(k),yy(k))=max(max(window0));
end
minwindow=min(min(window0));
[y,x]=find(minwindow~=window & window<=minf+alpha*(maxf-minf) & window>0);
[u,v]=find(minwindow==window);
if length(u)>1
for l=2:length(u)
fff(i-r-1+u(l),j-r-1+v(l))=0 ;
end
end
if length(x)~=0
for l=1:length(y)
fff(i-r-1+y(l),j-r-1+x(l))=0 ;
end
end
end
end
seuil=minf+alpha*(maxf-minf);
[u,v]=find(minf<=fff & fff<=seuil );
axes(handles.imagetr);
imshow(image);
hold on
plot(v,u,'.b','MarkerSize',10)
nombre_de_point_dinteret=length(v)

axes(handles.histotrt);
imhist(image,256);


handles.output = hObject;
guidata(hObject, handles);

% --------------------------------------------------------------------
function Haris_Callback(hObject, eventdata, handles)

im=handles.courant_data;
image=double(im);

[l c d]=size(im);

if(d==3)
    img=rgb2gray(im);
end
%==========================================================================
lambda=0.04;
sigma=1; seuil=200; r=6; w=5*sigma;
[m,n]=size(img);

imd=double(img);
dx=[-1 0 1
    -2 0 2
    -1 0 1]; % deriv?e horizontale : filtre de Sobel
dy=dx'; % deriv?e verticale : filtre de Sobel

g = fspecial('gaussian',max(1,fix(w)), sigma);
Ix=conv2(imd,dx,'same');
Iy=conv2(imd,dy,'same');
Ix2=conv2(Ix.^2, g, 'same');
Iy2=conv2(Iy.^2, g, 'same');
Ixy=conv2(Ix.*Iy, g,'same');

detM=Ix2.*Iy2-Ixy.^2;
trM=Ix2+Iy2;
R=detM-lambda*trM.^2;
%==========================================================================
R1=(1000/(1+max(max(R))))*R;
%==========================================================================          
[u,v]=find(R1<=seuil);
nb=length(u);
for k=1:nb
    R1(u(k),v(k))=0;
end
R11=zeros(m+2*r,n+2*r);
R11(r+1:m+r,r+1:n+r)=R1;
[m1,n1]=size(R11);

for i=r+1:m1-r
    for j=r+1:n1-r
        fenetre=R11(i-r:i+r,j-r:j+r);
        ma=max(max(fenetre));
        if fenetre(r+1,r+1)<ma
            R11(i,j)=0;
        end
    end
end

nv=uint8(img); 
axes(handles.imagetr);
imshow(nv);

hold on
R11=R11(r+1:m+r,r+1:n+r);
[x,y]=find(R11);
nb=length(x);

plot(y,x,'.r');

axes(handles.histotrt);
imhist(image,256);


handles.output = hObject;
guidata(hObject, handles);


% --------------------------------------------------------------------
function linesdet_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function haugh_Callback(hObject, eventdata, handles)
im=handles.courant_data;
image=double(im);

[l c d]=size(im);

if(d==3)
    img=rgb2gray(im);
end

theta_maximum = 90;
rho_maximum = floor(sqrt(rows^2 + cols^2)) - 1;
theta_range = -theta_maximum:theta_maximum - 1;
rho_range = -rho_maximum:rho_maximum;

Hough = zeros(length(rho_range), length(theta_range));

wb = waitbar(0, 'Hough Transform');
I = edge(I,'canny',0.2);
for row = 1:rows
    waitbar(row/rows, wb);
    for col = 1:cols
        if I(row, col) > 0
            x = col - 1;
            y = row - 1;
            for theta_ = theta_range
                rho_ = round((x * cosd(theta_)) + (y * sind(theta_)));
                rho_index = rho_ + rho_maximum + 1;
                theta_index = theta_ + theta_maximum + 1;
                Hough(rho_index, theta_index) = Hough(rho_index, theta_index) + 1;
            end
        end
    end
end

close(wb);
P  = houghpeaks(Hough,6,'threshold',0.1);
% Find lines and plot them
lines = houghlines(I,theta_range,rho_range,P,'FillGap',1e4,'MinLength',5);
axes(handles.imagetr);
imshow(img);
hold on;
for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
   % Plot beginnings and ends of lines
   plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
   plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
end

axes(handles.histotrt);
imhist(image,256);


handles.output = hObject;
guidata(hObject, handles);