function varargout = Termo_interfaz(varargin)
% TERMO_INTERFAZ MATLAB code for Termo_interfaz.fig
%      TERMO_INTERFAZ, by itself, creates a new TERMO_INTERFAZ or raises the existing
%      singleton*.
%
%      H = TERMO_INTERFAZ returns the handle to a new TERMO_INTERFAZ or the handle to
%      the existing singleton*.
%
%      TERMO_INTERFAZ('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TERMO_INTERFAZ.M with the given input arguments.
%
%      TERMO_INTERFAZ('Property','Value',...) creates a new TERMO_INTERFAZ or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Termo_interfaz_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Termo_interfaz_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Termo_interfaz

% Last Modified by GUIDE v2.5 09-Dec-2020 08:56:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Termo_interfaz_OpeningFcn, ...
                   'gui_OutputFcn',  @Termo_interfaz_OutputFcn, ...
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
% End initialization code - DO NOT EDIT


% --- Executes just before Termo_interfaz is made visible.
function Termo_interfaz_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Termo_interfaz (see VARARGIN)

% Choose default command line output for Termo_interfaz
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Termo_interfaz wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Termo_interfaz_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double
Tc=str2num(get(handles.edit1,'string')) ;
Pc=str2num(get(handles.edit2,'string')) ;
w=str2num(get(handles.edit3,'string')) ;
Cp=str2num(get(handles.edit4,'string')) ;

T1=str2num(get(handles.edit5,'string')) ;
P1=str2num(get(handles.edit6,'string')) ;
T2=str2num(get(handles.edit7,'string')) ;
P2=str2num(get(handles.edit8,'string')) ;

Tipo1=get(handles.radiobutton1,'Value') ;
Tipo2=get(handles.radiobutton2,'Value') ;

if Tipo1
    Tipo=Tipo1 ;
else
    Tipo=Tipo2 ;
end

A=str2num(get(handles.edit13,'string')) ;
B=str2num(get(handles.edit14,'string')) ;
C=str2num(get(handles.edit15,'string')) ;
F=str2num(get(handles.edit16,'string')) ;

Tc=Tc+273.15;

UnidadTC=get(handles.radiobutton5,'Value') ;
UnidadTK=get(handles.radiobutton6,'Value') ;

if UnidadTK
T1=T1+273.15;
T2=T2+273.15;
end

Pvap1=Pvapcam(T1,A,B,C,Tipo,F) ;
Pvap2=Pvapcam(T2,A,B,C,Tipo,F) ;

if UnidadTC
T1=T1+273.15;
T2=T2+273.15;
end


T1
T2

Tr1=T1/Tc;
Tr2=T2/Tc;
Pr1=P1/Pc;
Pr2=P2/Pc;
m=0.48+1.574*w-0.176*w^2;
alpha1=(1+m*(1-sqrt(Tr1)))^2;
alpha2=(1+m*(1-sqrt(Tr2)))^2;
A1=0.42748*(alpha1*Pr1/Tr1^2);
A2=0.42748*(alpha2*Pr2/Tr2^2);
B1=0.08664*(Pr1/Tr1);
B2=0.08664*(Pr2/Tr2);


f1=@(x) (x.^3-x.^2+x*(A1-B1-B1^2)-A1*B1);
f2=@(x) (x.^3-x.^2+x*(A2-B2-B2^2)-A2*B2);
syms x
df1=inline(diff(f1(x)));
df2=inline(diff(f2(x)));

ZL1=NR(f1,df1,0);
ZG1=NR(f1,df1,1);
ZL2=NR(f2,df2,0);
ZG2=NR(f2,df2,1);

if P1<=Pvap1+Pvap1/5 && P1>=Pvap1-Pvap1/5
 set(handles.text24,'Visible','On')
 set(handles.edit19,'Visible','On')
 else
 set(handles.text24,'Visible','Off')
 set(handles.edit19,'Visible','Off')
end

if P2<=Pvap2+Pvap2/5 && P2>=Pvap2-Pvap2/5
 set(handles.text25,'Visible','On')
 set(handles.edit20,'Visible','On')
else
 set(handles.text25,'Visible','Off')
 set(handles.edit20,'Visible','Off')
end


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
R=8.3145;
format short

Tc=str2num(get(handles.edit1,'string')) ;
Pc=str2num(get(handles.edit2,'string')) ;
w=str2num(get(handles.edit3,'string')) ;
Cp=str2num(get(handles.edit4,'string')) ;

T1=str2num(get(handles.edit5,'string')) ;
P1=str2num(get(handles.edit6,'string')) ;
T2=str2num(get(handles.edit7,'string')) ;
P2=str2num(get(handles.edit8,'string')) ;

Tipo1=get(handles.radiobutton1,'Value') ;
Tipo2=get(handles.radiobutton2,'Value') ;

if Tipo1
    Tipo=Tipo1 ;
else
    Tipo=Tipo2 ;
end

A=str2num(get(handles.edit13,'string')) ;
B=str2num(get(handles.edit14,'string')) ;
C=str2num(get(handles.edit15,'string')) ;
F=str2num(get(handles.edit16,'string')) ;



UnidadTC=get(handles.radiobutton5,'Value') ;
UnidadTK=get(handles.radiobutton6,'Value') ;
Tc=Tc+273.15;

if UnidadTK
T1=T1+273.15;
T2=T2+273.15;
end

Pvap1=Pvapcam(T1,A,B,C,Tipo,F) ;
Pvap2=Pvapcam(T2,A,B,C,Tipo,F) ;

if UnidadTC
T1=T1+273.15;
T2=T2+273.15;
end

Tr1=T1/Tc;
Tr2=T2/Tc;
Pr1=P1/Pc;
Pr2=P2/Pc;
m=0.48+1.574*w-0.176*w^2;
alpha1=(1+m*(1-sqrt(Tr1)))^2;
alpha2=(1+m*(1-sqrt(Tr2)))^2;
A1=0.42748*(alpha1*Pr1/Tr1^2);
A2=0.42748*(alpha2*Pr2/Tr2^2);
B1=0.08664*(Pr1/Tr1); 
B2=0.08664*(Pr2/Tr2);
a=0.42748*(R*Tc)^2/Pc;
b=0.08664*R*Tc/Pc;

f1=@(x) (x.^3-x.^2+x*(A1-B1-B1^2)-A1*B1);
f2=@(x) (x.^3-x.^2+x*(A2-B2-B2^2)-A2*B2);
syms x
df1=inline(diff(f1(x)));
df2=inline(diff(f2(x)));

ZL1=NR(f1,df1,0);
ZG1=NR(f1,df1,1);
ZL2=NR(f2,df2,0);
ZG2=NR(f2,df2,1);

vig1=R*T1/P1;
vig2=R*T2/P2;
vL1=R*T1*ZL1/P1; vG1=R*T1*ZG1/P1;
vL2=R*T2*ZL2/P2; vG2=R*T2*ZG2/P2;
dalpha1=-m*(((1+m)/(sqrt(Tr1)*Tc))-(m/Tc));
dalpha2=-m*(((1+m)/(sqrt(Tr2)*Tc))-(m/Tc)); 

if P1>Pvap1+Pvap1/5 || P1<Pvap1-Pvap1/5
    if P1>Pvap1+Pvap1/5
    v1=R*T1*ZL1/P1;
    else
    v1=R*T1*ZG1/P1;
    end
    
ur1=-(a*(alpha1-dalpha1*T1)/b)*log((v1+b)/v1);
sr1=R*log((v1-b)/vig1)+((dalpha1*T1)/b/R)*log((v1+b)/v1); 

else
    X1=str2num(get(handles.edit19,'string')) ;
    v1=R*T1/P1*(ZL1+X1*(ZG1-ZL1));
    
ur1L=-(a*(alpha1-dalpha1*T1)/b)*log((vL1+b)/vL1);
sr1L=R*log((vL1-b)/vig1)+((dalpha1*T1)/b/R)*log((vL1+b)/vL1);     
ur1G=-(a*(alpha1-dalpha1*T1)/b)*log((vG1+b)/vG1);
sr1G=R*log((vG1-b)/vig1)+((dalpha1*T1)/b/R)*log((vG1+b)/vG1);     

ur1=ur1L+X1*(ur1G-ur1L) ;
sr1=sr1L+X1*(sr1G-sr1L) ;

end

if P2>Pvap2+Pvap2/5 || P2<Pvap2-Pvap2/5
    if P2>Pvap2+Pvap2/5
    v2=R*T2*ZL2/P2;
    else
    v2=R*T2*ZG2/P2;
    end
    
ur2=-(a*(alpha2-dalpha2*T2)/b)*log((v2+b)/v2);
sr2=R*log((v2-b)/vig2)+((dalpha2*T2)/b/R)*log((v2+b)/v2);
else
    X2=str2num(get(handles.edit20,'string')) ;
    v2=R*T2/P2*(ZL2+X2*(ZG2-ZL2));

ur2L=-(a*(alpha2-dalpha2*T2)/b)*log((vL2+b)/vL2);
sr2L=R*log((vL2-b)/vig2)+((dalpha2*T2)/b/R)*log((vL2+b)/vL2);     
ur2G=-(a*(alpha2-dalpha2*T2)/b)*log((vG2+b)/vG2);
sr2G=R*log((vG2-b)/vig2)+((dalpha2*T2)/b/R)*log((vG2+b)/vG2);     

ur2=ur2L+X2*(ur2G-ur2L) ;
sr2=sr2L+X2*(sr2G-sr2L) ;  
end

%Hallar propiedades para gas ideal%

dh_ig=Cp*(T2-T1);
ds_ig=Cp*log(T2/T1)-R*log(P2/P1);
du_ig=(Cp-R)*(T2-T1);
dg_ig=dh_ig-ds_ig*(T2-T1);
da_ig=du_ig-ds_ig*(T2-T1);

%PROPIEDADES RESIDUALES

hr1=ur1+P1*v1-R*T1;
hr2=ur2+P2*v2-R*T2;

gr1=hr1-sr1*T1;
gr2=hr2-sr2*T2;

ar1=ur1-sr1*T1;
ar2=ur2-sr2*T2;

%RESULTADOS%

dh12=dh_ig+hr2-hr1 ;
ds12=ds_ig+sr2-sr1 ;
du12=du_ig+ur2-ur1 ;
dg12=dg_ig+gr2-gr1 ;
da12=da_ig+ar2-ar1 ;

%domo

   n=30 ; T=linspace(293.15,Tc,n) ;
   
      T=T-273.15 ;
    if Tipo==1
      P=exp(A-B./(T+C)) ;
    else
      P=10.^(A-B./(T+C)) ;
    end
   T=sort(T) ; P=sort(P) ;
   Psat=P'*F ;
   T=T'+273.15 ;
   
   Tr=T/Tc ;
   Pr=Psat/Pc ;      
   alfa=(1+m*(1-sqrt(Tr))).^2 ;
   A_=0.427480.*alfa.*Pr./Tr.^2 ;
   B_=0.086640.*Pr./Tr ;
   
   for i=1:n
       
    syms x
  f=x.^3-x.^2+(A_(i)-B_(i)-(B_(i))^2).*x-A_(i)*B_(i) ;
  f=inline(f) ;
  df=inline(diff(f(x))) ;
%Iteraciones
  pl0=0;
   ZLd(i,1)=NR(f,df,pl0) ;
  pg0=1;
   ZGd(i,1)=NR(f,df,pg0) ;
%Volumenes
   vLd(i,1)=ZLd(i,1)*R*T(i,1)/Psat(i,1) ;
   vGd(i,1)=ZGd(i,1)*R*T(i,1)/Psat(i,1) ;
   
   end

Domo=table(T,Psat,Tr,Pr,A_,B_,alfa,ZLd,ZGd,vLd,vGd)   ;
vGd=flip(vGd) ;
PGd=flip(Psat) ;
vdomo=[vLd ; vGd] ;
Pdomo=[Psat; PGd] ;

axes(handles.axes1)
cla
t = [0:0.01:1] ;
plot(t,f1(t),'b',t,f2(t),'r',ZL1,f1(ZL1),'ob',...
     ZL2,f2(ZL2),'or',ZG1,f1(ZG1),'ob',ZG2,f2(ZG2),'or')
legend('f1(Z)','f2(Z)')
xlabel('Z') ; ylabel('f(A,B,Z)') ;
axis tight
grid on

axes(handles.axes2)
cla
 P11=@(v) -a*alpha1./(v.*(v+b))+R*T1./(v-b) ;
 P22=@(v) -a*alpha2./(v.*(v+b))+R*T2./(v-b) ;

 v=[vLd(1):0.001:vGd(end)] ;
semilogx(v,P11(v),'b',v,P22(v),'r',v1,P1,'ob',v2,P2,'or',vdomo,Pdomo,'k')
xlabel('v(m3/kmol)') ; ylabel('P (kPa)') ;
legend('T1','T2','v1,P1','v2,P2','Domo sat.')
xlim([vLd(1) vGd(end)]) ; ylim([0 max(Pdomo(1:n))+max(Pdomo(1:n))/10]) ;
grid on

set(handles.edit17,'String',num2str(v1))
set(handles.edit18,'String',num2str(v2))

format short
u=[du_ig ur1 ur2 du12] ;
h=[dh_ig hr1 hr2 dh12] ;
s=[ds_ig sr1 sr2 ds12] ;
g=[dg_ig gr1 gr2 dg12] ;
ai=[da_ig ar1 ar2 da12] ;

Propiedades=[u ; h; s; g; ai] ;
set(handles.uitable1,'Data',Propiedades)


function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a double


% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit18_Callback(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit18 as text
%        str2double(get(hObject,'String')) returns contents of edit18 as a double


% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit19 as text
%        str2double(get(hObject,'String')) returns contents of edit19 as a double


% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double


% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
