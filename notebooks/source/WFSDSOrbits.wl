(* ::Package:: *)

BeginPackage["WFSdSOrbits`"];

sol::usage="sol[p0,e0,y] computes the evolution of osculating elements (p,e,f,\[Delta]) for initial p0 and e0, for a given SdS parameter y";
wfsdsorb::usage= "wfsdsorb[p0,e0,y0,tf0] plots the SdS orbit (for a given y0) vs a reference 2.5PN y=0 orbit, with fixed (p0, e0) up to time tf0";
wfsdsgwplot::usage="wfsdsgwplot[p0_,e0_,y0_,tf0_] plots the plus and cross polarization against each other";
wfsdsgwplusplot::usage="wfsdsgwplusplot[p0_,e0_,y0_,tf0_] plots the plus polarization against a reference 0pN waveform ";
wfsdsgwcrossplot::usage="wfsdsgwcrossplot[p0_,e0_,y0_,tf0_] plots the cross polarization against a reference 0pN waveform";
p::usage="semilatus rectum evolution"
e::usage="eccentricity evolution"
f::usage="true anomaly evolution"
\[Delta]::usage="periapsis evolution"
t::usage="time input in functions"
xpos::usage="xpos computes x position of emri at time t"
ypos::usage="ypos computes y position of emri at time t"
hplus::usage="hplus[x0] analytically derives the pN plus polarization at time x0 in terms of rates of orbital parameter evolution "
hcross::usage="hcross[x0] analytically derives the pN cross polarization at time x0 in terms of rates of orbital parameter evolution "
hplusSdS::usage="hplusSdS[y,x0] analytically derives the SdS (for a given y) plus polarization at time x0 in terms of rates of orbital parameter evolution "
hcrossSdS::usage="hcross[y,x0] analytically derives the SdS (for a given y) cross polarization at time x0 in terms of rates of orbital parameter evolution "

Begin["`Private`"];


(*Osculating parameter evolution*)

fdot[p_,e_,f_,y_]:=(Sqrt[1/p^3]-y/(2 Sqrt[1/p^3]))(1+e Cos[f])^2; (*frequency evolution*)
deldot[p_,e_,f_,y_]:=(Sqrt[1/p^3]-y/(2 Sqrt[1/p^3])) 3/p; (*1 pn periapsis evolution*)
pdot[p_,e_,y_]:=-((8 (1-e^2)^(3/2) (8+7 e^2))/(5 p^3))+(16 (-1+e^2)+(4 (-168-283 e^2+26 e^4))/(15 (1-e^2)^(3/2))) y; (*semilatus rectum evolution*)
edot[p_,e_,y_]:=-((e (1-e^2)^(3/2) (304+121 e^2))/(15 p^4))+((2 (-24-965 e^2-784 e^4+73 e^6))/(15 e (1-e^2)^(3/2) p)+(32 (1-e^2)^(3/2) (-1+Sqrt[1-e^2]))/(3 e p)) y;(*eccentricity evolution*)

(*x and y position of particle*)

xpos[t_]:=p[t]/(1+e[t]Cos[f[t]]) Cos[f[t]+\[Delta][t]];
ypos[t_]:=p[t]/(1+e[t]Cos[f[t]]) Sin[f[t]+\[Delta][t]];

(*solver of dissipative and conservative components*)

sol[p0_(*initial semi latus rectum*),e0_(*initial ecentricity*),y_(*sds parameter*)]:=Module[{pInit,eInit,yInit},
pInit=p0;
eInit=e0;
yInit=y;
NDSolve[{p'[t]==pdot[p[t],e[t],y],e'[t]==edot[p[t],e[t],y],f'[t]==fdot[p[t],e[t],f[t],y],\[Delta]'[t]==deldot[p[t],e[t],f[t],y],f[0]==0,e[0]==eInit,p[0]==pInit,\[Delta][0]==0, WhenEvent[p[t]==6+2e[t],"StopIntegration"](*stop at near separatrix*)},{p,e,f,\[Delta]},{t,0,10^6(*cycles*)}][[1]]];

(*hplus and hcross*)

hplus[x0_]:=Evaluate[D[xpos[t]^2-ypos[t]^2,{t,2}]]/.t->x0
hcross[x0_]:=Evaluate[D[2xpos[t]ypos[t],{t,2}]]/.t->x0
hplusSdS[y_,x0_]:=Evaluate[D[xpos[t]^2-ypos[t]^2,{t,2}]]+2Sqrt[y]D[xpos[t]^2-ypos[t]^2,{t,1}]/.t->x0
hcrossSdS[y_,x0_]:=Evaluate[D[2xpos[t]ypos[t],{t,2}]]+2 Sqrt[y]/3 D[2xpos[t]ypos[t],{t,1}]/.t->x0


(*plotter*)

Get["MaTeX`"];(*MaTeX for plot generation*)
SetOptions[MaTeX,FontSize->19,"Preamble"->{"\\usepackage{newtxtext,newtxmath}"},Magnification->1.2];
style=Directive[FontFamily->"Times",FontSize->23,Black];

wfsdsorb[p0_,e0_,y0_,tf0_]:=Module[{tf,p,e,y},
p=p0;
e=p0;
y=y0;
tf=tf0;
ParametricPlot[{Evaluate[{xpos[t],ypos[t]}/.sol[p0,e0,0]],Evaluate[{xpos[t],ypos[t]}/.sol[p0,e0,y]]},{t,0,tf},
Frame->True,ImageSize->600,
PlotStyle->{{Dashed,Thick,RGBColor[0.06274509803921569, 0.25882352941176473`, 0.5137254901960784]},{Thick,RGBColor[0.7450980392156863, 0.06666666666666667, 0.00392156862745098]}},
LabelStyle->style,FrameLabel->{MaTeX["r(\\psi) \\cos(\\phi)"],MaTeX["r(\\psi) \\sin(\\phi)"]},
Epilog->{{PointSize->0.015,RGBColor[0.06274509803921569, 0.25882352941176473`, 0.5137254901960784],
Point[{Evaluate[xpos[tf0]/.sol[p0,e0,0]],Evaluate[ypos[tf0]/.sol[p0,e0,0]]}]},{PointSize->0.015,RGBColor[0.7450980392156863, 0.06666666666666667, 0.00392156862745098],Point[{Evaluate[xpos[tf0]/.sol[p0,e0,y0]],Evaluate[ypos[tf0]/.sol[p0,e0,y0]]}]}},
PlotLegends->Placed[{MaTeX["\\lambda=0"],MaTeX["\\lambda="]TraditionalForm[y0//N]},{Left,Bottom}]]];


wfsdsgwplot[p0_,e0_,y0_,tf0_]:=Module[{tf,p,e,y},
p=p0;
e=p0;
y=y0;
tf=tf0;
Plot[{Evaluate[hplus[t]/.sol[p0,e0,0]],Evaluate[hcross[t]/.sol[p0,e0,0]]},{t,0,tf0},
Frame->True,AspectRatio->1/4,ImageSize->1000,PlotRange->{Automatic,Full},
LabelStyle->style,FrameLabel->{MaTeX["t"],MaTeX["h"]},
PlotLegends->Placed[{MaTeX["h_+"],MaTeX["h_\\times"]},{Left,Bottom}]]];

wfsdsgwplusplot[p0_,e0_,y0_,tf0_]:=Module[{tf,p,e,y},
p=p0;
e=p0;
y=y0;
tf=tf0;
Plot[{Evaluate[hplus[t]/.sol[p0,e0,0]],Evaluate[hplusSdS[y,t]/.sol[p0,e0,y]]},{t,0,tf0},
Frame->True,AspectRatio->1/4,ImageSize->1000,PlotRange->{Automatic,Full},
LabelStyle->style,FrameLabel->{MaTeX["t"],MaTeX["h_+"]},
PlotStyle->{{Dashed,Thick,RGBColor[0.06274509803921569, 0.25882352941176473`, 0.5137254901960784]},
{Thick,RGBColor[0.7450980392156863, 0.06666666666666667, 0.00392156862745098]}},
PlotLegends->Placed[{MaTeX["\\lambda=0"],MaTeX["\\lambda="]TraditionalForm[y0//N]},{Left,Bottom}]]];

wfsdsgwcrossplot[p0_,e0_,y0_,tf0_]:=Module[{tf,p,e,y},
p=p0;
e=p0;
y=y0;
tf=tf0;
Plot[{Evaluate[hcross[t]/.sol[p0,e0,0]],Evaluate[hcrossSdS[y,t]/.sol[p0,e0,y]]},{t,0,tf0},
Frame->True,AspectRatio->1/4,ImageSize->1000,PlotRange->{Automatic,Full},LabelStyle->style,
FrameLabel->{MaTeX["t"],MaTeX["h_\\times"]},
PlotStyle->{{Dashed,Thick,RGBColor[0.06274509803921569, 0.25882352941176473`, 0.5137254901960784]},
{Thick,RGBColor[0.7450980392156863, 0.06666666666666667, 0.00392156862745098]}},
PlotLegends->Placed[{MaTeX["\\lambda=0"],MaTeX["\\lambda="]TraditionalForm[y0//N]},{Left,Bottom}]]];


End[];
EndPackage[];

