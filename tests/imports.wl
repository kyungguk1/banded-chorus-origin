(* ::Package:: *)

Clear[importField]
importField[file:_String?FileExistsQ,xspan:_Span:;;]:=With[{(*file=metadata[["files",1]],*)group="/field",names={"B","E"}},
Module[{attr,time,dsets},
attr=Import[file,{"HDF5","Attributes",FileNameJoin[{group,First[names]},OperatingSystem->"Unix"]}];
time=attr["time"];
dsets=Import[file,{"HDF5","Datasets",FileNameJoin[{group,#},OperatingSystem->"Unix"]&/@names}];
Association[Append[Thread[Rule[StringJoin/@Tuples[{{"d"},names,Characters["123"]}],Part[#,xspan]&/@Apply[Join,Transpose/@dsets]]],"time"->time]]
]
]


Clear[importMoment]
importMoment[file_String?FileExistsQ]:=(Print["importMoment requires the second argument."];Abort[])
importMoment[file:_String?FileExistsQ,spid:_Integer?Positive,xspan:_Span:;;]:=With[{(*file=metadata[["files",1]],spid=1,xspan=;;,*)group="/moment"},
Module[{attr,time,parent,n,nV,nvv},
attr=Import[file,{"HDF5","Attributes",group}];
time=attr["time"];
parent=FileNameJoin[{group,ToString[spid-1]},OperatingSystem->"Unix"];
attr=KeyDrop[Import[file,{"HDF5","Attributes",parent}],Keys[attr]];
n=Part[Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"n"},OperatingSystem->"Unix"]}],xspan];
nV=Part[Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"nV"},OperatingSystem->"Unix"]}],xspan]//Transpose;
nvv=Part[Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"nvv"},OperatingSystem->"Unix"]}],xspan]//Transpose;
Merge[{attr,{"time"->time,"n"->n},
Thread[Rule[StringJoin/@Thread[{"nV",Characters["123"]}],nV]],
Thread[Rule[StringJoin/@Thread[{"nv",Characters["123"],"v",Characters["123"]}],nvv]]
},Last]
]
]


Clear[importRelativisticMoment]
importRelativisticMoment[file_String?FileExistsQ]:=(Print["importStressEnergyTensor requires the second argument."];Abort[])
importRelativisticMoment[file:_String?FileExistsQ,spid:_Integer?Positive,xspan:_Span:;;]:=With[{(*file=metadata[["files",1]],spid=1,xspan=;;,*)group="/moment"},
Module[{attr,c,time,parent,n,nV,Mij,M00,Mi0,nvv},
attr=Import[file,{"HDF5","Attributes",group}];
c=attr["c"];
time=attr["time"];
parent=FileNameJoin[{group,ToString[spid-1]},OperatingSystem->"Unix"];
attr=KeyDrop[Import[file,{"HDF5","Attributes",parent}],Keys[attr]];
n=Part[Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"n"},OperatingSystem->"Unix"]}],xspan];
nV=Part[Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"nV"},OperatingSystem->"Unix"]}],xspan]//Transpose;
Mij=Part[Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"Mij"},OperatingSystem->"Unix"]}],xspan]//Transpose;
M00=Mij[[1]];
Mi0=Mij[[2;;4]];
nvv=Mij[[5;;7]];
Merge[{attr,{"time"->time,"n"->n,"n\[Gamma]cc"->M00},
Thread[Rule[StringJoin/@Thread[{"nV",Characters["123"]}],nV]],
Thread[Rule[StringJoin/@Thread[{"n\[Gamma]v",Characters["123"],"c"}],Mi0]],
Thread[Rule[StringJoin/@Thread[{"n\[Gamma]v",Characters["123"],"v",Characters["123"]}],nvv]]
},Last]
]
]


Clear[importParticle]
importParticle[file_String?FileExistsQ]:=(Print["importParticle requires the second argument."];Abort[])
importParticle[file:_String?FileExistsQ,spid:_Integer?Positive]:=With[{(*file=metadata[["files",1]],spid=1,*)group="/particle"},
Module[{attr,time,parent,n,nV,nvv,id,x,v,psd},
attr=Import[file,{"HDF5","Attributes",group}];
time=attr["time"];
parent=FileNameJoin[{group,ToString[spid-1]},OperatingSystem->"Unix"];
attr=KeyDrop[Import[file,{"HDF5","Attributes",parent}],Keys[attr]];
n=Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"n"},OperatingSystem->"Unix"]}];
nV=Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"nV"},OperatingSystem->"Unix"]}]//Transpose;
nvv=Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"nvv"},OperatingSystem->"Unix"]}]//Transpose;
id=Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"id"},OperatingSystem->"Unix"]}];
x=Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"pos"},OperatingSystem->"Unix"]}];
v=Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"vel"},OperatingSystem->"Unix"]}]//Transpose;
psd=Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"psd"},OperatingSystem->"Unix"]}]//Transpose;
Merge[{attr,{"time"->time,"n"->n,"id"->id,"q1"->x},
Thread[Rule[StringJoin/@Thread[{"nV",Characters["123"]}],nV]],
Thread[Rule[StringJoin/@Thread[{"nv",Characters["123"],"v",Characters["123"]}],nvv]],
Thread[Rule[StringJoin/@Thread[{"v",Characters["123"]}],v]],
Thread[Rule[{"w","f","g"},psd]]
},Last]
]
]


Clear[importRelativisticParticle]
importRelativisticParticle[file_String?FileExistsQ]:=(Print["importRelativisticParticle requires the second argument."];Abort[])
importRelativisticParticle[file:_String?FileExistsQ,spid:_Integer?Positive]:=With[{(*file=metadata[["files",1]],spid=1,*)group="/particle"},
Module[{attr,c,time,parent,n,nV,Mij,M00,Mi0,nvv,id,x,\[Gamma]c\[Gamma]v,psd},
attr=Import[file,{"HDF5","Attributes",group}];
c=attr["c"];
time=attr["time"];
parent=FileNameJoin[{group,ToString[spid-1]},OperatingSystem->"Unix"];
attr=KeyDrop[Import[file,{"HDF5","Attributes",parent}],Keys[attr]];
n=Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"n"},OperatingSystem->"Unix"]}];
nV=Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"nV"},OperatingSystem->"Unix"]}]//Transpose;
Mij=Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"Mij"},OperatingSystem->"Unix"]}]//Transpose;
M00=Mij[[1]];
Mi0=Mij[[2;;4]];
nvv=Mij[[5;;7]];
id=Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"id"},OperatingSystem->"Unix"]}];
x=Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"pos"},OperatingSystem->"Unix"]}];
\[Gamma]c\[Gamma]v=Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"gcgvel"},OperatingSystem->"Unix"]}]//Transpose;
psd=Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"psd"},OperatingSystem->"Unix"]}]//Transpose;
Merge[{attr,{"time"->time,"n"->n,"id"->id,"q1"->x,"n\[Gamma]cc"->M00,"\[Gamma]c"->First[\[Gamma]c\[Gamma]v]},
Thread[Rule[StringJoin/@Thread[{"nV",Characters["123"]}],nV]],
Thread[Rule[StringJoin/@Thread[{"n\[Gamma]v",Characters["123"],"c"}],Mi0]],
Thread[Rule[StringJoin/@Thread[{"n\[Gamma]v",Characters["123"],"v",Characters["123"]}],nvv]],
Thread[Rule[StringJoin/@Thread[{"\[Gamma]v",Characters["123"]}],Rest[\[Gamma]c\[Gamma]v]]],
Thread[Rule[{"w","f","g"},psd]]
},Last]
]
]


Clear[calVelocityHistogram]
calVelocityHistogram[ptl_?AssociationQ,{v1lim_Interval,Nv1_Integer?Positive},{v2lim_Interval,Nv2_Integer?Positive},q1lim:_Interval:Interval[{-\[Infinity],\[Infinity]}]]:=
With[{v1spec=Append[MinMax[v1lim],Divide[-Subtract@@MinMax[v1lim],Nv1]],v2spec=Append[MinMax[v2lim],Divide[-Subtract@@MinMax[v2lim],Nv2]]},
Module[{extent,interpV1,interpV2,interpV3,q1,pos,v1,v2,v3,binning,bins,g,f,w,df,g1,f1,df1,g2,f2,df2},
extent=Array[N,ptl["Nx"]+2,ptl["half_grid_domain_extent"]+{-1,0}];
interpV1=Function[Interpolation[Thread[{extent,Join[Take[#,-1],#,Take[#,1]]}]]][ptl["nV1"]/ptl["n"]];
interpV2=Function[Interpolation[Thread[{extent,Join[Take[#,-1],#,Take[#,1]]}]]][ptl["nV2"]/ptl["n"]];
interpV3=Function[Interpolation[Thread[{extent,Join[Take[#,-1],#,Take[#,1]]}]]][ptl["nV3"]/ptl["n"]];
q1=ptl["q1"];
pos=Position[q1,$_Real/;IntervalMemberQ[q1lim,$]];
pos=If[Length[pos]==0,(Print["no particles selected"];Abort[]),pos//Extract];
q1=pos[q1];
w=pos[ptl["w"]];
f=pos[ptl["f"]];
g=pos[ptl["g"]];
v1=pos[ptl["v1"]]-interpV1[q1];
v2=pos[ptl["v2"]]-interpV2[q1];
v3=pos[ptl["v3"]]-interpV3[q1];
v2=Sqrt[v2^2+v3^2];
v3=.;
binning=Compile[{{v1,_Real},{v2,_Real},{w,_Real},{f,_Real}},
Module[{i1=0,i2=0},
i1=Floor[(v1-v1spec[[1]])/Last[v1spec]]+1;
i2=Floor[(v2-v2spec[[1]])/Last[v2spec]]+1;
{{i2,i1},{w,f}}
],
Parallelization->True,
RuntimeAttributes->{Listable}
];
bins=binning[v1,v2,w,f/g];
v2=MovingAverage[Array[N,Nv2+1,MinMax[v2lim]],2];
v1=MovingAverage[Array[N,Nv1+1,MinMax[v1lim]],2];
f=Cases[
Map[Rule[Round[Part[#,1,1]],Prepend@@Through[{Total,Length}[Last/@#]]]&,SplitBy[SortBy[bins,Most],Most]],
HoldPattern[Rule[{$2_Integer?Positive/;$2<=Length[v2],$1_Integer?Positive/;$1<=Length[v1]},_List]]
];
f=ReplacePart[Table[{0,0,0},{v2},{v1}],f];
{g1,df1,f1}=Total[f]\[Transpose];
{g2,df2,f2}=Map[Total,f]\[Transpose];
{g,df,f}=Flatten[f,{{3},{1},{2}}];
g1/=Length[q1]Last[v1spec];
f1/=Length[q1]Last[v1spec];
df1/=Length[q1]Last[v1spec];
g2/=Length[q1]Last[v2spec](2Pi v2);
f2/=Length[q1]Last[v2spec](2Pi v2);
df2/=Length[q1]Last[v2spec](2Pi v2);
g/=Length[q1]Last[v1spec]Last[v2spec](2Pi v2);
f/=Length[q1]Last[v1spec]Last[v2spec](2Pi v2);
df/=Length[q1]Last[v1spec]Last[v2spec](2Pi v2);
(*return*)
Association["v1"->v1,"v2"->v2,"f"->f,"g"->g,"df"->df,"f1"->f1,"g1"->g1,"df1"->df1,"f2"->f2,"g2"->g2,"df2"->df2,"time"->ptl["time"],"q1lim"->MinMax[q1]]
]
]


Clear[calMomentumHistogram]
calMomentumHistogram[ptl_?AssociationQ,{\[Gamma]v1lim_Interval,N\[Gamma]v1_Integer?Positive},{\[Gamma]v2lim_Interval,N\[Gamma]v2_Integer?Positive},q1lim:_Interval:Interval[{-\[Infinity],\[Infinity]}]]:=
With[{\[Gamma]v1spec=Append[MinMax[\[Gamma]v1lim],Divide[-Subtract@@MinMax[\[Gamma]v1lim],N\[Gamma]v1]],\[Gamma]v2spec=Append[MinMax[\[Gamma]v2lim],Divide[-Subtract@@MinMax[\[Gamma]v2lim],N\[Gamma]v2]]},
Module[{extent,interpV1,interpV2,interpV3,q1,pos,w,boost,\[Gamma]v1,\[Gamma]v2,\[Gamma]v3,binning,bins,g,f,df,g1,f1,df1,g2,f2,df2},
extent=Array[N,ptl["Nx"]+2,ptl["half_grid_domain_extent"]+{-1,0}];
interpV1=Function[Interpolation[Thread[{extent,Join[Take[#,-1],#,Take[#,1]]}]]][ptl["nV1"]/ptl["n"]];
interpV2=Function[Interpolation[Thread[{extent,Join[Take[#,-1],#,Take[#,1]]}]]][ptl["nV2"]/ptl["n"]];
interpV3=Function[Interpolation[Thread[{extent,Join[Take[#,-1],#,Take[#,1]]}]]][ptl["nV3"]/ptl["n"]];
q1=ptl["q1"];
pos=Position[q1,$_Real/;IntervalMemberQ[q1lim,$]];
pos=If[Length[pos]==0,(Print["no particles selected"];Abort[]),pos//Extract];
q1=pos[q1];
w=pos[ptl["w"]];
f=pos[ptl["f"]];
g=pos[ptl["g"]];
\[Gamma]v1=ptl["\[Gamma]v1"];
\[Gamma]v2=ptl["\[Gamma]v2"];
\[Gamma]v3=ptl["\[Gamma]v3"];
boost=With[{c=ptl["c"]},
Compile[{{\[Gamma]v1,_Real},{\[Gamma]v2,_Real},{\[Gamma]v3,_Real},{V1,_Real},{V2,_Real},{V3,_Real}},
Module[{V,\[DoubleStruckN],\[Gamma]d,\[Gamma]u,\[Gamma]\[DoubleStruckV]={\[Gamma]v1,\[Gamma]v2,\[Gamma]v3}},
V=Sqrt[V1^2+V2^2+V3^2];
If[V < 1.0*^-10 c, \[Gamma]\[DoubleStruckV],
\[Gamma]d=c/Sqrt[(c-V)(c+V)];
\[DoubleStruckN]={V1,V2,V3}/V;
\[Gamma]u=Sqrt[1+\[Gamma]\[DoubleStruckV] . \[Gamma]\[DoubleStruckV]/c^2];
\[Gamma]\[DoubleStruckV]+(\[Gamma]d-1)Dot[\[Gamma]\[DoubleStruckV],\[DoubleStruckN]]\[DoubleStruckN]-\[Gamma]u \[Gamma]d V \[DoubleStruckN]
]
],
Parallelization->True,
RuntimeAttributes->{Listable}
]
];
\[Gamma]v3=boost[\[Gamma]v1,\[Gamma]v2,\[Gamma]v3,interpV1[q1],interpV2[q1],interpV3[q1]];
\[Gamma]v1=\[Gamma]v3[[All,1]];
\[Gamma]v2=\[Gamma]v3[[All,2]];
\[Gamma]v3=\[Gamma]v3[[All,3]];
\[Gamma]v2=Sqrt[\[Gamma]v2^2+\[Gamma]v3^2];
\[Gamma]v3=.;
binning=Compile[{{\[Gamma]v1,_Real},{\[Gamma]v2,_Real},{w,_Real},{f,_Real}},
Module[{i1=0,i2=0},
i1=Floor[(\[Gamma]v1-\[Gamma]v1spec[[1]])/Last[\[Gamma]v1spec]]+1;
i2=Floor[(\[Gamma]v2-\[Gamma]v2spec[[1]])/Last[\[Gamma]v2spec]]+1;
{{i2,i1},{w,f}}
],
Parallelization->True,
RuntimeAttributes->{Listable}
];
bins=binning[\[Gamma]v1,\[Gamma]v2,w,f/g];
\[Gamma]v2=MovingAverage[Array[N,N\[Gamma]v2+1,MinMax[\[Gamma]v2lim]],2];
\[Gamma]v1=MovingAverage[Array[N,N\[Gamma]v1+1,MinMax[\[Gamma]v1lim]],2];
f=Cases[
Map[Rule[Round[Part[#,1,1]],Prepend@@Through[{Total,Length}[Last/@#]]]&,SplitBy[SortBy[bins,Most],Most]],
HoldPattern[Rule[{$2_Integer?Positive/;$2<=Length[\[Gamma]v2],$1_Integer?Positive/;$1<=Length[\[Gamma]v1]},_List]]
];
f=ReplacePart[Table[{0,0,0},{\[Gamma]v2},{\[Gamma]v1}],f];
{g1,df1,f1}=Total[f]\[Transpose];
{g2,df2,f2}=Map[Total,f]\[Transpose];
{g,df,f}=Flatten[f,{{3},{1},{2}}];
g1/=Length[q1]Last[\[Gamma]v1spec];
f1/=Length[q1]Last[\[Gamma]v1spec];
df1/=Length[q1]Last[\[Gamma]v1spec];
g2/=Length[q1]Last[\[Gamma]v2spec](2Pi \[Gamma]v2);
f2/=Length[q1]Last[\[Gamma]v2spec](2Pi \[Gamma]v2);
df2/=Length[q1]Last[\[Gamma]v2spec](2Pi \[Gamma]v2);
g/=Length[q1]Last[\[Gamma]v1spec]Last[\[Gamma]v2spec](2Pi \[Gamma]v2);
f/=Length[q1]Last[\[Gamma]v1spec]Last[\[Gamma]v2spec](2Pi \[Gamma]v2);
df/=Length[q1]Last[\[Gamma]v1spec]Last[\[Gamma]v2spec](2Pi \[Gamma]v2);
(*return*)
Association["\[Gamma]v1"->\[Gamma]v1,"\[Gamma]v2"->\[Gamma]v2,"f"->f,"g"->g,"df"->df,"f1"->f1,"g1"->g1,"df1"->df1,"f2"->f2,"g2"->g2,"df2"->df2,"time"->ptl["time"],"q1lim"->MinMax[q1]]
]
]


Clear[importVHistogram]
importVHistogram[file_String?FileExistsQ]:=(Print["importVHistogram requires the second argument."];Abort[])
importVHistogram[file:_String?FileExistsQ,spid:_Integer?Positive]:=With[{(*file=metadata[["files",1]],spid=1,*)group="/vhist2d"},
Module[{attr,time,parent,idx,ghist,whist,fhist,v1,v2,dV,psd,gpsd,wpsd,fpsd},
attr=Import[file,{"HDF5","Attributes",group}];
time=attr["time"];
parent=FileNameJoin[{group,ToString[spid-1]},OperatingSystem->"Unix"];
attr=KeyDrop[Import[file,{"HDF5","Attributes",parent}],Keys[attr]];
idx=Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"idx"},OperatingSystem->"Unix"]}];
psd=Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"psd"},OperatingSystem->"Unix"]}]//Transpose;
{ghist,whist,fhist}=psd;
ghist=Normal@SparseArray[idx->ghist,attr["vdims"]];
whist=Normal@SparseArray[idx->whist,attr["vdims"]];
fhist=Normal@SparseArray[idx->fhist,attr["vdims"]];
v1=Array[N,attr[["vdims",1]]+1,attr["v1lim"]]~MovingAverage~2;
v2=Array[N,attr[["vdims",2]]+1,attr["v2lim"]]~MovingAverage~2;
dV=2\[Pi] v2(v1[[2]]-v1[[1]])(v2[[2]]-v2[[1]]);
gpsd=Divide[#,dV]&/@ghist;
wpsd=Divide[#,dV]&/@whist;
fpsd=Divide[#,dV]&/@fhist;
Merge[{attr,{"time"->time,"v1"->v1,"v2"->v2,"vhist"->ghist,"whist"->whist,"fhist"->fhist,"vpsd"->gpsd,"wpsd"->wpsd,"fpsd"->fpsd}},Last]
]
]


Clear[importRelativisticVHistogram]
importRelativisticVHistogram[file_String?FileExistsQ]:=(Print["importVHistogram requires the second argument."];Abort[])
importRelativisticVHistogram[file:_String?FileExistsQ,spid:_Integer?Positive]:=With[{(*file=metadata[["files",1]],spid=1,*)group="/gvhist2d"},
Module[{attr,time,parent,idx,ghist,whist,fhist,gv1,gv2,dV,psd,gpsd,wpsd,fpsd},
attr=Import[file,{"HDF5","Attributes",group}];
time=attr["time"];
parent=FileNameJoin[{group,ToString[spid-1]},OperatingSystem->"Unix"];
attr=KeyDrop[Import[file,{"HDF5","Attributes",parent}],Keys[attr]];
idx=Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"idx"},OperatingSystem->"Unix"]}];
psd=Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"psd"},OperatingSystem->"Unix"]}]//Transpose;
{ghist,whist,fhist}=psd;
ghist=Normal@SparseArray[idx->ghist,attr["gvdims"]];
whist=Normal@SparseArray[idx->whist,attr["gvdims"]];
fhist=Normal@SparseArray[idx->fhist,attr["gvdims"]];
gv1=Array[N,attr[["gvdims",1]]+1,attr["gv1lim"]]~MovingAverage~2;
gv2=Array[N,attr[["gvdims",2]]+1,attr["gv2lim"]]~MovingAverage~2;
dV=2\[Pi] gv2(gv1[[2]]-gv1[[1]])(gv2[[2]]-gv2[[1]]);
gpsd=Divide[#,dV]&/@ghist;
wpsd=Divide[#,dV]&/@whist;
fpsd=Divide[#,dV]&/@fhist;
Merge[{attr,{"time"->time,"\[Gamma]v1"->gv1,"\[Gamma]v2"->gv2,"vhist"->ghist,"whist"->whist,"fhist"->fhist,"vpsd"->gpsd,"wpsd"->wpsd,"fpsd"->fpsd}},Last]
]
]
