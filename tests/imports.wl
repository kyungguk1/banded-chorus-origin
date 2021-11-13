(* ::Package:: *)

Clear[importField]
importField[file:_String?FileExistsQ,xspan:_Span:;;]:=With[{(*file=metadata[["files",1]],*)group="/field",names={"B","E"}},
Module[{attr,time,dsets},
attr=Import[file,{"HDF5","Attributes",FileNameJoin[{group,First[names]}]}];
time=attr["time"];
dsets=Import[file,{"HDF5","Datasets",FileNameJoin[{group,#}]&/@names}];
Association[Append[Thread[Rule[StringJoin/@Tuples[{{"d"},names,Characters["123"]}],Part[#,xspan]&/@Apply[Join,Transpose/@dsets]]],"time"->time]]
]
]


Clear[importMoment]
importMoment[file_String?FileExistsQ]:=(Print["importMoment requires the second argument."];Abort[])
importMoment[file:_String?FileExistsQ,spid:_Integer?Positive,xspan:_Span:;;]:=With[{(*file=metadata[["files",1]],spid=1,xspan=;;,*)group="/moment"},
Module[{attr,time,parent,n,nV,nvv},
attr=Import[file,{"HDF5","Attributes",group}];
time=attr["time"];
parent=FileNameJoin[{group,ToString[spid-1]}];
attr=KeyDrop[Import[file,{"HDF5","Attributes",parent}],Keys[attr]];
n=Part[Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"n"}]}],xspan];
nV=Part[Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"nV"}]}],xspan]//Transpose;
nvv=Part[Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"nvv"}]}],xspan]//Transpose;
Merge[{attr,{"time"->time,"n"->n},
Thread[Rule[StringJoin/@Thread[{"nV",Characters["123"]}],nV]],
Thread[Rule[StringJoin/@Thread[{"nv",Characters["123"],"v",Characters["123"]}],nvv]]
},Last]
]
]


Clear[importStressEnergyTensor]
importStressEnergyTensor[file_String?FileExistsQ]:=(Print["importStressEnergyTensor requires the second argument."];Abort[])
importStressEnergyTensor[file:_String?FileExistsQ,spid:_Integer?Positive,xspan:_Span:;;]:=With[{(*file=metadata[["files",1]],spid=1,xspan=;;,*)group="/moment"},
Module[{attr,c,time,parent,n,nV,Mij,M00,Mi0,nvv},
attr=Import[file,{"HDF5","Attributes",group}];
c=attr["c"];
time=attr["time"];
parent=FileNameJoin[{group,ToString[spid-1]}];
attr=KeyDrop[Import[file,{"HDF5","Attributes",parent}],Keys[attr]];
n=Part[Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"n"}]}],xspan];
nV=Part[Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"nV"}]}],xspan]//Transpose;
Mij=Part[Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"Mij"}]}],xspan]//Transpose;
M00=Mij[[1]];
Mi0=Mij[[2;;4]];
nvv=Mij[[5;;7]];
Merge[{attr,{"time"->time,"n"->n,"n\[Gamma]c2"->M00},
Thread[Rule[StringJoin/@Thread[{"nV",Characters["123"]}],nV]],
Thread[Rule[StringJoin/@Thread[{"nU",Characters["123"]}],Mi0/c]],
Thread[Rule[StringJoin/@Thread[{"nu",Characters["123"],"v",Characters["123"]}],nvv]]
},Last]
]
]


Clear[importParticle]
importParticle[file_String?FileExistsQ]:=(Print["importParticle requires the second argument."];Abort[])
importParticle[file:_String?FileExistsQ,spid:_Integer?Positive]:=With[{(*file=metadata[["files",1]],spid=1,*)group="/particle"},
Module[{attr,time,parent,n,nV,nvv,id,x,v,psd},
attr=Import[file,{"HDF5","Attributes",group}];
time=attr["time"];
parent=FileNameJoin[{group,ToString[spid-1]}];
attr=KeyDrop[Import[file,{"HDF5","Attributes",parent}],Keys[attr]];
n=Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"n"}]}];
nV=Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"nV"}]}]//Transpose;
nvv=Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"nvv"}]}]//Transpose;
id=Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"id"}]}];
x=Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"pos"}]}];
v=Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"vel"}]}]//Transpose;
psd=Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"psd"}]}]//Transpose;
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
Module[{attr,c,time,parent,n,nV,Mij,M00,Mi0,nvv,id,x,v,psd},
attr=Import[file,{"HDF5","Attributes",group}];
c=attr["c"];
time=attr["time"];
parent=FileNameJoin[{group,ToString[spid-1]}];
attr=KeyDrop[Import[file,{"HDF5","Attributes",parent}],Keys[attr]];
n=Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"n"}]}];
nV=Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"nV"}]}]//Transpose;
Mij=Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"Mij"}]}]//Transpose;
M00=Mij[[1]];
Mi0=Mij[[2;;4]];
nvv=Mij[[5;;7]];
id=Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"id"}]}];
x=Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"pos"}]}];
v=Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"vel"}]}]//Transpose;
psd=Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"psd"}]}]//Transpose;
Merge[{attr,{"time"->time,"n"->n,"id"->id,"q1"->x,"n\[Gamma]c2"->M00},
Thread[Rule[StringJoin/@Thread[{"nV",Characters["123"]}],nV]],
Thread[Rule[StringJoin/@Thread[{"nU",Characters["123"]}],Mi0/c]],
Thread[Rule[StringJoin/@Thread[{"nu",Characters["123"],"v",Characters["123"]}],nvv]],
Thread[Rule[StringJoin/@Thread[{"v",Characters["123"]}],v]],
Thread[Rule[{"w","f","g"},psd]]
},Last]
]
]


Clear[calHist]
calHist[ptl_Association,vlim_Interval,Nv_Integer?Positive]:=(*With[{ptl=particle,vlim=.4Interval[{-1,1}],Nv=150},*)
With[{v1spec=Append[MinMax[vlim],Divide[-Subtract@@MinMax[vlim],Nv]],v2spec={0.,Max[vlim],Divide[-Subtract@@MinMax[vlim],Nv]}},
Module[{extent,interpV1,interpV2,interpV3,q1,v1,v2,v3,w,binning,bins,f,df,f1,df1,f2,df2},
extent=Array[N,ptl["Nx"]+1,ptl["full_grid_domain_extent"]];
interpV1=Function[Interpolation[Thread[{extent,Append[#,First[#]]}]]][ptl["nV1"]/ptl["n"]];
interpV2=Function[Interpolation[Thread[{extent,Append[#,First[#]]}]]][ptl["nV2"]/ptl["n"]];
interpV3=Function[Interpolation[Thread[{extent,Append[#,First[#]]}]]][ptl["nV3"]/ptl["n"]];
w=ptl["w"];
q1=ptl["q1"];
v1=ptl["v1"]-interpV1[q1];
v2=ptl["v2"]-interpV2[q1];
v3=ptl["v3"]-interpV3[q1];
v2=Sqrt[v2^2+v3^2];
v3=.;
binning=Compile[{{v1,_Real},{v2,_Real},{w,_Real}},
Module[{i1=0,i2=0},
i1=Floor[(v1-v1spec[[1]])/Last[v1spec]]+1;
i2=Floor[(v2-v2spec[[1]])/Last[v2spec]]+1;
{i1,w,i2}
],
Parallelization->True,
RuntimeAttributes->{Listable}
];
bins=binning[v1,v2,w];
v2=N@MovingAverage[Range@@v2spec,2];
v1=N@MovingAverage[Range@@v1spec,2];
f=Cases[
Rule[Round[Most@Part[#,1]],Through[{Length,Total}[Last/@#]]]&/@SplitBy[SortBy[(*{i2,i1,w}*)RotateRight/@bins,Most],Most],
HoldPattern[Rule[{$2_Integer?Positive/;$2<=Length[v2],$1_Integer?Positive/;$1<=Length[v1]},_List]]
];
f=ReplacePart[Table[{0,0},{v2},{v1}],f];
{f1,df1}=Total[f]\[Transpose];
{f2,df2}=Map[Total,f]\[Transpose];
{f,df}=Flatten[f,{{3},{1},{2}}];
f1/=Length[w]Last[v1spec];
df1/=Length[w]Last[v1spec];
f2/=Length[w]Last[v2spec](2Pi v2);
df2/=Length[w]Last[v2spec](2Pi v2);
f/=Length[w]Last[v1spec]Last[v2spec](2Pi v2);
df/=Length[w]Last[v1spec]Last[v2spec](2Pi v2);
(*return*)
Association["v1"->v1,"v2"->v2,"f"->f,"df"->df,"f1"->f1,"df1"->df1,"f2"->f2,"df2"->df2,"time"->ptl["time"]]
]
]


Clear[importVHistogram]
importVHistogram[file_String?FileExistsQ]:=(Print["importVHistogram requires the second argument."];Abort[])
importVHistogram[file:_String?FileExistsQ,spid:_Integer?Positive]:=With[{(*file=metadata[["files",1]],spid=1,*)group="/vhist2d"},
Module[{attr,time,parent,idx,vhist,whist,v1,v2,dV,vpsd,wpsd},
attr=Import[file,{"HDF5","Attributes",group}];
time=attr["time"];
parent=FileNameJoin[{group,ToString[spid-1]}];
attr=KeyDrop[Import[file,{"HDF5","Attributes",parent}],Keys[attr]];
idx=Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"idx"}]}];
{vhist,whist}=Import[file,{"HDF5","Datasets",FileNameJoin[{parent,"psd"}]}]//Transpose;
vhist=Normal@SparseArray[idx->vhist,attr["vdims"]];
whist=Normal@SparseArray[idx->whist,attr["vdims"]];
v1=Array[N,attr[["vdims",1]]+1,attr["v1lim"]]~MovingAverage~2;
v2=Array[N,attr[["vdims",2]]+1,attr["v2lim"]]~MovingAverage~2;
dV=2\[Pi] v2(Subtract@@v1[[{2,1}]])(Subtract@@v2[[{2,1}]]);
vpsd=Divide[#,dV]&/@vhist;
wpsd=Divide[#,dV]&/@whist;
Merge[{attr,{"time"->time,"v1"->v1,"v2"->v2,"vhist"->vhist,"whist"->whist,"vpsd"->vpsd,"wpsd"->wpsd}},Last]
]
]
