(* ::Package:: *)

Options[BoutCollect]={xind->Null,yind->Null,zind->Null,tind->Null,path->".",yguards->False,info->True,prefix->"BOUT.dmp"};

BoutCollect[varname_,OptionsPattern[]]:=
Module[
	{Xind,Yind,Zind,Tind,Path,Yguards,Info,Prefix,
		varnameissymbol,vars,position,dimensions,nxpe,nype,mxsub,mysub,mxg,myg,mz,tarray,files,nfiles,data,tempdata,ts,te,xs,xe,ys,ye,zs,ze,localx,localy,import,lxs,lxe,lys,lye},

    Throw["This is currently broken for BOUT++ > v4.0.0. See issue #394"]

	Xind=OptionValue[xind];
	Yind=OptionValue[yind];
	Zind=OptionValue[zind];
	Tind=OptionValue[tind];
	Path=OptionValue[path];
	Yguards=OptionValue[yguards];
	Info=OptionValue[info];
	Prefix=OptionValue[prefix];
	varnameissymbol=False;
	If[ValueQ[varname] || StringQ[varname],1,varname=SymbolName[varname];varnameissymbol=True];
	If[StringQ[varname],1,Print["Variable name error"];Return[]];
	files=FileNames[StringJoin[Prefix,"*.nc"],Path];
	nfiles=Length[files];
	vars=Import[files[[1]]];
	position=Position[vars,varname][[1,1]];
	nxpe=Import[files[[1]],{"Data",Position[vars,"NXPE"][[1,1]]}];
	nype=Import[files[[1]],{"Data",Position[vars,"NYPE"][[1,1]]}];
	mxsub=Import[files[[1]],{"Data",Position[vars,"MXSUB"][[1,1]]}];
	mysub=Import[files[[1]],{"Data",Position[vars,"MYSUB"][[1,1]]}];
	mxg=Import[files[[1]],{"Data",Position[vars,"MXG"][[1,1]]}];
	myg=Import[files[[1]],{"Data",Position[vars,"MYG"][[1,1]]}];
	mz=Import[files[[1]],{"Data",Position[vars,"MZ"][[1,1]]}];
	tarray=Import[files[[1]],{"Data",Position[vars,"t_array"][[1,1]]}];
	If[ListQ[Tind],ts=Tind[[1]];te=Tind[[2]],If[NumberQ[Tind],ts=Tind;te=Tind,ts=1;te=Length[tarray]]];
	If[ListQ[Xind],xs=Xind[[1]];xe=Xind[[2]],If[NumberQ[Xind],xs=Xind;xe=Xind,xs=1;xe=nxpe mxsub+2mxg]];
	If[ListQ[Yind],ys=Yind[[1]];ye=Yind[[2]],If[NumberQ[Yind],ys=Yind;ye=Yind,If[Yguards,ys=1;ye=nype mysub+2myg,ys=myg+1;ye=nype mysub+myg]]];
	If[ListQ[Zind],zs=Zind[[1]];ze=Zind[[2]],If[NumberQ[Zind],zs=Zind;ze=Zind,zs=1;ze=-2]];
	localx[jx_,xp_]:=Module[{temp},temp=jx-xp mxsub;Which[temp<=mxg&&xp!=0,-Infinity,temp>mxsub+mxg&&xp!=nxpe-1,Infinity,True,temp]];
	localy[jy_,yp_]:=Module[{temp},temp=jy-yp mysub;Which[temp<=myg&&yp!=0,-Infinity,temp>mysub+myg&&yp!=nype-1,Infinity,True,temp]];
	import[jxs_,jxe_,jys_,jye_,xp_,yp_,jts_:Null,jte_:Null]:=Module[{indata},
		indata=Import[StringJoin[Path,"/",Prefix,".",ToString[xp nype+yp],".nc"],{"Data",position}];
		If[NumberQ[jts]&&NumberQ[jte],indata[[jts;;jte,jxs;;jxe,jys;;jye]],indata[[jxs;;jxe,jys;;jye]]]
	];
	dimensions=Import[files[[1]],{"Dimensions",position}];
	If[nfiles==1,
		data=Import[files[[1]],{"Data",position}];
		Switch[Length[dimensions],
			0,1,
			1,data=data[[ts;;te]],
			2,data=data[[xs;;xe,ys;;ye]],
			3,If[Length[tarray]==mz&&mxsub+2mxg==mz&&mysub+2myg==mz,Print["Cannot tell if array is [T,X,Y] or [X,Y,Z]"];Return[]]
				Switch[dimensions,
					{Length[tarray],mxsub+2mxg,mysub+2myg},data=data[[ts;;te,xs;;xe,ys;;ye]],
					{mxsub+2mxg,mysub+2myg,mz},data=data[[xs;;xe,ys;;ye,zs;;ze]]
				],
			4,data=data[[ts;;te,xs;;xe,ys;;ye,zs;;ze]]
		]
	,
		Switch[Length[dimensions],
			0,data=Import[files[[1]],{"Data",position}],
			1,data=Import[files[[1]],{"Data",position}][[ts;;te]],
			2,data={};
				Do[
					tempdata={{}};
					Do[
						lxs=localx[xs,xp];
						lxe=localx[xe,xp];
						lys=localy[ys,yp];
						lye=localy[ye,yp];
						If[lxs!=Infinity&&lxe!=-Infinity&&lys!=Infinity&&lye!=-Infinity,
							If[lxs==-Infinity,lxs=mxg+1];
							If[lxe==Infinity,lxe=mxsub+mxg];
							If[lys==-Infinity,lys=myg+1];
							If[lye==Infinity,lye=mysub+myg];
							tempdata=Join[tempdata,import[lxs,lxe,lys,lye,xp,yp],2]
						]
					,{yp,0,nype-1}];
					data=Join[data,tempdata];
				,{xp,0,nxpe-1}],
			3,If[Length[tarray]==mz&&mxsub+2mxg==mz&&mysub+2myg==mz,Print["Cannot tell if array is [T,X,Y] or [X,Y,Z]"];Return[]];
				Switch[dimensions,
					{Length[tarray],mxsub+2mxg,mysub+2myg},
						data={{}};
						Do[
							tempdata={{{}}};
							Do[
								lxs=localx[xs,xp];
								lxe=localx[xe,xp];
								lys=localy[ys,yp];
								lye=localy[ye,yp];
								If[lxs!=Infinity&&lxe!=-Infinity&&lys!=Infinity&&lye!=-Infinity,
									If[lxs==-Infinity,lxs=mxg+1];
									If[lxe==Infinity,lxe=mxsub+mxg];
									If[lys==-Infinity,lys=myg+1];
									If[lye==Infinity,lye=mysub+myg];
									tempdata=Join[tempdata,Evaluate[import[lxs,lxe,lys,lye,xp,yp,ts,te]],3]
								]
							,{yp,0,nype-1}];
							data=Join[data,tempdata,2];
						,{xp,0,nxpe-1}],
					{mxsub+2mxg,mysub+2myg,mz},
						data={};
						Do[
							tempdata={{}};
							Do[
								lxs=localx[xs,xp];
								lxe=localx[xe,xp];
								lys=localy[ys,yp];
								lye=localy[ye,yp];
								If[lxs!=Infinity&&lxe!=-Infinity&&lys!=Infinity&&lye!=-Infinity,
									If[lxs==-Infinity,lxs=mxg+1];
									If[lxe==Infinity,lxe=mxsub+mxg];
									If[lys==-Infinity,lys=myg+1];
									If[lye==Infinity,lye=mysub+myg];
									tempdata=Join[tempdata,Evaluate[import[lxs,lxe,lys,lye,xp,yp]],2]
								]
							,{yp,0,nype-1}];
							data=Join[data,tempdata];
						,{xp,0,nxpe-1}];
						data=data[[;;,;;,zs;;ze]];
				],
			4,data={{}};
				Do[
					tempdata={{{}}};
					Do[
						lxs=localx[xs,xp];
						lxe=localx[xe,xp];
						lys=localy[ys,yp];
						lye=localy[ye,yp];
						If[lxs!=Infinity&&lxe!=-Infinity&&lys!=Infinity&&lye!=-Infinity,
							If[lxs==-Infinity,lxs=mxg+1];
							If[lxe==Infinity,lxe=mxsub+mxg];
							If[lys==-Infinity,lys=myg+1];
							If[lye==Infinity,lye=mysub+myg];
							tempdata=Join[tempdata,import[lxs,lxe,lys,lye,xp,yp,ts,te],3]
						]
					,{yp,0,nype-1}];
					data=Join[data,tempdata,2];
				,{xp,0,nxpe-1}];
				data=data[[;;,;;,;;,zs;;ze]];
		];
	];
	If[varnameissymbol,Clear[varname]];
	data
]

bc=BoutCollect
