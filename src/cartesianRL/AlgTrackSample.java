//
package cartesianRL;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.function.Predicate;
import fltDispersion.FLTUtils;
import fltDispersion.FltParticle;
import miniufo.application.basic.DynamicMethodsInCTS;
import miniufo.application.contour.ContourCartesianSpatialModel;
import miniufo.application.statisticsModel.LagrangianStatisticsByDavis;
import miniufo.application.statisticsModel.SingleParticleStatResult;
import miniufo.basic.ArrayUtil;
import miniufo.concurrent.ConcurrentUtil;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.CartesianSpatialModel;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.MDate;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.io.CtlDataWriteStream;
import miniufo.io.DataIOFactory;
import miniufo.io.DataWrite;
import miniufo.io.IOUtil;
import miniufo.lagrangian.LagrangianUtil;
import miniufo.lagrangian.Particle;
import miniufo.lagrangian.Record;
import miniufo.statistics.FilterModel;
import miniufo.util.GridDataFetcher;
import miniufo.util.TicToc;

//
public final class AlgTrackSample{
	//
	private static int varGroups=4;
	
	private static String[] uvNames =new String[]{"u","v"};
	private static String[] trcNames=new String[]{"tr1","tr2","tr3","tr4","tr5","tr6","tr7","tr8","tr9","tr10"};
	private static String[] varNames=new String[uvNames.length+trcNames.length];
	private static String[] grdNames=new String[uvNames.length+trcNames.length];
	private static String[] YeqNames=new String[uvNames.length+trcNames.length];
	private static String[] dspNames=new String[uvNames.length+trcNames.length];
	
	private static final String path=Parameters.path;
	
	static{
		for(int i=0,I=uvNames.length;i<I;i++){
			varNames[i]=uvNames[i];
			grdNames[i]=uvNames[i]+"grd";
			YeqNames[i]=uvNames[i]+"Yeq";
			dspNames[i]=uvNames[i]+"dsp";
		}
		for(int i=0,I=trcNames.length;i<I;i++){
			varNames[i+uvNames.length]=trcNames[i];
			grdNames[i+uvNames.length]=trcNames[i]+"grd";
			YeqNames[i+uvNames.length]=trcNames[i]+"Yeq";
			dspNames[i+uvNames.length]=trcNames[i]+"dsp";
		}
	}
	
	
	//
	public static void main(String[] args){
		String test="runK150";
		
		
		/*
		List<FltParticle> fps=FLTUtils.readFLTTrajectory(path+"dispInCC/"+test+"/fltOutput/",new MDate(2000,1,1,0,0),rec->true);
		System.out.println("finish reading "+fps.size()+" FltParticles");
		try(ObjectOutputStream oos=new ObjectOutputStream(new BufferedOutputStream(new FileOutputStream(path+"dispInCC/"+test+"/fltOutput/FltParticle.bin"),819200))){
			for(FltParticle p:fps) oos.writeObject(p);
			oos.writeObject(null);
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
		
		System.exit(0);
		
		List<Particle> ls=new ArrayList<>(55530);
		
		try(ObjectInputStream ois=new ObjectInputStream(new BufferedInputStream(new FileInputStream(path+"dispInCC/"+test+"/fltOutput/FltParticle.bin"),819200))){
			while(true){
				FltParticle p=(FltParticle)(ois.readObject());
				if(p==null) break;
				
				ls.add(FLTUtils.toParticle(p,
					2+varNames.length*varGroups,
					ArrayUtil.concatAll(String.class,new String[]{"uvel","vvel"},uvNames,trcNames,grdNames,dspNames),
					false
				));
			}
			
		}catch(IOException|ClassNotFoundException e){ e.printStackTrace(); System.exit(0);}
		
		samplePTracers(ls,"H:/cartRL_advSchemes/Leith1_k0/Stat.cts");
		tracerValueToEqvY(ls,"H:/cartRL_advSchemes/Leith1_k0/Stat.cts");
		cCrossContourDisplacement(ls);
		
		LagrangianUtil.asRecordStream(ls).forEach(r->{
			r.setData(2,r.getXPos());
			r.setData(3,r.getYPos());
		});
		
		ConcurrentUtil.initDefaultExecutor(2);
		for(int i=1;i<=trcNames.length;i++){
			cYPositionByDisplacement(ls,i);
			calDiffMap(ls,i);
		}
		ConcurrentUtil.shutdown(); System.exit(0);*/
		
		
		
		
		LagrangianStat.regions.forEach(region->{
			/**/
			List<Particle> ps=Parameters.getParticlesDeployedInRegion(
				path+"fltInit_11km_All.bin",path+"dispInCC/"+test+"/fltOutput/",region,2+varNames.length*varGroups,
				ArrayUtil.concatAll(String.class,new String[]{"uvel","vvel"},uvNames,trcNames,grdNames,dspNames)
			);
			
			samplePTracers(ps,"H:/cartRL_advSchemes/Leith1_k0/Stat.cts");
			tracerValueToEqvY(ps,"H:/cartRL_advSchemes/Leith1_k0/Stat.cts");
			cCrossContourDisplacement(ps);
			toAlongTrackFile(ps,path+"dispInCC/"+test+"/algTrack/algTrackData"+region.getName()+".dat");
			//LagrangianUtil.writeTrajecories(ps,path+"dispInCC/"+test+"/TXT/",true,p->true);
			
			LagrangianUtil.writeAsBinaryFile(ps,path+"dispInCC/"+test+"/algTrack/Particles"+region.getName()+".bin");
			
			//List<Particle> ps=LagrangianUtil.readFromBinaryFile(path+"dispInCC/"+test+"/algTrack/Particles"+region.getName()+".bin");
			
			//for(int i=1;i<10;i++){
			//	cYPositionByDisplacement(ps,i);
			//	toSPStatFile(ps,path+"dispInCC/"+test+"/LagStat/CTCoord/"+region.getName(),i);
			//}
		});
	}
	
	
	static void samplePTracers(List<Particle> ps,String cts){
		DiagnosisFactory df=DiagnosisFactory.parseFile(cts);
		DataDescriptor dd=df.getDataDescriptor();
		
		int tLen=ps.get(0).getTCount();
		int ttag=dd.getTNum(ps.get(0).getTime(0))+1;
		
		CartesianSpatialModel csm=new CartesianSpatialModel(dd);
		GridDataFetcher gdf=new GridDataFetcher(dd);
		
		DynamicMethodsInCTS dm=new DynamicMethodsInCTS(csm);
		
		System.out.println("ttag for initial time ("+ps.get(0).getTime(0)+"): "+ttag);
		
		for(int m=0,M=varNames.length;m<M;m++){
			TicToc.tic("samping "+varNames[m]);
			Variable buffer=gdf.prepareXYTBuffer(varNames[m],1,ttag,tLen,0);buffer.replaceUndefData(Record.undef);
			Variable bufGrd=dm.c2DGradientMagnitude(buffer);
			
			final int mm=m+2;
			ps.stream().flatMap(p->p.stream()).forEach(r->{
				float tr=gdf.fetchXYTBuffer(r.getXPos(),r.getYPos(),r.getTime(),buffer);
				float gd=gdf.fetchXYTBuffer(r.getXPos(),r.getYPos(),r.getTime(),bufGrd);
				
				r.setData(mm                ,tr);
				r.setData(mm+varNames.length,gd);
			});
			TicToc.toc(TimeUnit.SECONDS);
		}
		
		gdf.closeFile();
	}
	
	static void tracerValueToEqvY(List<Particle> ps,String cts){
		DiagnosisFactory df=DiagnosisFactory.parseFile(cts);df.setPrinting(false);
		DataDescriptor dd=df.getDataDescriptor();
		
		int tLen=ps.get(0).getTCount();
		int ttag=dd.getTNum(ps.get(0).getTime(0));
		System.out.println("ttag: "+ttag);
		
		ContourCartesianSpatialModel ccsm=new ContourCartesianSpatialModel(dd);
		//Variable[] trs=df.getVariables(new Range("t("+(ttag+1)+","+(ttag+1)+")",dd),varNames);
		
		for(int l=0;l<tLen;l++){
			Variable[] trs=df.getVariables(new Range("t("+(ttag+l+1)+","+(ttag+l+1)+")",dd),varNames);
			
			for(int k=0,K=trs.length;k<K;k++){
				trs[k].replaceUndefData(Record.undef);
				
				float[][] tdata=trs[k].getData()[0][0];
				float[] extrema=ArrayUtil.getExtrema(tdata,trs[k].getUndef());
				
				double[] cVals=new double[ps.size()];
				
				for(int i=0,I=ps.size();i<I;i++){
					Particle p=ps.get(i);
					
					cVals[i]=p.getRecord(l).getDataValue(2+k);
					if(cVals[i]!=Record.undef&&cVals[i]<extrema[0]){
						System.out.println(p.getID()+", z: "+k+", l: "+i+",   cVals[l]: "+cVals[i]+",  min: "+extrema[0]+", "+varNames[k]);
						cVals[i]=extrema[0];
					}
					if(cVals[i]!=Record.undef&&cVals[i]>extrema[1]){
						System.out.println(p.getID()+", z: "+k+", l: "+i+",   cVals[l]: "+cVals[i]+",  max: "+extrema[1]+", "+varNames[k]);
						cVals[i]=extrema[1];
					}
				}
				
				double[] Ys=ccsm.computeEquivalentYs(trs[k],cVals,2);
				
				for(int i=0,I=ps.size();i<I;i++) ps.get(i).getRecord(l).setData(2+k+2*varNames.length,(float)Ys[i]);
			}
		}
		
		/**
		for(int k=0,K=trs.length;k<K;k++){
			trs[k].replaceUndefData(Record.undef);
			
			float[][] tdata=trs[k].getData()[0][0];
			float[] extrema=ArrayUtil.getExtrema(tdata,trs[k].getUndef());
			
			for(Particle p:ps){
				int C=p.getTCount();
				
				double[] cVals=new double[C];
				
				for(int l=0;l<C;l++){
					cVals[l]=p.getRecord(l).getDataValue(k+offset);
					if(cVals[l]!=Record.undef&&cVals[l]<extrema[0]){
						System.out.println(p.getID()+", z: "+k+", l: "+l+",   cVals[l]: "+cVals[l]+",  min: "+extrema[0]);
						cVals[l]=extrema[0];
					}
					if(cVals[l]!=Record.undef&&cVals[l]>extrema[1]){
						System.out.println(p.getID()+", z: "+k+", l: "+l+",   cVals[l]: "+cVals[l]+",  max: "+extrema[1]);
						cVals[l]=extrema[1];
					}
				}
				
				double[] Ys=ccsm.computeEquivalentYs(trs[k],cVals,2);
				
				for(int l=0;l<C;l++) p.getRecord(l).setData(k+offset,(float)Ys[l]);
			}
		}*/
	}
	
	static void cCrossContourDisplacement(List<Particle> ps){
		System.out.println("cross-contour displacement");
		ps.stream().forEach(p->{
			for(int i=1,I=varNames.length-2;i<=I;i++) cDisplacement(p,i);
		});
	}
	
	static void cDisplacement(Particle p,int trNum){
		int tlen=p.getTCount();
		
		float undef=Record.undef;
		
		// get crossing contour displacement, then integrate to get position for a new Particle
		float[] tr=p.getAttachedData(trNum-1+4);
		float[] gd=p.getAttachedData(trNum-1+4+varNames.length);
		
		tr=FilterModel.runningMean(tr,3,undef);
		gd=FilterModel.runningMean(gd,3,undef);
		
		float[] dsplcmt=new float[tlen-1];
		
		for(int l=1;l<tlen;l++) if(tr[l]!=undef&&tr[l-1]!=undef){
			float deltaTr=tr[l]-tr[l-1];
			
			if(deltaTr==0) dsplcmt[l-1]=0;
			else{
				if(gd[l]!=undef&&gd[l-1]!=undef) dsplcmt[l-1]=(deltaTr)/((gd[l]+gd[l-1])/2f);
				else{
					dsplcmt[l-1]=undef;
					//System.out.println("grd is undef: "+gd[l]+" "+gd[l-1]);
				}
			}
			
			if(Float.isNaN(dsplcmt[l-1])||Float.isInfinite(dsplcmt[l-1]))
			System.out.println(dsplcmt[l-1]+" "+tr[l]+" "+tr[l-1]+" "+gd[l]+" "+gd[l-1]);
			
		}else{
			dsplcmt[l-1]=undef;
			//System.out.println("tr is undef: "+tr[l]+" "+tr[l-1]);
		}
		
		for(int l=1;l<tlen;l++){
			Record o=p.getRecord(l);
			
			o.setData(trNum-1+4+varNames.length*3,dsplcmt[l-1]);
		}
	}
	
	static void cYPositionByDisplacement(List<Particle> ps,int trNum){
		ps.stream().forEach(p->toPositionParticle(p,trNum));
	}
	
	static void toPositionParticle(Particle p,int trNum){
		int tlen=p.getTCount();
		
		boolean out=false;
		
		Record init=p.getRecord(0);
		
		init.setXPos(0);
		init.setYPos(0);
		
		for(int l=1;l<tlen;l++){
			Record o=p.getRecord(l);
			
			float disp=o.getDataValue(trNum-1+4+3*varNames.length);
			
			if(disp==Record.undef){ disp=0; out=true;}
			
			o.setYPos(p.getRecord(l-1).getYPos()+disp);
			o.setXPos(0);
		}
		
		if(out) System.out.println(p.getID());
		
		p.cVelocityByPosition();
	}
	
	
	static void calDiffMap(List<Particle> ps,int trNum,String test){
		final DataDescriptor dd=DiagnosisFactory.parseContent(
			"dset ^Stat.dat\r\n" + 
			"options big_endian\r\n" + 
			"undef -9999\r\n" + 
			"title mitgcm\r\n" + 
			"xdef  140 linear 0 22000.0\r\n" + 
			"ydef  100 linear 0 22000.0\r\n" + 
			"zdef    1 linear    0 1\r\n" + 
			"tdef    1 linear 00:00Z02Jan2000 1dy\r\n" + 
			"f0   4.988021E-5\r\n" + 
			"Beta 2.151004E-11\r\n" + 
			"vars 2\r\n" + 
			"u      1 99 u\r\n" + 
			"v      1 99 v\r\n" + 
			"endvars"
		).getDataDescriptor();
		
		LagrangianStatisticsByDavis lstat=new LagrangianStatisticsByDavis(ps,dd);
		
		int tRad=16;
		int str =0;
		int end =15;
		int minT=10;
		
		float bRad=55000f;
		
		Variable[] stats1=lstat.cMaxStatisticsMapByDavisTheory     (tRad,bRad,str,end,minT);
		Variable[] stats2=lstat.cMaxStatisticsMapByTaylorTheory    (tRad,bRad,str,end,minT);
		Variable[] stats3=lstat.cMaxStatisticsMapByDispersionTheory(tRad,bRad,str,end,minT);
		
		for(int i=0;i<stats1.length;i++){
			stats1[i].setName(stats1[i].getName()+"1");
			stats2[i].setName(stats2[i].getName()+"2");
			stats3[i].setName(stats3[i].getName()+"3");
		}
		
		DataWrite dw=DataIOFactory.getDataWrite(dd,path+"dispInCC/"+test+"/LagStat/CTCoord/diffMapTr"+trNum+"InCC"+str+"_"+end+".dat");
		dw.writeData(dd,ArrayUtil.concatAll(Variable.class,stats1,stats2,stats3));
		dw.closeFile();
	}
	
	
	static void toAlongTrackFile(List<Particle> ps,String out){
		Variable[] vs=ps.stream().map(p->toTimeSeries(p)).toArray(Variable[]::new);
		
		CtlDataWriteStream cdws=new CtlDataWriteStream(out);
		cdws.writeData(vs); cdws.closeFile();
		
		StringBuilder sb=new StringBuilder();
		
		Particle p=ps.get(0); MDate tstr=new MDate(p.getTime(0));
		
		int dt=tstr.getDT(new MDate(p.getTime(1)));
		
		if(dt%86400!=0) throw new IllegalArgumentException("dt ("+dt+") cannot be divided by 86400");
		
		sb.append("dset ^"+IOUtil.getFileName(out)+"\n");
		sb.append("undef -9999\n");
		sb.append("title along-track sampling\n");
		sb.append("xdef   1 linear 0 1\n");
		sb.append("ydef   1 linear 0 1\n");
		sb.append("zdef  "+vs[0].getZCount()+" linear 0 1\n");
		sb.append("tdef "+vs[0].getTCount()+" linear "+tstr.toGradsDate()+" "+dt/86400+"dy\n");
		sb.append("vars "+vs.length+"\n");
		for(int i=0;i<vs.length;i++)
		sb.append(String.format("%-11s "+vs[i].getZCount()+" 99  %s\n","v"+(i+1),vs[i].getName()));
		sb.append("endvars\n");
		
		try(FileWriter fw=new FileWriter(IOUtil.getCompleteFileNameWithoutExtension(out)+".ctl")){ fw.write(sb.toString());}
		catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	static void toSPStatFile(List<Particle> ps,String out,int trNum){
		DataDescriptor dd=DiagnosisFactory.getDataDescriptor("H:/cartRL_advSchemes/Leith1_k0/Stat.cts");
		
		LagrangianStatisticsByDavis lstat=new LagrangianStatisticsByDavis(ps,dd);
		
		long initTime=ps.get(0).getTime(0);
		
		Predicate<Record> psudoTrackCond=r->r.getTime()==initTime;
		
		SingleParticleStatResult r1=lstat.cStatisticsByDavisTheory     (psudoTrackCond,50);
		SingleParticleStatResult r2=lstat.cStatisticsByTaylorTheory    (psudoTrackCond,50);
		SingleParticleStatResult r3=lstat.cStatisticsByDispersionTheory(psudoTrackCond,50);
		
		r1.toFile(out+"TR"+trNum+"SPInCC1.txt");
		r2.toFile(out+"TR"+trNum+"SPInCC2.txt");
		r3.toFile(out+"TR"+trNum+"SPInCC3.txt");
	}
	
	static Variable toTimeSeries(Particle p){
		int t=p.getTCount(),z=2+varNames.length*varGroups;
		
		Variable v=new Variable("v"+p.getID(),new Range(t,z,1,1));
		
		float[][][][] vdata=v.getData();
		
		for(int k=0;k<z;k++){
			float[] data=p.getAttachedData(k);
			
			for(int l=0;l<t;l++) vdata[k][0][0][l]=data[l];
		};
		
		return v;
	}
}
