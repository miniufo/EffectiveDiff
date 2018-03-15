//
package cartesianRL;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Predicate;
import fltDispersion.FLTUtils;
import fltDispersion.FltParticle;
import miniufo.application.contour.ContourCartesianSpatialModel;
import miniufo.application.statisticsModel.LagrangianStatisticsByDavis;
import miniufo.application.statisticsModel.SingleParticleStatResult;
import miniufo.basic.ArrayUtil;
import miniufo.concurrent.ConcurrentUtil;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.MDate;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.io.CtlDataWriteStream;
import miniufo.io.DataIOFactory;
import miniufo.io.DataWrite;
import miniufo.io.IOUtil;
import miniufo.lagrangian.AttachedMeta;
import miniufo.lagrangian.LagrangianSampling;
import miniufo.lagrangian.LagrangianUtil;
import miniufo.lagrangian.Particle;
import miniufo.lagrangian.Record;
import miniufo.statistics.FilterModel;

//
public final class AlgTrackSample{
	//
	static final AttachedMeta UVEL=new AttachedMeta("uvel",0);
	static final AttachedMeta VVEL=new AttachedMeta("vvel",1);
	static final AttachedMeta USPL=new AttachedMeta("uspl",2);
	static final AttachedMeta VSPL=new AttachedMeta("vspl",3);
	static final AttachedMeta TRS =new AttachedMeta("tr9" ,4);
	static final AttachedMeta GRD =new AttachedMeta("grd" ,5);
	static final AttachedMeta Yeq =new AttachedMeta("Yeq" ,6);
	static final AttachedMeta DSP =new AttachedMeta("DSP" ,7);
	
	
	//
	public static void main(String[] args){
		String test="runH50";
		
		//FLTUtils.combineFLTOutIntoBin("H:/dispInCC/"+test+"/fltOutput/",Parameters.baseTime); System.exit(0);
		
		main_regions(test); System.exit(0);
		
		List<Particle> ls=new ArrayList<>(55530);
		
		try(ObjectInputStream ois=new ObjectInputStream(new BufferedInputStream(new FileInputStream("H:/dispInCC/"+test+"/fltOutput/float_trajAll.bin"),819200))){
			while(true){
				FltParticle p=(FltParticle)(ois.readObject());
				if(p==null) break;
				
				ls.add(FLTUtils.toParticle(p,false,UVEL,VVEL,USPL,VSPL,TRS,GRD,Yeq,DSP));
			}
			
		}catch(IOException|ClassNotFoundException e){ e.printStackTrace(); System.exit(0);}
		
		new LagrangianSampling(ls).sampleVariables("H:/dispInCC/"+test+"/Stat.cts",TRS,GRD);
		tracerValueToEqvY(ls,"H:/dispInCC/"+test+"/Stat.cts");
		cCrossContourDisplacement(ls);
		
		System.out.println(ls.get(1000).toString());
		
		LagrangianUtil.asRecordStream(ls).forEach(r->{
			r.setData(USPL,r.getXPos());
			r.setData(VSPL,r.getYPos());
		});
		
		ConcurrentUtil.initDefaultExecutor(2);
		cYPositionByDisplacement(ls);
		calDiffMap(ls,test);
		ConcurrentUtil.shutdown(); System.exit(0);
	}
	
	static void main_regions(String test){
		LagrangianStat.regions.forEach(region->{
			/*
			List<Particle> ps=Parameters.getParticlesDeployedInRegion(
				Parameters.path+"fltInit_11km_All.bin",
				"H:/dispInCC/"+test+"/fltOutput/",
				region,2+varNames.length*varGroups,
				ArrayUtil.concatAll(String.class,new String[]{"uvel","vvel"},uvNames,trcNames,grdNames,YeqNames,dspNames)
			);
			
			samplePTracerAndGradient(ps,"H:/dispInCC/"+test+"/Stat.cts");
			tracerValueToEqvY(ps,"H:/dispInCC/"+test+"/Stat.cts");
			cCrossContourDisplacement(ps);
			toAlongTrackFile(ps,Parameters.path+"dispInCC/"+test+"/algTrack/algTrackData"+region.getName()+".dat");
			//LagrangianUtil.writeTrajecories(ps,path+"dispInCC/"+test+"/TXT/",true,p->true);
			
			LagrangianUtil.writeAsBinaryFile(ps,Parameters.path+"dispInCC/"+test+"/algTrack/Particles"+region.getName()+".bin");*/
			
			List<Particle> ps=LagrangianUtil.readFromBinaryFile(Parameters.path+"dispInCC/"+test+"/algTrack/Particles"+region.getName()+".bin");
			
			cYPositionByDisplacement(ps);
			toSPStatFile(ps,Parameters.path+"dispInCC/"+test+"/LagStat/CTCoord/"+region.getName());
		});
	}
	
	
	static void tracerValueToEqvY(List<Particle> ps,String cts){
		DiagnosisFactory df=DiagnosisFactory.parseFile(cts);df.setPrinting(false);
		DataDescriptor dd=df.getDataDescriptor();
		
		int tLen=ps.get(0).getTCount();
		int ttag=dd.getTNum(ps.get(0).getTime(0))+1;
		
		System.out.println("ttag for initial time in tracerValueToEqvY() ("+ps.get(0).getTime(0)+"): "+ttag);
		
		ContourCartesianSpatialModel ccsm=new ContourCartesianSpatialModel(dd);
		//Variable[] trs=df.getVariables(new Range("t("+(ttag+1)+","+(ttag+1)+")",dd),varNames);
		
		for(int l=0;l<tLen;l++){
			if(l%10==0) System.out.println("tracerToEqvY at "+l);
			
			Variable tr=df.getVariables(new Range("t("+(ttag+l)+","+(ttag+l)+")",dd),TRS.name)[0];
			
			tr.replaceUndefData(Record.undef);
			
			float[][] tdata=tr.getData()[0][0];
			float[] extrema=ArrayUtil.getExtrema(tdata,tr.getUndef());
			
			double[] cVals=new double[ps.size()];
			
			for(int i=0,I=ps.size();i<I;i++){
				Particle p=ps.get(i);
				
				cVals[i]=p.getRecord(l).getDataValue(TRS);
				if(cVals[i]!=Record.undef&&cVals[i]<extrema[0]){
					System.out.println(p.getID()+", i: "+i+",   cVals[i]: "+cVals[i]+",  min: "+extrema[0]+", "+tr.getName());
					cVals[i]=extrema[0];
				}
				if(cVals[i]!=Record.undef&&cVals[i]>extrema[1]){
					System.out.println(p.getID()+", i: "+i+",   cVals[i]: "+cVals[i]+",  max: "+extrema[1]+", "+tr.getName());
					cVals[i]=extrema[1];
				}
			}
			
			double[] Ys=ccsm.computeEquivalentYs(tr,cVals,2,false);
			
			for(int i=0,I=ps.size();i<I;i++) ps.get(i).getRecord(l).setData(Yeq,(float)Ys[i]);
		}
	}
	
	static void cCrossContourDisplacement(List<Particle> ps){
		System.out.println("cross-contour displacement");
		ps.stream().forEach(p->cDisplacement(p));
	}
	
	static void cDisplacement(Particle p){
		int tlen=p.getTCount();
		
		float undef=Record.undef;
		
		// get crossing contour displacement, then integrate to get position for a new Particle
		float[] tr=p.getAttachedData(TRS);
		float[] gd=p.getAttachedData(GRD);
		
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
			
			o.setData(DSP,dsplcmt[l-1]);
		}
	}
	
	static void cYPositionByDisplacement(List<Particle> ps){
		ps.stream().forEach(p->toPositionParticle(p));
	}
	
	static void toPositionParticle(Particle p){
		int tlen=p.getTCount();
		
		boolean out=false;
		
		Record init=p.getRecord(0);
		
		init.setXPos(0);
		init.setYPos(0);
		
		for(int l=1;l<tlen;l++){
			Record o=p.getRecord(l);
			
			float disp=o.getDataValue(DSP);
			
			if(disp==Record.undef){ disp=0; out=true;}
			
			o.setYPos(p.getRecord(l-1).getYPos()+disp);
			o.setXPos(0);
		}
		
		if(out) System.out.println(p.getID());
		
		p.cVelocityByPosition();
	}
	
	
	static void calDiffMap(List<Particle> ps,String test){
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
		
		DataWrite dw=DataIOFactory.getDataWrite(dd,Parameters.path+"dispInCC/"+test+"/LagStat/CTCoord/diffMap"+TRS.name+"InCC"+str+"_"+end+".dat");
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
	
	static void toSPStatFile(List<Particle> ps,String out){
		DataDescriptor dd=DiagnosisFactory.getDataDescriptor("H:/cartRL_advSchemes/Leith1_k0/Stat.cts");
		
		LagrangianStatisticsByDavis lstat=new LagrangianStatisticsByDavis(ps,dd);
		
		long initTime=ps.get(0).getTime(0);
		
		Predicate<Record> psudoTrackCond=r->r.getTime()==initTime;
		
		SingleParticleStatResult r1=lstat.cStatisticsByDavisTheory     (psudoTrackCond,50);
		SingleParticleStatResult r2=lstat.cStatisticsByTaylorTheory    (psudoTrackCond,50);
		SingleParticleStatResult r3=lstat.cStatisticsByDispersionTheory(psudoTrackCond,50);
		
		r1.toFile(out+TRS.name+"SPInCC1.txt");
		r2.toFile(out+TRS.name+"SPInCC2.txt");
		r3.toFile(out+TRS.name+"SPInCC3.txt");
	}
	
	static Variable toTimeSeries(Particle p){
		AttachedMeta[] meta=p.getAttachedMeta();
		
		int t=p.getTCount(),z=meta.length;
		
		Variable v=new Variable("v"+p.getID(),new Range(t,z,1,1));
		
		float[][][][] vdata=v.getData();
		
		for(int k=0;k<z;k++){
			float[] data=p.getAttachedData(meta[k]);
			
			for(int l=0;l<t;l++) vdata[k][0][0][l]=data[l];
		};
		
		return v;
	}
}
