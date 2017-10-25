//
package cartesianRL;

import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import fltDispersion.FLTUtils;
import fltDispersion.FltParticle;
import fltDispersion.FltRecord;
import miniufo.application.basic.DynamicMethodsInCTS;
import miniufo.application.contour.ContourCartesianSpatialModel;
import miniufo.application.statisticsModel.LagrangianStatisticsByDavis;
import miniufo.application.statisticsModel.SingleParticleStatResult;
import miniufo.basic.ArrayUtil;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.CartesianSpatialModel;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.MDate;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.io.CtlDataWriteStream;
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
	private static int varNumInCts=8;
	
	private static String[] varNames=new String[]{
		// == varNumInCts, sampled offline
		"u","v","tr1","tr2","tr3","tr4","tr5","tr6"
	};
	
	private static final String path=Grids.path;
	
	
	//
	public static void main(String[] args){
		String test="runH150";
		List<Particle> ps=readFromFltFile(path+"dispInCC/"+test+"/");
		
		samplePTracers(ps,path+"dispInCC/"+test+"/Stat.cts");
		
		toAlongTrackFile(ps,path+"dispInCC/"+test+"/algTrackPTracer.dat");
		LagrangianUtil.writeTrajecories(ps,path+"dispInCC/"+test+"/TXT/",true,p->true);
		
		tracerValueToEqvY(ps,path+"dispInCC/"+test+"/Stat.cts");
		toAlongTrackFile(ps,path+"dispInCC/"+test+"/algTrackEquvYs.dat"); System.exit(0);
		
		List<Particle> patch1=ps.stream().filter(p->p.getID().startsWith("1")).map(p->toPositionParticle(p,6)).collect(Collectors.toList());
		List<Particle> patch2=ps.stream().filter(p->p.getID().startsWith("2")).map(p->toPositionParticle(p,6)).collect(Collectors.toList());
		List<Particle> patch3=ps.stream().filter(p->p.getID().startsWith("3")).map(p->toPositionParticle(p,6)).collect(Collectors.toList());
		List<Particle> patch4=ps.stream().filter(p->p.getID().startsWith("4")).map(p->toPositionParticle(p,6)).collect(Collectors.toList());
		List<Particle> patch5=ps.stream().filter(p->p.getID().startsWith("5")).map(p->toPositionParticle(p,6)).collect(Collectors.toList());
		List<Particle> patch6=ps.stream().filter(p->p.getID().startsWith("6")).map(p->toPositionParticle(p,6)).collect(Collectors.toList());
		
		toSPStatFile(patch1,path+"dispInCC/"+test+"/",1);
		toSPStatFile(patch2,path+"dispInCC/"+test+"/",2);
		toSPStatFile(patch3,path+"dispInCC/"+test+"/",3);
		toSPStatFile(patch4,path+"dispInCC/"+test+"/",4);
		toSPStatFile(patch5,path+"dispInCC/"+test+"/",5);
		toSPStatFile(patch6,path+"dispInCC/"+test+"/",6);
	}
	
	static List<Particle> readFromFltFile(String folderName){
		List<FltParticle> fps=FLTUtils.readFLTTrajectory(folderName,new MDate(2000,1,1,0,0));
		
		return fps.stream().map(fp->toParticle(fp)).collect(Collectors.toList());
	}
	
	static void toSPStatFile(List<Particle> ps,String out,int trNum){
		DataDescriptor dd=DiagnosisFactory.DFHalf.getDataDescriptor();
		
		LagrangianStatisticsByDavis lstat=new LagrangianStatisticsByDavis(ps,dd);
		
		Predicate<Record> psudoTrackCond=r->true;
		
		SingleParticleStatResult r1=lstat.cStatisticsByDavisTheory1(psudoTrackCond,80);
		SingleParticleStatResult r2=lstat.cStatisticsByDavisTheory2(psudoTrackCond,80);
		SingleParticleStatResult r3=lstat.cStatisticsByDavisTheory3(psudoTrackCond,80);
		
		r1.toFile(out+"patch"+trNum+"SPInCC1.txt");
		r2.toFile(out+"patch"+trNum+"SPInCC2.txt");
		r3.toFile(out+"patch"+trNum+"SPInCC3.txt");
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
	
	static void tracerValueToEqvY(List<Particle> ps,String ctl){
		DiagnosisFactory df=DiagnosisFactory.parseFile(ctl);df.setPrinting(false);
		DataDescriptor dd=df.getDataDescriptor();
		
		int tLen=ps.get(0).getTCount();
		int ttag=dd.getTNum(ps.get(0).getTime(0));
		System.out.println("ttag: "+ttag);
		
		ContourCartesianSpatialModel ccsm=new ContourCartesianSpatialModel(dd);
		Variable[] trs=df.getVariables(new Range("t("+(ttag+1)+","+(ttag+1)+")",dd),"tr1","tr2","tr3","tr4","tr5","tr6");
		
		for(int l=0;l<tLen;l++){
			//Variable[] trs=df.getVariables(new Range("t("+(ttag+l+1)+","+(ttag+l+1)+")",dd),"tr1","tr2","tr3","tr4","tr5","tr6");
			
			int offset=2+varNumInCts-trs.length;
			
			for(int k=0,K=trs.length;k<K;k++){
				trs[k].replaceUndefData(Record.undef);
				
				float[][] tdata=trs[k].getData()[0][0];
				float[] extrema=ArrayUtil.getExtrema(tdata,trs[k].getUndef());
				
				double[] cVals=new double[ps.size()];
				
				for(int i=0,I=ps.size();i<I;i++){
					Particle p=ps.get(i);
					
					cVals[i]=p.getRecord(l).getDataValue(k+offset);
					if(cVals[i]!=Record.undef&&cVals[i]<extrema[0]){
						System.out.println(p.getID()+", z: "+k+", l: "+i+",   cVals[l]: "+cVals[i]+",  min: "+extrema[0]);
						cVals[i]=extrema[0];
					}
					if(cVals[i]!=Record.undef&&cVals[i]>extrema[1]){
						System.out.println(p.getID()+", z: "+k+", l: "+i+",   cVals[l]: "+cVals[i]+",  max: "+extrema[1]);
						cVals[i]=extrema[1];
					}
				}
				
				double[] Ys=ccsm.computeEquivalentYs(trs[k],cVals,2);
				
				for(int i=0,I=ps.size();i<I;i++) ps.get(i).getRecord(l).setData(k+offset,(float)Ys[i]);
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
	
	static void samplePTracers(List<Particle> ps,String cts){
		DiagnosisFactory df=DiagnosisFactory.parseFile(cts);
		DataDescriptor dd=df.getDataDescriptor();
		
		CartesianSpatialModel csm=new CartesianSpatialModel(dd);
		GridDataFetcher gdf=new GridDataFetcher(dd);
		
		DynamicMethodsInCTS dm=new DynamicMethodsInCTS(csm);
		
		System.out.println("ttag for initial time ("+ps.get(0).getTime(0)+"): "+dd.getTNum(ps.get(0).getTime(0)));
		
		for(int m=0;m<varNumInCts;m++) if(varNames[m].equalsIgnoreCase("tr6")){
			TicToc.tic("samping "+varNames[m]);
			Variable buffer=gdf.prepareXYTBuffer(varNames[m],1);buffer.replaceUndefData(Record.undef);
			Variable bufGrd=dm.c2DGradientMagnitude(buffer);
			
			final int mm=m+2;
			ps.stream().flatMap(p->p.stream()).forEach(r->{
				float tr=gdf.fetchXYTBuffer(r.getXPos(),r.getYPos(),r.getTime(),buffer);
				float gd=gdf.fetchXYTBuffer(r.getXPos(),r.getYPos(),r.getTime(),bufGrd);
				
				r.setData(mm            ,tr);
				r.setData(mm+varNumInCts,gd);
			});
			TicToc.toc(TimeUnit.SECONDS);
		}
		
		gdf.closeFile();
	}
	
	static Particle toPositionParticle(Particle p,int trNum){
		int trIdx=3+trNum;	// tr6 index in attachedData
		int tlen=p.getTCount();
		
		float undef=Record.undef;
		
		// get crossing contour displacement, then integrate to get position for a new Particle
		float[] tr=p.getAttachedData(trIdx  );
		float[] gd=p.getAttachedData(trIdx+8);
		
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
		
		int attLen=p.getRecord(0).getDataLength();
		
		Particle r=new Particle(p.getID(),tlen,attLen);
		
		r.addRecord(new Record(p.getRecord(0).getTime(),0,0,attLen));
		
		boolean out=false;
		
		for(int l=1;l<tlen;l++){
			Record o=p.getRecord(l);
			Record n=new Record(o.getTime(),0,0,attLen);
			
			float amount=dsplcmt[l-1];
			if(amount==undef){ amount=0; out=true;}
			n.setYPos(r.getRecord(l-1).getYPos()+amount);
			//n.setYPos(amount);
			
			//n.setData(0,o.getDataValue(0));
			//n.setData(1,o.getDataValue(1));
			n.setData(trIdx,o.getDataValue(trIdx));
			
			r.addRecord(n);
		}
		
		if(out) System.out.println(p.getID());
		
		r.cVelocityByPosition();
		
		return r;
	}
	
	static Variable toTimeSeries(Particle p){
		int t=p.getTCount(),z=2+varNumInCts*2;
		
		Variable v=new Variable("v"+p.getID(),new Range(t,z,1,1));
		
		float[][][][] vdata=v.getData();
		
		for(int k=0;k<z;k++){
			float[] data=p.getAttachedData(k);
			
			for(int l=0;l<t;l++) vdata[k][0][0][l]=data[l];
		};
		
		return v;
	}
	
	static Particle toParticle(FltParticle fp){
		Particle p=new Particle(Integer.toString(fp.id),fp.recs.size(),2+varNumInCts*2);
		
		for(FltRecord fr:fp.recs) p.addRecord(toRecord(fr));
		
		p.setAttachedDataNames(ArrayUtil.concatAll(String.class,new String[]{"uflt","vflt"},varNames));
		
		return p;
	}
	
	static Record toRecord(FltRecord fr){
		long time=fr.getTime();
		
		float xpos=fr.getXPos();
		float ypos=fr.getYPos();
		
		// attached data: uflt, vflt, u, v, tr1, ..., trN, ugrd, vgrd, tr1grd, ..., trNgrd
		Record r=new Record(time,xpos,ypos,2+varNumInCts*2);
		
		r.setData(0,fr.getUVel()[0]);
		r.setData(1,fr.getVVel()[0]);
		
		return r;
	}
}
