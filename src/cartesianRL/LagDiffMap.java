//
package cartesianRL;

import java.util.List;
import java.util.stream.Collectors;
import fltDispersion.FLTUtils;
import fltDispersion.FltParticle;
import fltDispersion.FltRecord;
import miniufo.concurrent.ConcurrentUtil;
import miniufo.application.statisticsModel.LagrangianStatisticsByDavis;
import miniufo.basic.ArrayUtil;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Variable;
import miniufo.io.DataIOFactory;
import miniufo.io.DataWrite;
import miniufo.lagrangian.Particle;
import miniufo.lagrangian.Record;


//
public final class LagDiffMap{
	//
	private static final String path=Parameters.path;
	private static final String testPath="dispInCC/ctrl/";
	
	private static final DataDescriptor dd=DiagnosisFactory.parseContent(
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
	
	
	//
	public static void main(String[] args){
		List<Particle> ps=getParticles(path+testPath+"fltOutput/");
		
		LagrangianStat.removeMeanFlow(path+testPath+"/EulerianStat.cts",ps);
		
		ConcurrentUtil.initDefaultExecutor(2);
		calDiffMap(ps);
		ConcurrentUtil.shutdown();
	}
	
	static void calDiffMap(List<Particle> ls){
		LagrangianStatisticsByDavis lstat=new LagrangianStatisticsByDavis(ls,dd);
		
		int tRad=31;
		int str =20;
		int end =30;
		int minT=1000;
		
		float bRad=44000f;
		
		Variable[] stats1=lstat.cMeanStatisticsMapByDavisTheory     (tRad,bRad,str,end,minT);
		Variable[] stats2=lstat.cMeanStatisticsMapByTaylorTheory    (tRad,bRad,str,end,minT);
		Variable[] stats3=lstat.cMeanStatisticsMapByDispersionTheory(tRad,bRad,str,end,minT);
		
		for(int i=0;i<stats1.length;i++){
			stats1[i].setName(stats1[i].getName()+"1");
			stats2[i].setName(stats2[i].getName()+"2");
			stats3[i].setName(stats3[i].getName()+"3");
		}
		
		DataWrite dw=DataIOFactory.getDataWrite(dd,path+"dispInCC/ctrl/LagStat/XYCoord/diffMapRes"+str+"_"+end+".dat");
		dw.writeData(dd,ArrayUtil.concatAll(Variable.class,stats1,stats2,stats3));
		dw.closeFile();
	}
	
	static List<Particle> getParticles(String outPath){
		List<FltParticle> fps=FLTUtils.readFLTTrajectory(outPath,Parameters.baseTime,rec->true);
		
		return fps.stream().map(fp->toParticle(fp)).collect(Collectors.toList());
	}
	
	
	static Particle toParticle(FltParticle fp){
		Particle p=new Particle(Integer.toString(fp.id),fp.recs.size(),2,false);
		
		for(FltRecord fr:fp.recs) p.addRecord(toRecord(fr));
		
		p.setAttachedMeta(Particle.UVEL,Particle.VVEL);
		
		return p;
	}
	
	static Record toRecord(FltRecord fr){
		long time=fr.getTime();
		
		float xpos=fr.getXPos();
		float ypos=fr.getYPos();
		
		// attached data: u, v
		Record r=new Record(time,xpos,ypos,2);
		
		r.setData(Particle.UVEL,fr.getUVel()[0]);
		r.setData(Particle.VVEL,fr.getVVel()[0]);
		
		return r;
	}
}
