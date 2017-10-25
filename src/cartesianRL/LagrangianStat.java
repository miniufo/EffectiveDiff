//
package cartesianRL;

import java.util.List;
import java.util.stream.Collectors;
import fltDispersion.FLTUtils;
import fltDispersion.FltParticle;
import fltDispersion.FltRecord;
import miniufo.application.statisticsModel.SingleParticleStatResult;
import miniufo.application.statisticsModel.TwoParticleStatResult;
import miniufo.application.statisticsModel.TwoParticleStatistics;
import miniufo.application.statisticsModel.LagrangianStatisticsByDavis;
import miniufo.application.statisticsModel.ParticlePair;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.MDate;
import miniufo.lagrangian.Particle;
import miniufo.lagrangian.Record;

//
public final class LagrangianStat{
	//
	private static final String path=Grids.path;
	
	
	//
	public static void main(String[] args){
		List<FltParticle> fps=FLTUtils.readFLTTrajectory(path+"dispInCC/ctrl/",new MDate(2000,1,1,0,0));
		FLTUtils.writeTrajAndGS(fps,path+"dispInCC/ctrl/TXT/",Grids.region); System.exit(0);
		
		List<Particle> ps=fps.stream().map(fp->toParticle(fp)).collect(Collectors.toList());
		System.out.println(ps.size());
		
		//singleParticleStat(ps);
		multiParticleStat(ps);
	}
	
	
	static void multiParticleStat(List<Particle> ps){
		DataDescriptor dd=DiagnosisFactory.DFHalf.getDataDescriptor();
		
		TwoParticleStatistics tps=new TwoParticleStatistics(ps,dd);
		
		List<ParticlePair> pairs=tps.selectPairs(pair->pair.cDistance(0)<=6000);
		
		//for(ParticlePair pp:pairs) System.out.println(pp);
		
		TwoParticleStatResult result=tps.cRelativeDispersion(pairs);
		
		result.toFile(path+"dispInCC/TestRnd0/TP1.txt");
	}
	
	static void singleParticleStat(List<Particle> ps){
		DataDescriptor dd=DiagnosisFactory.DFHalf.getDataDescriptor();
		
		LagrangianStatisticsByDavis lstat=new LagrangianStatisticsByDavis(ps,dd);
		
		SingleParticleStatResult r1=lstat.cStatisticsByDavisTheory1(r->true, 60);
		SingleParticleStatResult r2=lstat.cStatisticsByDavisTheory2(r->true, 60);
		SingleParticleStatResult r3=lstat.cStatisticsByDavisTheory3(r->true, 60);
		
		r1.toFile(path+"dispInCC/TestRnd0/SP1.txt");
		r2.toFile(path+"dispInCC/TestRnd0/SP2.txt");
		r3.toFile(path+"dispInCC/TestRnd0/SP3.txt");
	}
	
	static Particle toParticle(FltParticle fp){
		Particle p=new Particle(Integer.toString(fp.id),fp.recs.size(),2);
		
		for(FltRecord fr:fp.recs) p.addRecord(toRecord(fr));
		
		p.setAttachedDataNames("uvel","vvel");
		
		return p;
	}
	
	static Record toRecord(FltRecord fr){
		long time=fr.getTime();
		
		float xpos=fr.getXPos();
		float ypos=fr.getYPos();
		
		// attached data: u, v
		Record r=new Record(time,xpos,ypos,2);
		
		r.setData(0,fr.getUVel()[0]);
		r.setData(1,fr.getVVel()[0]);
		
		return r;
	}
}
