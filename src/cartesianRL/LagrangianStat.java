//
package cartesianRL;

import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.function.Predicate;
import java.util.stream.Stream;
import miniufo.application.statisticsModel.SingleParticleStatResult;
import miniufo.application.statisticsModel.TwoParticleStatResult;
import miniufo.application.statisticsModel.TwoParticleStatistics;
import miniufo.basic.InterpolationModel;
import miniufo.application.statisticsModel.LagrangianStatisticsByDavis;
import miniufo.application.statisticsModel.ParticlePair;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Variable;
import miniufo.lagrangian.LagrangianSampling;
import miniufo.lagrangian.Particle;
import miniufo.lagrangian.Record;
import miniufo.util.GridDataFetcher;
import miniufo.util.Region2D;
import miniufo.util.TicToc;

//
public final class LagrangianStat{
	//
	private static final String path=Parameters.path;
	private static final String testPath="dispInCC/ctrl/";
	
	public static final Region2D r1 =new Region2D( 750e3f,2000e3f, 900e3f,2150e3f,"r1");
	public static final Region2D r2 =new Region2D(2000e3f,2000e3f,2150e3f,2150e3f,"r2");
	
	public static final Region2D r3 =new Region2D(  40e3f,1550e3f, 190e3f,1700e3f,"r3");
	public static final Region2D r4 =new Region2D( 900e3f,1550e3f,1050e3f,1700e3f,"r4");
	public static final Region2D r5 =new Region2D(1900e3f,1550e3f,2050e3f,1700e3f,"r5");
	public static final Region2D r6 =new Region2D(2900e3f,1550e3f,3050e3f,1700e3f,"r6");
	
	public static final Region2D r7 =new Region2D( 500e3f,1030e3f, 650e3f,1180e3f,"r7");
	public static final Region2D r8 =new Region2D(1450e3f,1030e3f,1600e3f,1180e3f,"r8");
	public static final Region2D r9 =new Region2D(2500e3f,1030e3f,2650e3f,1180e3f,"r9");
	
	public static final Region2D r10=new Region2D(  40e3f, 470e3f, 190e3f, 620e3f,"r10");
	public static final Region2D r11=new Region2D( 900e3f, 470e3f,1050e3f, 620e3f,"r11");
	public static final Region2D r12=new Region2D(1900e3f, 470e3f,2050e3f, 620e3f,"r12");
	public static final Region2D r13=new Region2D(2900e3f, 470e3f,3050e3f, 620e3f,"r13");
	
	public static final Region2D r14=new Region2D( 750e3f,  50e3f, 900e3f, 200e3f,"r14");
	public static final Region2D r15=new Region2D(2000e3f,  50e3f,2150e3f, 200e3f,"r15");
	
	public static final Stream<Region2D> regions=Stream.of(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15);
	
	
	//
	public static void main(String[] args){
		regions.forEach(region->{
		//Stream.of(r1).forEach(region->{
			System.out.println("computing "+region.getName()+" width: "+(region.getXMax()-region.getXMin())+" height: "+(region.getYMax()-region.getYMin()));
			
			List<Particle> ps=Parameters.getParticlesDeployedInRegion(path+"fltInit_11km_All.bin",path+testPath+"fltOutput/",region);
			
			new LagrangianSampling(ps).sampleVariables("H:/cartRL_advSchemes/Leith1_k0/Stat.cts");
			
			//LagrangianUtil.writeTrajecories(ps,path+testPath+"TXT/TXT"+region.getName()+"/",true,r->true);
			
			//removeMeanFlow("H:/cartRL_advSchemes/Leith1_k0/meanFlow.cts",ps);
			
			//LagrangianUtil.writeTrajecories(ps,path+testPath+"TXT/TXT"+region.getName()+"NoMean/",true,r->true);
			
			long iniTime=ps.get(0).getTime(0);
			
			singleParticleStat(ps,r->r.getTime()==iniTime,region.getName());
			
			System.out.println();
		});
	}
	
	
	static void removeMeanFlow(String meanFile,List<Particle> ps){
		TicToc.tic("removing the mean flow");
		
		DiagnosisFactory df=DiagnosisFactory.parseFile(meanFile);
		DataDescriptor dd=df.getDataDescriptor();
		
		GridDataFetcher gdsU=new GridDataFetcher(dd);
		GridDataFetcher gdsV=new GridDataFetcher(dd);
		
		Variable bufU=gdsU.prepareXYBuffer("u1",1,1);
		Variable bufV=gdsV.prepareXYBuffer("v1",1,1);
		
		float undef=bufU.getUndef();
		float disX=dd.getDXDef()[0];
		float disY=dd.getDYDef()[0];
		float xmin=dd.getXDef().getMin();
		float xmax=dd.getXDef().getMax()-disX;
		float ymin=dd.getYDef().getMin();
		float ymax=dd.getYDef().getMax()-disY;
		
		for(Particle p:ps)
		for(int ll=0,LL=p.getTCount();ll<LL;ll++){
			Record r=p.getRecord(ll);
			
			float xpos=r.getXPos();
			float ypos=r.getYPos();
			
			float urec=gdsU.fetchXYBuffer(xpos,ypos,bufU);
			float vrec=gdsV.fetchXYBuffer(xpos,ypos,bufV);
			
			float ucurr=r.getDataValue(Particle.UVEL);
			float vcurr=r.getDataValue(Particle.VVEL);
			
			if(urec==undef){
				if(xpos<xmin&&ypos>ymin){
					float uinner=gdsU.fetchXYBuffer(xmin,ypos,bufU);
					urec=InterpolationModel.linearInterpolation(0,uinner,xpos/xmin);
				}
				if(ypos<ymin&&xpos>xmin){
					float uinner=gdsU.fetchXYBuffer(xpos,ymin,bufU);
					urec=InterpolationModel.linearInterpolation(0,uinner,ypos/ymin);
				}
				if(ypos<ymin&&xpos<xmin){
					float uinner=gdsU.fetchXYBuffer(xmin,ymin,bufU);
					urec=InterpolationModel.bilinearInterpolation(0,0,0,uinner,xpos/xmin,ypos/ymin);
				}
				if(xpos>xmax&&ypos<ymax){
					float uinner=gdsU.fetchXYBuffer(xmax,ypos,bufU);
					urec=InterpolationModel.linearInterpolation(uinner,0,(xpos-xmax)/disX);
				}
				if(xpos<xmax&&ypos>ymax){
					float uinner=gdsU.fetchXYBuffer(xpos,ymax,bufU);
					urec=InterpolationModel.linearInterpolation(uinner,0,(ypos-ymax)/disY);
				}
				if(urec==undef)
				System.out.println("at "+r.getTime()+" um ("+urec+") undef for xpos:"+xpos+" ("+(xpos>xmax)+"), ypos:"+ypos+" ("+(ypos>ymax)+")"+xmax+" "+ymax+" "+gdsU.fetchXYBuffer(xpos,ymax,bufU)+" "+(ypos-ymax)/disY);
			}
			
			if(vrec==undef){
				if(xpos<xmin&&ypos>ymin){
					float vinner=gdsV.fetchXYBuffer(xmin,ypos,bufV);
					vrec=InterpolationModel.linearInterpolation(0,vinner,xpos/xmin);
				}
				if(ypos<ymin&&xpos>xmin){
					float vinner=gdsV.fetchXYBuffer(xpos,ymin,bufV);
					vrec=InterpolationModel.linearInterpolation(0,vinner,ypos/ymin);
				}
				if(ypos<ymin&&xpos<xmin){
					float vinner=gdsV.fetchXYBuffer(xmin,ymin,bufV);
					vrec=InterpolationModel.bilinearInterpolation(0,0,0,vinner,xpos/xmin,ypos/ymin);
				}
				if(xpos>xmax&&ypos<ymax){
					float vinner=gdsV.fetchXYBuffer(xmax,ypos,bufV);
					vrec=InterpolationModel.linearInterpolation(vinner,0,(xpos-xmax)/disX);
				}
				if(xpos<xmax&&ypos>ymax){
					float vinner=gdsV.fetchXYBuffer(xpos,ymax,bufV);
					vrec=InterpolationModel.linearInterpolation(vinner,0,(ypos-ymax)/disY);
				}
				if(vrec==undef)
				System.out.println("at "+r.getTime()+" vm ("+urec+") undef for xpos:"+xpos+" ("+(xpos>xmax)+"), ypos:"+ypos+" ("+(ypos>ymax)+")");
			}
			
			r.setData(Particle.UVEL,ucurr-urec);
			r.setData(Particle.VVEL,vcurr-vrec);
		}
		
		TicToc.toc(TimeUnit.SECONDS);
	}
	
	
	static void multiParticleStat(List<Particle> ps){
		DataDescriptor dd=DiagnosisFactory.DFHalf.getDataDescriptor();
		
		TwoParticleStatistics tps=new TwoParticleStatistics(ps,dd);
		
		List<ParticlePair> pairs=tps.selectPairs(pair->pair.cDistance(0)<=6000);
		
		//for(ParticlePair pp:pairs) System.out.println(pp);
		
		TwoParticleStatResult result=tps.cRelativeDispersion(pairs);
		
		result.toFile(path+"dispInCC/TestRnd0/TP1.txt");
	}
	
	static void singleParticleStat(List<Particle> ps,Predicate<Record> psudoTrack,String outPrefix){
		DataDescriptor dd=DiagnosisFactory.getDataDescriptor("H:/cartRL_advSchemes/Leith1_k0/Stat.cts");
		
		LagrangianStatisticsByDavis lstat=new LagrangianStatisticsByDavis(ps,dd);
		
		SingleParticleStatResult r1=lstat.cStatisticsByDavisTheory     (psudoTrack, 60); r1.toFile(path+testPath+"LagStat/"+outPrefix+"SP1.txt");
		SingleParticleStatResult r2=lstat.cStatisticsByTaylorTheory    (psudoTrack, 60); r2.toFile(path+testPath+"LagStat/"+outPrefix+"SP2.txt");
		SingleParticleStatResult r3=lstat.cStatisticsByDispersionTheory(psudoTrack, 60); r3.toFile(path+testPath+"LagStat/"+outPrefix+"SP3.txt");
	}
}
