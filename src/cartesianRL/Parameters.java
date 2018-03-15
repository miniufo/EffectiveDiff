//
package cartesianRL;

import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import common.MITgcmUtil.DataPrec;
import fltDispersion.FLTUtils;
import fltDispersion.FltInitData;
import fltDispersion.FltParticle;
import miniufo.diagnosis.MDate;
import miniufo.lagrangian.AttachedMeta;
import miniufo.lagrangian.Particle;
import miniufo.util.Region2D;

//
public final class Parameters{
	//
	public static final int y=400;
	public static final int x=560;
	
	public static final float res=5500;
	
	public static final float depth=3500;
	public static final float tauMax=0.5f;
	
	public static final MDate baseTime=new MDate(2000,1,1,0,0);
	
	public static final String path="D:/Data/MITgcm/barotropicDG/BetaCartRL/";
	
	public static final Region2D region=new Region2D(0,0,res*(x-1),res*(y-1));
	
	
	public static List<Particle> getParticlesDeployedInRegion(String initFile,String outPath,Region2D r){
		return getParticlesDeployedInRegion(initFile,outPath,r,Particle.UVEL,Particle.VVEL);
	}
	
	public static List<Particle> getParticlesDeployedInRegion(String initFile,String outPath,Region2D r,AttachedMeta... meta){
		Predicate<FltInitData> inRegion=fid->{
			return	fid.xpart>r.getXMin()&&fid.xpart<r.getXMax()&&
					fid.ypart>r.getYMin()&&fid.ypart<r.getYMax();
		};
		
		List<FltInitData> inis=FLTUtils.readFLTInitFile(initFile,DataPrec.float32,inRegion);
		
		float[] ids=new float[inis.size()];
		
		for(int i=0,I=ids.length;i<I;i++) ids[i]=inis.get(i).npart;
		
		List<FltParticle> fps=FLTUtils.readFLTTrajectory(outPath,baseTime,rec->IDInList(ids,rec.getID()));
		
		return fps.stream().map(fp->FLTUtils.toParticle(fp,false,meta)).collect(Collectors.toList());
	}
	
	
	static boolean IDInList(float[] ids,float id){
		for(float oneID:ids) if(oneID==id) return true;
		return false;
	}
}
