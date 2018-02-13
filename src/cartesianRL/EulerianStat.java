//
package cartesianRL;

import miniufo.basic.ArrayUtil;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.io.DataIOFactory;
import miniufo.io.DataWrite;

//
public final class EulerianStat{
	//
	private static final String path=Parameters.path;
	private static final String testPath="dispInCC/ctrl/";
	
	
	//
	public static void main(String[] args){
		DiagnosisFactory df=DiagnosisFactory.parseFile("H:/cartRL_advSchemes/Leith1_k0/Stat.cts");
		DataDescriptor dd=df.getDataDescriptor();  df.setPrinting(false);
		
		int tstr=721;
		int tend=3600;
		
		Variable um1=cMean(df,tstr,tend,"u");
		Variable vm1=cMean(df,tstr,tend,"v");
		Variable[] all1=cEKEAndEllipse(df,tstr,tend,um1,vm1);
		
		tstr=1081;
		tend=1440;
		
		Variable um2=cMean(df,tstr,tend,"u");
		Variable vm2=cMean(df,tstr,tend,"v");
		Variable[] all2=cEKEAndEllipse(df,tstr,tend,um2,vm2);
		
		um1.setName("u1"); vm1.setName("v1"); for(Variable v:all1) v.setName(v.getName()+"1");
		um2.setName("u2"); vm2.setName("v2"); for(Variable v:all2) v.setName(v.getName()+"2");
		
		um1.setUndef(0); vm1.setUndef(0); for(Variable v:all1) v.setUndef(0);
		um2.setUndef(0); vm2.setUndef(0); for(Variable v:all2) v.setUndef(0);
		
		DataWrite dw=DataIOFactory.getDataWrite(dd,path+testPath+"EulerianStat.dat");
		dw.writeData(dd,ArrayUtil.concatAll(Variable.class,new Variable[]{um1,vm1,um2,vm2},all1,all2));
		dw.closeFile();
	}
	
	static Variable[] cEKEAndEllipse(DiagnosisFactory df,int tstr,int tend,Variable um,Variable vm){
		Variable eke  =new Variable("eke"  ,um);
		Variable uvar =new Variable("uvar" ,um);
		Variable vvar =new Variable("vvar" ,um);
		Variable cvar =new Variable("cvar" ,um);
		Variable major=new Variable("major",um);
		Variable minor=new Variable("minor",um);
		Variable theta=new Variable("theta",um);
		
		  eke.setCommentAndUnit("EKE relative to the average from "+tstr+" to "+tend+" (none)");
		 uvar.setCommentAndUnit("u-variance averaged from "+tstr+" to "+tend+" (none)");
		 vvar.setCommentAndUnit("v-variance averaged from "+tstr+" to "+tend+" (none)");
		 cvar.setCommentAndUnit("uv covariance averaged from "+tstr+" to "+tend+" (none)");
		major.setCommentAndUnit("major component of variance ellipse from "+tstr+" to "+tend+" (none)");
		minor.setCommentAndUnit("minor component of variance ellipse from "+tstr+" to "+tend+" (none)");
		theta.setCommentAndUnit("angle of major component of variance ellipse from "+tstr+" to "+tend+" (none)");
		
		df.getVariablesTimeByTime(tstr,tend,um.getName(),vm.getName()).forEach(vs->{
			Variable ua=vs[0].minus(um);
			Variable va=vs[1].minus(vm);
			
			Variable ua2=ua.square();
			Variable va2=va.square();
			
			uvar.plusEq(ua2);
			vvar.plusEq(va2);
			cvar.plusEq(ua.multiply(va));
			eke.plusEq(ua2).plusEq(va2);
		});
		
		eke.divideEq((tend-tstr+1)*2f);
		uvar.divideEq(tend-tstr+1);
		vvar.divideEq(tend-tstr+1);
		cvar.divideEq(tend-tstr+1);
		
		float[][] madata=major.getData()[0][0];
		float[][] midata=minor.getData()[0][0];
		float[][] thdata=theta.getData()[0][0];
		float[][] urdata= uvar.getData()[0][0];
		float[][] vrdata= vvar.getData()[0][0];
		float[][] crdata= cvar.getData()[0][0];
		
		for(int j=0,J=um.getYCount();j<J;j++)
		for(int i=0,I=um.getXCount();i<I;i++){
			float u=urdata[j][i];
			float v=vrdata[j][i];
			float c=crdata[j][i];
			
			madata[j][i]=(float)((u+v+Math.sqrt((u-v)*(u-v)+4.0*c*c))/2.0);
			midata[j][i]=u+v-madata[j][i];
			thdata[j][i]=(float)Math.atan2(madata[j][i]-u,c);
		}
		
		return new Variable[]{uvar,vvar,cvar,major,minor,theta,eke};
	}
	
	static Variable cMean(DiagnosisFactory df,int tstr,int tend,String vname){
		DataDescriptor dd=df.getDataDescriptor();
		
		Variable mean=new Variable(vname,new Range(1,1,dd.getYCount(),dd.getXCount()));
		mean.setCommentAndUnit("average from "+tstr+" to "+tend+" of "+vname+" (none)");
		
		df.getVariableTimeByTime(tstr,tend,vname).reduce(mean,(v1,v2)->v1.plusEq(v2));
		
		mean.divideEq(tend-tstr+1);
		
		return mean;
	}
}
