//
package cartesianRL;

import java.util.List;
import java.util.stream.Stream;
import miniufo.application.contour.ContourCartesianSpatialModel;
import miniufo.basic.ArrayUtil;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.lagrangian.AttachedMeta;
import miniufo.lagrangian.LagrangianSampling;
import miniufo.lagrangian.LagrangianUtil;
import miniufo.lagrangian.Particle;
import miniufo.lagrangian.Record;
import miniufo.util.Region2D;

//
public final class SampleTrajs{
	//
	static final Region2D r1=new Region2D( 300e3f,1050e3f, 400e3f,1150e3f,"r1");
	static final Region2D r2=new Region2D(2600e3f, 500e3f,2700e3f, 600e3f,"r2");
	static final Region2D r3=new Region2D(1000e3f, 550e3f,1100e3f, 650e3f,"r3");
	
	static final AttachedMeta UVEL=new AttachedMeta("uvel",0);
	static final AttachedMeta VVEL=new AttachedMeta("vvel",1);
	static final AttachedMeta USPL=new AttachedMeta("uspl",2);
	static final AttachedMeta VSPL=new AttachedMeta("vspl",3);
	static final AttachedMeta TR9 =new AttachedMeta("tr9" ,4);
	static final AttachedMeta Yeq1=new AttachedMeta("Yeq1",5);
	static final AttachedMeta Yeq2=new AttachedMeta("Yeq2",6);
	
	
	//
	public static void main(String[] args){
		String test="ctrl";
		Stream.of(r1).forEach(region->{
			/**/
			List<Particle> ps=Parameters.getParticlesDeployedInRegion(
				Parameters.path+"fltInit_11km_All.bin",
				"H:/dispInCC/"+test+"/fltOutput/",
				region,
				UVEL,VVEL,USPL,VSPL,TR9,Yeq1,Yeq2
			);
			
			LagrangianSampling smpl=new LagrangianSampling(ps);
			
			smpl.sampleVariables("H:/dispInCC/"+test+"/Stat.cts",TR9);
			smpl.sampleVariables("H:/dispInCC/"+test+"/Keff.cts",Yeq2);
			tracerValueToEqvY(ps,"H:/dispInCC/"+test+"/Stat.cts",TR9,Yeq1);
			LagrangianUtil.writeTrajecories(ps,Parameters.path+"dispInCC/"+test+"/TXT/",true,p->true);
		});
	}
	
	
	static void tracerValueToEqvY(List<Particle> ps,String cts,AttachedMeta tracer,AttachedMeta Y1){
		DiagnosisFactory df=DiagnosisFactory.parseFile(cts);df.setPrinting(false);
		DataDescriptor dd=df.getDataDescriptor();
		
		int tLen=ps.get(0).getTCount();
		int ttag=dd.getTNum(ps.get(0).getTime(0))+1;
		
		System.out.println("ttag for initial time in tracerValueToEqvY() ("+ps.get(0).getTime(0)+"): "+ttag);
		
		ContourCartesianSpatialModel ccsm=new ContourCartesianSpatialModel(dd);
		
		for(int l=0;l<tLen;l++){
			if(l%10==0) System.out.println("tracerValueToEqvY at "+l);
			Variable tr=df.getVariables(new Range("t("+(ttag+l)+","+(ttag+l)+")",dd),tracer.name)[0];
			
			tr.replaceUndefData(Record.undef);
			
			float[][] tdata=tr.getData()[0][0];
			float[] extrema=ArrayUtil.getExtrema(tdata,tr.getUndef());
			
			double[] cVals=new double[ps.size()];
			
			for(int i=0,I=ps.size();i<I;i++){
				Particle p=ps.get(i);
				
				cVals[i]=p.getRecord(l).getDataValue(tracer);
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
			
			for(int i=0,I=ps.size();i<I;i++) ps.get(i).getRecord(l).setData(Y1,(float)Ys[i]);
		}
	}
}
